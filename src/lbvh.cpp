#include "lbvh.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "math_utils.h"
#include "types.h"

namespace LBVH {

void compute_morton_and_sort_indices(size_t num_particles,
                                     const Vec3& domain_size,
                                     const std::vector<double>& pos_x,
                                     const std::vector<double>& pos_y,
                                     const std::vector<double>& pos_z,
                                     std::vector<uint64_t>& morton_codes,
                                     std::vector<int>& sorted_indices) {
    if (num_particles == 0) return;

    morton_codes.resize(num_particles);
    sorted_indices.resize(num_particles);

    double inv_domain_x = 1.0 / domain_size.x;
    double inv_domain_y = 1.0 / domain_size.y;
    double inv_domain_z = 1.0 / domain_size.z;
    // We use 21 bits per dimension (2^21 = 2097152) to fit in a 64-bit int
    double bound = 2097152.0;

// Compute Morton codes for all particles
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        // Normalize coordinates to [0, 1) and scale to the 21-bit integer range
        uint32_t x = static_cast<uint32_t>(
            fmod(pos_x[i] * inv_domain_x + 1.0, 1.0) * bound);
        uint32_t y = static_cast<uint32_t>(
            fmod(pos_y[i] * inv_domain_y + 1.0, 1.0) * bound);
        uint32_t z = static_cast<uint32_t>(
            fmod(pos_z[i] * inv_domain_z + 1.0, 1.0) * bound);

        morton_codes[i] = morton3D(x, y, z);
    }

    // Initialize indices and sort them based on the Morton codes
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&](int a, int b) { return morton_codes[a] < morton_codes[b]; });

    // Ensure the morton_codes array is also sorted so it is ready for the
    // Karras tree topology algorithm without needing extra shuffles later.
    std::vector<uint64_t> sorted_morton(num_particles);
    for (size_t i = 0; i < num_particles; ++i) {
        sorted_morton[i] = morton_codes[sorted_indices[i]];
    }
    morton_codes = std::move(sorted_morton);
}

void build_topology_and_aggregate(size_t num_particles,
                                  const std::vector<double>& sorted_pos_x,
                                  const std::vector<double>& sorted_pos_y,
                                  const std::vector<double>& sorted_pos_z,
                                  const std::vector<double>& sorted_mass,
                                  const std::vector<double>* sorted_h,
                                  std::vector<uint64_t>& morton_codes,
                                  std::vector<BVHNode>& bvh_nodes) {
    if (num_particles == 0) return;

    bvh_nodes.resize(2 * num_particles - 1);
    const int n = static_cast<int>(num_particles);

    // Utility lambda to find the longest common prefix between two Morton codes
    auto delta = [&](int i, int j) -> int {
        if (j < 0 || j >= n) return -1;
        uint64_t code_i = morton_codes[i];
        uint64_t code_j = morton_codes[j];

        if (code_i == code_j) {
            // Tie-breaker for identical coordinates using the original index.
            // __builtin_clzll(unsigned long long) returns the number of leading
            // 0-bits starting at the most significant bit position.
            // If arg == 0 the result is undefined
            return 64 + __builtin_clzll(static_cast<unsigned long long>(i ^ j));
        }
        return __builtin_clzll(code_i ^ code_j);
    };

// Initialize Leaf Nodes
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        int leaf_idx = n - 1 + i;
        bvh_nodes[leaf_idx].particle_idx = i;
        bvh_nodes[leaf_idx].left_child = -1;
        bvh_nodes[leaf_idx].right_child = -1;
        // The parent will be set by the internal node that points to this leaf
    }

// Construct Internal Nodes in parallel (Karras 2012)
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n - 1; ++i) {
        // Determine direction of the range (+1 or -1)
        int d = (delta(i, i + 1) - delta(i, i - 1)) > 0 ? 1 : -1;

        // Compute upper bound for the length of the range
        int delta_min = delta(i, i - d);
        int l_max = 2;
        while (delta(i, i + l_max * d) > delta_min) {
            l_max *= 2;
        }

        // Find the other end of the range using binary search
        int l = 0;
        for (int t = l_max / 2; t >= 1; t /= 2) {
            if (delta(i, i + (l + t) * d) > delta_min) {
                l += t;
            }
        }
        int j = i + l * d;

        // Find the split position using binary search
        int delta_node = delta(i, j);
        int s = 0;
        int t = l;
        do {
            t = (t + 1) >> 1;  // ceil(t/2)
            if (s + t < l && delta(i, i + (s + t) * d) > delta_node) {
                s += t;
            }
        } while (t > 1);

        int split = i + s * d + std::min(d, 0);
        int min_idx = std::min(i, j);
        int max_idx = std::max(i, j);

        // Assign children
        int left_child, right_child;

        if (min_idx == split) {
            left_child = n - 1 + split;  // Points to leaf
        } else {
            left_child = split;  // Points to internal node
        }

        if (max_idx == split + 1) {
            right_child = n - 1 + split + 1;  // Points to leaf
        } else {
            right_child = split + 1;  // Points to internal node
        }

        bvh_nodes[i].left_child = left_child;
        bvh_nodes[i].right_child = right_child;
        bvh_nodes[i].particle_idx = -1;  // -1 indicates an internal node

        // Assign parent pointers to children
        bvh_nodes[left_child].parent = i;
        bvh_nodes[right_child].parent = i;
    }

    // Set the root node's parent to -1
    bvh_nodes[0].parent = -1;

    // BOTTOM-UP AGGREGATION (Bounding Boxes & Center of Mass)
    // Counter for each internal node to track when both children are processed
    std::vector<int> atomic_flags(n - 1, 0);

// Initialize Leaf Nodes and trigger the walk up
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        int leaf_idx = n - 1 + i;
        double px = sorted_pos_x[i], py = sorted_pos_y[i], pz = sorted_pos_z[i];

        // Extract the smoothing length if the pointer is valid
        double radius = 0.0;
        if (sorted_h != nullptr) {
            radius = (*sorted_h)[i];
        }

        bvh_nodes[leaf_idx].bbox.min_x = px - radius;
        bvh_nodes[leaf_idx].bbox.max_x = px + radius;
        bvh_nodes[leaf_idx].bbox.min_y = py - radius;
        bvh_nodes[leaf_idx].bbox.max_y = py + radius;
        bvh_nodes[leaf_idx].bbox.min_z = pz - radius;
        bvh_nodes[leaf_idx].bbox.max_z = pz + radius;

        bvh_nodes[leaf_idx].mass = sorted_mass[i];

        // Store mass-weighted positions temporarily to make summation easy
        bvh_nodes[leaf_idx].com_x = px * sorted_mass[i];
        bvh_nodes[leaf_idx].com_y = py * sorted_mass[i];
        bvh_nodes[leaf_idx].com_z = pz * sorted_mass[i];

        // Walk up the tree
        int curr = bvh_nodes[leaf_idx].parent;
        while (curr != -1) {
            int old_flag;
#pragma omp atomic capture
            old_flag = atomic_flags[curr]++;

            if (old_flag == 0) {
                // First thread to arrive. The other child isn't ready yet.
                // Terminate
                break;
            }

            // Second thread to arrive. Both children are ready. Compute parent
            int left = bvh_nodes[curr].left_child;
            int right = bvh_nodes[curr].right_child;

            // Combine Bounding Boxes
            bvh_nodes[curr].bbox.min_x = std::min(bvh_nodes[left].bbox.min_x,
                                                  bvh_nodes[right].bbox.min_x);
            bvh_nodes[curr].bbox.max_x = std::max(bvh_nodes[left].bbox.max_x,
                                                  bvh_nodes[right].bbox.max_x);
            bvh_nodes[curr].bbox.min_y = std::min(bvh_nodes[left].bbox.min_y,
                                                  bvh_nodes[right].bbox.min_y);
            bvh_nodes[curr].bbox.max_y = std::max(bvh_nodes[left].bbox.max_y,
                                                  bvh_nodes[right].bbox.max_y);
            bvh_nodes[curr].bbox.min_z = std::min(bvh_nodes[left].bbox.min_z,
                                                  bvh_nodes[right].bbox.min_z);
            bvh_nodes[curr].bbox.max_z = std::max(bvh_nodes[left].bbox.max_z,
                                                  bvh_nodes[right].bbox.max_z);

            // Sum Mass and mass-weighted positions
            bvh_nodes[curr].mass = bvh_nodes[left].mass + bvh_nodes[right].mass;
            bvh_nodes[curr].com_x =
                bvh_nodes[left].com_x + bvh_nodes[right].com_x;
            bvh_nodes[curr].com_y =
                bvh_nodes[left].com_y + bvh_nodes[right].com_y;
            bvh_nodes[curr].com_z =
                bvh_nodes[left].com_z + bvh_nodes[right].com_z;

            // Move up to the next parent
            curr = bvh_nodes[curr].parent;
        }
    }

// Normalize Center of Mass
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 2 * n - 1; ++i) {
        if (bvh_nodes[i].mass > 0.0) {
            bvh_nodes[i].com_x /= bvh_nodes[i].mass;
            bvh_nodes[i].com_y /= bvh_nodes[i].mass;
            bvh_nodes[i].com_z /= bvh_nodes[i].mass;
        }
    }
}

void update_h_and_bboxes(size_t num_particles, const std::vector<double>& pos_x,
                         const std::vector<double>& pos_y,
                         const std::vector<double>& pos_z,
                         const std::vector<double>& h,
                         std::vector<BVHNode>& bvh_nodes) {
    if (num_particles == 0) return;
    const int n = static_cast<int>(num_particles);
    std::vector<int> atomic_flags(n - 1, 0);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        int leaf_idx = n - 1 + i;
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];
        double radius = h[i];

        // Update the Leaf
        bvh_nodes[leaf_idx].bbox.min_x = px - radius;
        bvh_nodes[leaf_idx].bbox.max_x = px + radius;
        bvh_nodes[leaf_idx].bbox.min_y = py - radius;
        bvh_nodes[leaf_idx].bbox.max_y = py + radius;
        bvh_nodes[leaf_idx].bbox.min_z = pz - radius;
        bvh_nodes[leaf_idx].bbox.max_z = pz + radius;

        // Walk up the tree
        int curr = bvh_nodes[leaf_idx].parent;
        while (curr != -1) {
            int old_flag;
#pragma omp atomic capture
            old_flag = atomic_flags[curr]++;

            if (old_flag == 0) {
                // First thread to arrive. The other child isn't ready yet
                break;
            }

            // Second thread to arrive. Both children are ready. Compute parent
            int left = bvh_nodes[curr].left_child;
            int right = bvh_nodes[curr].right_child;

            bvh_nodes[curr].bbox.min_x = std::min(bvh_nodes[left].bbox.min_x,
                                                  bvh_nodes[right].bbox.min_x);
            bvh_nodes[curr].bbox.max_x = std::max(bvh_nodes[left].bbox.max_x,
                                                  bvh_nodes[right].bbox.max_x);
            bvh_nodes[curr].bbox.min_y = std::min(bvh_nodes[left].bbox.min_y,
                                                  bvh_nodes[right].bbox.min_y);
            bvh_nodes[curr].bbox.max_y = std::max(bvh_nodes[left].bbox.max_y,
                                                  bvh_nodes[right].bbox.max_y);
            bvh_nodes[curr].bbox.min_z = std::min(bvh_nodes[left].bbox.min_z,
                                                  bvh_nodes[right].bbox.min_z);
            bvh_nodes[curr].bbox.max_z = std::max(bvh_nodes[left].bbox.max_z,
                                                  bvh_nodes[right].bbox.max_z);

            // Move up to the next parent
            curr = bvh_nodes[curr].parent;
        }
    }
}

void compute_morton_and_sort_indices_gpu(
    size_t num_particles, const Vec3& domain_size, const double* pos_x,
    const double* pos_y, const double* pos_z, uint64_t* morton_codes,
    int* sorted_indices) {
    if (num_particles == 0) return;

    double inv_domain_x = 1.0 / domain_size.x;
    double inv_domain_y = 1.0 / domain_size.y;
    double inv_domain_z = 1.0 / domain_size.z;
    double bound = 2097152.0;

    // Generate Morton codes and initialize indices on the GPU
#pragma omp target teams distribute parallel for map(                  \
        present : pos_x[0 : num_particles], pos_y[0 : num_particles],  \
            pos_z[0 : num_particles], morton_codes[0 : num_particles], \
            sorted_indices[0 : num_particles])
    for (size_t i = 0; i < num_particles; ++i) {
        uint32_t x = static_cast<uint32_t>(
            fmod(pos_x[i] * inv_domain_x + 1.0, 1.0) * bound);
        uint32_t y = static_cast<uint32_t>(
            fmod(pos_y[i] * inv_domain_y + 1.0, 1.0) * bound);
        uint32_t z = static_cast<uint32_t>(
            fmod(pos_z[i] * inv_domain_z + 1.0, 1.0) * bound);
        morton_codes[i] = morton3D(x, y, z);
        sorted_indices[i] = static_cast<int>(i);
    }

    // Pad to the next power of 2 (Required for Bitonic Sort networks)
    size_t padded_n = 1;
    while (padded_n < num_particles) padded_n <<= 1;

    // Fill virtual out-of-bounds space with maximum values so they sort to the
    // end
    std::vector<uint64_t> padded_morton(padded_n, UINT64_MAX);
    std::vector<int> padded_indices(padded_n, -1);

    uint64_t* d_padded_morton = padded_morton.data();
    int* d_padded_indices = padded_indices.data();

    // Map padded arrays to device
#pragma omp target enter data map(to : d_padded_morton[0 : padded_n], \
                                      d_padded_indices[0 : padded_n])

    // Copy actual particles into the padded target buffers
#pragma omp target teams distribute parallel for map(                         \
        present : morton_codes[0 : num_particles],                            \
            sorted_indices[0 : num_particles], d_padded_morton[0 : padded_n], \
            d_padded_indices[0 : padded_n])
    for (size_t i = 0; i < num_particles; ++i) {
        d_padded_morton[i] = morton_codes[i];
        d_padded_indices[i] = sorted_indices[i];
    }

    // Parallel Bitonic Sort (Executes in-place on VRAM without atomics)
    for (size_t k = 2; k <= padded_n; k <<= 1) {
        for (size_t j = k >> 1; j > 0; j >>= 1) {
#pragma omp target teams distribute parallel for map( \
        present : d_padded_morton[0 : padded_n],      \
            d_padded_indices[0 : padded_n])
            for (size_t i = 0; i < padded_n; ++i) {
                size_t ixj = i ^ j;

                // Ensure only one thread per pair does the comparison and swap
                if (ixj > i) {
                    bool ascending = ((i & k) == 0);
                    uint64_t val_i = d_padded_morton[i];
                    uint64_t val_ixj = d_padded_morton[ixj];

                    if ((val_i > val_ixj) == ascending) {
                        // Swap morton codes
                        d_padded_morton[i] = val_ixj;
                        d_padded_morton[ixj] = val_i;

                        // Swap original indices
                        int temp_idx = d_padded_indices[i];
                        d_padded_indices[i] = d_padded_indices[ixj];
                        d_padded_indices[ixj] = temp_idx;
                    }
                }
            }
        }
    }

    // Retrieve the valid sorted section back into the persistent arrays
#pragma omp target teams distribute parallel for map(                         \
        present : morton_codes[0 : num_particles],                            \
            sorted_indices[0 : num_particles], d_padded_morton[0 : padded_n], \
            d_padded_indices[0 : padded_n])
    for (size_t i = 0; i < num_particles; ++i) {
        morton_codes[i] = d_padded_morton[i];
        sorted_indices[i] = d_padded_indices[i];
    }

    // Clean up temporary padded buffers
#pragma omp target exit data map(delete : d_padded_morton[0 : padded_n], \
                                     d_padded_indices[0 : padded_n])
}

void build_topology_and_aggregate_gpu(size_t num_particles,
                                      const double* sorted_pos_x,
                                      const double* sorted_pos_y,
                                      const double* sorted_pos_z,
                                      const double* sorted_mass,
                                      const uint64_t* morton_codes,
                                      BVHNode* bvh_nodes) {
    if (num_particles == 0) return;
    const int n = static_cast<int>(num_particles);
    const int num_nodes = 2 * n - 1;

    // Initialize Leaf Nodes
#pragma omp target teams distribute parallel for map( \
        present : bvh_nodes[0 : num_nodes])
    for (int i = 0; i < n; ++i) {
        int leaf_idx = n - 1 + i;
        bvh_nodes[leaf_idx].particle_idx = i;
        bvh_nodes[leaf_idx].left_child = -1;
        bvh_nodes[leaf_idx].right_child = -1;
    }

    if (n < 2) return;

    // Construct Internal Nodes (Karras 2012)
#pragma omp target teams distribute parallel for map( \
        present : morton_codes[0 : n], bvh_nodes[0 : num_nodes])
    for (int i = 0; i < n - 1; ++i) {
        // Device-compatible delta function
        auto delta = [=](int idx_a, int idx_b) -> int {
            if (idx_b < 0 || idx_b >= n) return -1;
            uint64_t code_a = morton_codes[idx_a];
            uint64_t code_b = morton_codes[idx_b];
            if (code_a == code_b) {
                return 64 + __builtin_clzll(
                                static_cast<unsigned long long>(idx_a ^ idx_b));
            }
            return __builtin_clzll(code_a ^ code_b);
        };

        int d = (delta(i, i + 1) - delta(i, i - 1)) > 0 ? 1 : -1;
        int delta_min = delta(i, i - d);
        int l_max = 2;
        while (delta(i, i + l_max * d) > delta_min) l_max *= 2;

        int l = 0;
        for (int t = l_max / 2; t >= 1; t /= 2) {
            if (delta(i, i + (l + t) * d) > delta_min) l += t;
        }
        int j = i + l * d;
        int delta_node = delta(i, j);
        int s = 0;
        int t = l;
        do {
            t = (t + 1) >> 1;
            if (s + t < l && delta(i, i + (s + t) * d) > delta_node) s += t;
        } while (t > 1);

        int split = i + s * d + std::min(d, 0);
        int min_idx = std::min(i, j);
        int max_idx = std::max(i, j);

        int left_child = (min_idx == split) ? (n - 1 + split) : split;
        int right_child =
            (max_idx == split + 1) ? (n - 1 + split + 1) : (split + 1);

        bvh_nodes[i].left_child = left_child;
        bvh_nodes[i].right_child = right_child;
        bvh_nodes[i].particle_idx = -1;
        bvh_nodes[left_child].parent = i;
        bvh_nodes[right_child].parent = i;
    }

#pragma omp target map(present : bvh_nodes[0 : num_nodes])
    {
        bvh_nodes[0].parent = -1;
    }

    // BOTTOM-UP AGGREGATION
    // Allocate the atomic flags array on the device
    std::vector<int> host_atomic_flags(n - 1, 0);
    int* atomic_flags = host_atomic_flags.data();
    int flags_size = n - 1;

#pragma omp target enter data map(to : atomic_flags[0 : flags_size])

#pragma omp target teams distribute parallel for map(                          \
        present : sorted_pos_x[0 : n], sorted_pos_y[0 : n],                    \
            sorted_pos_z[0 : n], sorted_mass[0 : n], bvh_nodes[0 : num_nodes], \
            atomic_flags[0 : flags_size])
    for (int i = 0; i < n; ++i) {
        int leaf_idx = n - 1 + i;
        double px = sorted_pos_x[i], py = sorted_pos_y[i], pz = sorted_pos_z[i];

        bvh_nodes[leaf_idx].bbox.min_x = px;
        bvh_nodes[leaf_idx].bbox.max_x = px;
        bvh_nodes[leaf_idx].bbox.min_y = py;
        bvh_nodes[leaf_idx].bbox.max_y = py;
        bvh_nodes[leaf_idx].bbox.min_z = pz;
        bvh_nodes[leaf_idx].bbox.max_z = pz;

        bvh_nodes[leaf_idx].mass = sorted_mass[i];
        bvh_nodes[leaf_idx].com_x = px * sorted_mass[i];
        bvh_nodes[leaf_idx].com_y = py * sorted_mass[i];
        bvh_nodes[leaf_idx].com_z = pz * sorted_mass[i];

        int curr = bvh_nodes[leaf_idx].parent;
        while (curr != -1) {
            int old_flag;
#pragma omp atomic capture
            old_flag = atomic_flags[curr]++;

            if (old_flag == 0) break;  // First thread terminates

            int left = bvh_nodes[curr].left_child;
            int right = bvh_nodes[curr].right_child;

            bvh_nodes[curr].bbox.min_x =
                (bvh_nodes[left].bbox.min_x < bvh_nodes[right].bbox.min_x)
                    ? bvh_nodes[left].bbox.min_x
                    : bvh_nodes[right].bbox.min_x;
            bvh_nodes[curr].bbox.max_x =
                (bvh_nodes[left].bbox.max_x > bvh_nodes[right].bbox.max_x)
                    ? bvh_nodes[left].bbox.max_x
                    : bvh_nodes[right].bbox.max_x;
            bvh_nodes[curr].bbox.min_y =
                (bvh_nodes[left].bbox.min_y < bvh_nodes[right].bbox.min_y)
                    ? bvh_nodes[left].bbox.min_y
                    : bvh_nodes[right].bbox.min_y;
            bvh_nodes[curr].bbox.max_y =
                (bvh_nodes[left].bbox.max_y > bvh_nodes[right].bbox.max_y)
                    ? bvh_nodes[left].bbox.max_y
                    : bvh_nodes[right].bbox.max_y;
            bvh_nodes[curr].bbox.min_z =
                (bvh_nodes[left].bbox.min_z < bvh_nodes[right].bbox.min_z)
                    ? bvh_nodes[left].bbox.min_z
                    : bvh_nodes[right].bbox.min_z;
            bvh_nodes[curr].bbox.max_z =
                (bvh_nodes[left].bbox.max_z > bvh_nodes[right].bbox.max_z)
                    ? bvh_nodes[left].bbox.max_z
                    : bvh_nodes[right].bbox.max_z;

            bvh_nodes[curr].mass = bvh_nodes[left].mass + bvh_nodes[right].mass;
            bvh_nodes[curr].com_x =
                bvh_nodes[left].com_x + bvh_nodes[right].com_x;
            bvh_nodes[curr].com_y =
                bvh_nodes[left].com_y + bvh_nodes[right].com_y;
            bvh_nodes[curr].com_z =
                bvh_nodes[left].com_z + bvh_nodes[right].com_z;

            curr = bvh_nodes[curr].parent;
        }
    }

#pragma omp target exit data map(delete : atomic_flags[0 : flags_size])

// Normalize Center of Mass on GPU
#pragma omp target teams distribute parallel for map( \
        present : bvh_nodes[0 : num_nodes])
    for (int i = 0; i < num_nodes; ++i) {
        if (bvh_nodes[i].mass > 0.0) {
            bvh_nodes[i].com_x /= bvh_nodes[i].mass;
            bvh_nodes[i].com_y /= bvh_nodes[i].mass;
            bvh_nodes[i].com_z /= bvh_nodes[i].mass;
        }
    }
}

}  // namespace LBVH
