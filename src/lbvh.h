#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include "math_utils.h"
#include "types.h"

struct BoundingBox {
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
};

struct BVHNode {
    int parent;
    int left_child;
    int right_child;
    int particle_idx;  // -1 if internal node

    BoundingBox bbox;

    double mass;
    double com_x, com_y, com_z;
};

namespace LBVH {

// Step 1: Computes Morton codes based on 3D positions and generates the sorting
// index map
void compute_morton_and_sort_indices(size_t num_particles,
                                     const Vec3& domain_size,
                                     const std::vector<double>& pos_x,
                                     const std::vector<double>& pos_y,
                                     const std::vector<double>& pos_z,
                                     std::vector<uint64_t>& morton_codes,
                                     std::vector<int>& sorted_indices);

// Step 2: Wires the Karras tree topology and aggregates bounding boxes/masses
// bottom-up Note: The pos, mass, and h arrays passed here MUST be the newly
// sorted arrays
void build_topology_and_aggregate(
    size_t num_particles, const std::vector<double>& sorted_pos_x,
    const std::vector<double>& sorted_pos_y,
    const std::vector<double>& sorted_pos_z,
    const std::vector<double>& sorted_mass,
    const std::vector<double>* sorted_h,  // Pass nullptr for DM
    std::vector<uint64_t>& morton_codes, std::vector<BVHNode>& bvh_nodes);

void update_h_and_bboxes(size_t num_particles, const std::vector<double>& pos_x,
                         const std::vector<double>& pos_y,
                         const std::vector<double>& pos_z,
                         const std::vector<double>& h,
                         std::vector<BVHNode>& bvh_nodes);

#ifdef USE_GPU
#pragma omp declare target
#endif
template <typename LeafOperation>
inline void traverse(const BVHNode* bvh_nodes, const double* target_pos_x,
                     const double* target_pos_y, const double* target_pos_z,
                     double px, double py, double pz, double search_sq,
                     const Vec3& domain_size, LeafOperation&& leaf_op) {
    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = bvh_nodes[node_idx];

        double dist_sq = min_periodic_dist_sq(px, node.bbox.min_x,
                                              node.bbox.max_x, domain_size.x) +
                         min_periodic_dist_sq(py, node.bbox.min_y,
                                              node.bbox.max_y, domain_size.y) +
                         min_periodic_dist_sq(pz, node.bbox.min_z,
                                              node.bbox.max_z, domain_size.z);

        if (dist_sq > search_sq) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;

            double dx =
                periodic_displacement(target_pos_x[j] - px, domain_size.x);
            double dy =
                periodic_displacement(target_pos_y[j] - py, domain_size.y);
            double dz =
                periodic_displacement(target_pos_z[j] - pz, domain_size.z);
            double r2 = dx * dx + dy * dy + dz * dz;

            // Call the lambda
            leaf_op(j, dx, dy, dz, r2);
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }
}
#ifdef USE_GPU
#pragma omp end declare target
#endif

// GPU specific versions
#ifdef GRAVITY_GPU_MEMORY

class Workspace {
   public:
    size_t current_padded_capacity = 0;
    uint64_t* d_padded_morton = nullptr;
    int* d_padded_indices = nullptr;

    size_t current_flags_capacity = 0;
    int* d_atomic_flags = nullptr;

    Workspace() = default;
    ~Workspace();

    // Prevent accidental copying which would cause double-free errors
    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;
};

void compute_morton_and_sort_indices_gpu(
    Workspace& workspace, size_t num_particles, const Vec3& domain_size,
    const double* pos_x, const double* pos_y, const double* pos_z,
    uint64_t* morton_codes, int* sorted_indices);

void build_topology_and_aggregate_gpu(
    Workspace& workspace, size_t num_particles, const double* sorted_pos_x,
    const double* sorted_pos_y, const double* sorted_pos_z,
    const double* sorted_mass, const uint64_t* morton_codes,
    BVHNode* bvh_nodes);
#endif
}  // namespace LBVH
