#pragma once
#include <cstddef>  // Adds size_t
#include <cstdint>
#include <vector>

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
    double max_h;

    double mass;
    double com_x, com_y, com_z;
};

namespace LBVH {

// Step 1: Computes Morton codes based on 3D positions and generates the sorting
// index map
void compute_morton_and_sort_indices(size_t num_particles, double domain_size,
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
}  // namespace LBVH