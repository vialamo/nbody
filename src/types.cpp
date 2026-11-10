#include "types.h"

#include <omp.h>

Grid3D Grid3D::roll(int shift, int axis) const {
    Grid3D res(Vec3i{nx, ny, nz});

// Collapse flattens the 3 nested loops into 1 for thread distribution
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                if (axis == 0) {
                    int new_i = ((i + shift) % nx + nx) % nx;
                    res(new_i, j, k) = (*this)(i, j, k);
                } else if (axis == 1) {
                    int new_j = ((j + shift) % ny + ny) % ny;
                    res(i, new_j, k) = (*this)(i, j, k);
                } else if (axis == 2) {
                    int new_k = ((k + shift) % nz + nz) % nz;
                    res(i, j, new_k) = (*this)(i, j, k);
                }
            }
        }
    }
    return res;
}
