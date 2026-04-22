#include "types.h"
#include <omp.h>

Grid3D Grid3D::roll(int shift, int axis) {
    Grid3D res(n);
    
    // Collapse flattens the 3 nested loops into 1 for perfect thread distribution
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (axis == 0)
                    res((i + shift + n) % n, j, k) = (*this)(i, j, k);
                else if (axis == 1)
                    res(i, (j + shift + n) % n, k) = (*this)(i, j, k);
                else if (axis == 2)
                    res(i, j, (k + shift + n) % n) = (*this)(i, j, k);
            }
        }
    }
    return res;
}