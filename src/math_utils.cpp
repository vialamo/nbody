#include "math_utils.h"

void QuadraticInterpolator::Set(const Eigen::Vector2d& p1,
                                const Eigen::Vector2d& p2,
                                const Eigen::Vector2d& p3) {
    Eigen::Matrix3d A;
    A << p1.x() * p1.x(), p1.x(), 1.0, p2.x() * p2.x(), p2.x(), 1.0,
        p3.x() * p3.x(), p3.x(), 1.0;

    Eigen::Vector3d y;
    y << p1.y(), p2.y(), p3.y();

    // Solve the linear system A * coeffs = y
    Eigen::Vector3d coeffs = A.colPivHouseholderQr().solve(y);

    a = coeffs(0);
    b = coeffs(1);
    c = coeffs(2);
}
