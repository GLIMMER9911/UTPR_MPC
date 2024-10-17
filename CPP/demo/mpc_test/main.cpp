#include <Eigen/Dense>
#include <cmath>
#include <iostream>

// Define constants
const double g = 9.81;
const double l1 = 0.325;
const double lg1 = 0.245436;
const double l2 = 0.2;
const double lg2 = 0.165468;
const double l3 = 0.2;
const double lg3 = 0.127337;
const double m1 = 1.915;
const double m2 = 1.469;
const double m3 = 1.141;
const double J1 = 0.025720752;
const double J2 = 0.007667511;
const double J3 = 0.004949014;

// Global matrices
Eigen::MatrixXd A(6, 6);
Eigen::MatrixXd B(6, 2);
double T;

// SI function equivalent to sin(x) / x with singularity handling
double si(double x) {
    return (x == 0) ? 1.0 : std::sin(x) / x;
}

void utprLpv(const Eigen::VectorXd& x) {
    // Dissect states
    Eigen::Vector3d theta = x.segment<3>(0);
    Eigen::Vector3d omega = x.segment<3>(3);

    // Intermediate variables
    double a1 = J1 + m1 * lg1 * lg1 + (m2 + m3) * l1 * l1;
    double a2 = J2 + m2 * lg2 * lg2 + m3 * l2 * l2;
    double a3 = (m2 * lg2 + m3 * l2) * l1;
    double a4 = J3 + m3 * lg3 * lg3;
    double a5 = m3 * l1 * lg3;
    double a6 = m3 * l2 * lg3;

    double b1 = (m1 * lg1 + m2 * l1 + m3 * l1) * g;
    double b2 = (m2 * lg2 + m3 * l2) * g;
    double b3 = (m3 * lg3) * g;

    // Build preliminary matrices
    Eigen::Matrix3d M;
    M << a1 + a2 + a4 + 2 * a5 * std::cos(theta(1) + theta(2)) + 2 * a3 * std::cos(theta(1)) + 2 * a6 * std::cos(theta(2)),
        a2 + a4 + a3 * std::cos(theta(1)) + a5 * std::cos(theta(1) + theta(2)) + 2 * a6 * std::cos(theta(2)),
        a4 + a5 * std::cos(theta(1) + theta(2)) + a6 * std::cos(theta(2)),
        a2 + a4 + a3 * std::cos(theta(1)) + a5 * std::cos(theta(1) + theta(2)) + 2 * a6 * std::cos(theta(2)),
        a2 + a4 + 2 * a6 * std::cos(theta(2)),
        a4 + a6 * std::cos(theta(2)),
        a4 + a5 * std::cos(theta(1) + theta(2)) + a6 * std::cos(theta(2)),
        a4 + a6 * std::cos(theta(2)),
        a4;

    Eigen::Matrix3d C;
    C << -a5 * (omega(1) + omega(2)) * std::sin(theta(1) + theta(2)) - a3 * omega(1) * std::sin(theta(1)) - a6 * omega(2) * std::sin(theta(2)),
        -a5 * (omega(0) + omega(1) + omega(2)) * std::sin(theta(1) + theta(2)) - a3 * (omega(0) + omega(1)) * std::sin(theta(1)) - a6 * omega(2) * std::sin(theta(2)),
        -(omega(0) + omega(1) + omega(2)) * (a5 * std::sin(theta(1) + theta(2)) + a6 * std::sin(theta(2))),
        a5* omega(0)* std::sin(theta(1) + theta(2)) + a3 * omega(0) * std::sin(theta(1)) - a6 * omega(2) * std::sin(theta(2)),
        -a6 * omega(2) * std::sin(theta(2)),
        -a6 * (omega(0) + omega(1) + omega(2)) * std::sin(theta(2)),
        a5* omega(0)* std::sin(theta(1) + theta(2)) + a6 * (omega(0) + omega(1)) * std::sin(theta(2)),
        a6* (omega(0) + omega(1))* std::sin(theta(2)),
        0;

    Eigen::Matrix3d K;
    K << -b1 * si(theta(0)) - b2 * si(theta(0) + theta(1)) - b3 * si(theta(0) + theta(1) + theta(2)),
        -b2 * si(theta(0) + theta(1)) - b3 * si(theta(0) + theta(1) + theta(2)),
        -b3 * si(theta(0) + theta(1) + theta(2)),
        -b2 * si(theta(0) + theta(1)) - b3 * si(theta(0) + theta(1) + theta(2)),
        -b2 * si(theta(0) + theta(1)) - b3 * si(theta(0) + theta(1) + theta(2)),
        -b3 * si(theta(0) + theta(1) + theta(2)),
        -b3 * si(theta(0) + theta(1) + theta(2)),
        -b3 * si(theta(0) + theta(1) + theta(2)),
        -b3 * si(theta(0) + theta(1) + theta(2));

    // Compute the inverse of M
    Eigen::Matrix3d Mi = M.inverse();

    // Assemble model matrices A & B
    A.block<3, 3>(0, 0) = Eigen::Matrix3d::Zero();
    A.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    A.block<3, 3>(3, 0) = -Mi * K;
    A.block<3, 3>(3, 3) = -Mi * C;

    B.block<3, 1>(0, 0) = Eigen::Vector3d::Zero();
    B.block<3, 1>(0, 1) = Eigen::Vector3d::Zero();
    B.block<3, 1>(3, 0) = Mi.col(1);
    B.block<3, 1>(3, 1) = Mi.col(2);

    // Discretize if sampling time T is given
    if (T > 0) {
        A = Eigen::MatrixXd::Identity(6, 6) + T * A;
        B = T * B;
    }
}

int main() {
    // Example usage
    Eigen::VectorXd x(6);
    x << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
    T = 0.01;

    utprLpv(x);

    std::cout << "A matrix:\n" << A << std::endl;
    std::cout << "B matrix:\n" << B << std::endl;

    return 0;
}
