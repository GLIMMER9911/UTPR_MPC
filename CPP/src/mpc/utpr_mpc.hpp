#ifndef UTPR_MPC_HPP
#define UTPR_MPC_HPP

#include <array>
#include <vector>
#include <iostream>

#define EIGEN_USE_THREADS
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#include <unsupported/Eigen/KroneckerProduct> 
#include <aris.hpp>


#include "../solver/qpsolver.hpp"

using namespace Eigen;

//#define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(x) do { std::cout << x << std::endl; } while (0)
#else
#define DEBUG_PRINT(x) do {} while (0)
#endif

class UtprMpc {
public:
    UtprMpc();
    void initializeParameters();
    void InitialMPCParameters();
    void setInitialState(const std::array<double, 6>& state);
    auto compute(const std::array<double, 6>& state) -> std::vector<double>;

    void setQ(const Eigen::MatrixXd& Q);
    void setR(const Eigen::MatrixXd& R);
    void setUBar(const Eigen::VectorXd& u_bar);
    void setref(const Eigen::VectorXd& ref);

    void setT(double T);
    void setN(int N);
    void setIter_max(int iter_max);

private:
    void buildModelMatrices(const Eigen::VectorXd& x);
    void solveMpc(const Eigen::VectorXd& x);
    double si(double x);

    void fprintMatrix(int index, const Eigen::MatrixXd& H, const Eigen::MatrixXd& g);
    void fprintRho(int index, const Eigen::MatrixXd& rho);
    void fprintResult(int index, const Eigen::MatrixXd& result);



    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;
    Eigen::MatrixXd Phi;
    Eigen::MatrixXd Gamma;
    Eigen::VectorXd x0;
    Eigen::VectorXd ref;


    Eigen::MatrixXd g;
    Eigen::MatrixXd H;
    Eigen::MatrixXd Q_hat;
    Eigen::MatrixXd R_hat;
    Eigen::MatrixXd Acons;
    Eigen::MatrixXd Acons2;
    Eigen::MatrixXd bcons;
    Eigen::MatrixXd ucons;
    Eigen::MatrixXd lcons;

    
    Eigen::VectorXd rho;
    
    double T;
    int nx;
    int nu;
    int N;
    int iter_max;
    unsigned int count_ = 0;
    Eigen::VectorXd u_bar;
    std::vector<double> torque;


};

#endif // UTPR_MPC_HPP
