#include <filesystem>
#include <iostream>
#include <aris.hpp>
#include <cmath>
#include <chrono>


#include "../../src/mpc/utpr_mpc.hpp"
#include "../../src/zmq/zmqmsg.hpp"

const double M_PI = 3.141592653589793;


double deg2rad(double degrees);

int main() {

    Zmqmsg zmqmsg;
    zmqmsg.init();

    UtprMpc utpr;

    int nx = 6;
    int nu = 2;
    double T = 0.01;
    int N = 100;
    int iter_max = 2;

    //Eigen::setNbThreads(4);

    Eigen::VectorXd u_bar(2,1);
    u_bar << 20, 20;             //

    Eigen::MatrixXd Q(nx, nx);
    Q.setZero();
    Q.diagonal() << 2000, 1000, 500, 2000, 100, 100;
    
    Eigen::MatrixXd R(nu, nu);
    R.setZero();
    R.diagonal() << 10, 8;

    Eigen::VectorXd ref(nx, 1);
    ref << 2* M_PI, 0, 0, 0, 0, 0;

    utpr.setIter_max(iter_max);
    utpr.setT(T);
    utpr.setN(N);
    utpr.setUBar(u_bar);
    utpr.setQ(Q);
    utpr.setR(R);
    utpr.setref(ref);

    std::array<double, 6> initialState = { deg2rad(0.6), deg2rad(0.6), deg2rad(0.6), 0, 0, 0 };
    utpr.setInitialState(initialState);
    utpr.InitialMPCParameters();

    std::array<double, 6> state{0};

    std::vector<double> data{0};
    std::vector<double> torque{0., 0.};

    while (true) {
        data = zmqmsg.get_request();
        //std::cout << " --------state------- " << std::endl;

        for (int i = 0; i < 6; ++i) {
            state[i] = data[i];
            //std::cout << state[i] << " " << std::endl;
        }

        torque = utpr.compute(state);

        zmqmsg.send_msg(torque);

    }

    return 0;
}


double deg2rad(double degrees) {
    return degrees * M_PI / 180.0;
}
