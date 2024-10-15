#include <vector>
#include <cmath>

// utpr_lpv 函数
void utpr_lpv(const std::vector<double>& x, double T, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B) {
    // 定义物理参数
    double l_1 = 0.34;
    double lg_1 = 0.17;
    double l_2 = 0.29;
    double lg_2 = 0.145;
    double l_3 = 0.52;
    double lg_3 = 0.26;
    double m1 = 1.258;
    double m2 = 5.686;
    double m3 = 2.162;
    double J1 = 0.0121;
    double J2 = 0.0398;
    double J3 = 0.0487;
    double g = 9.81;

    // 从状态变量中提取角度和角速度
    std::vector<double> theta = {x[0], x[1], x[2]};
    std::vector<double> omega = {x[3], x[4], x[5]};

    // 计算中间变量
    double a1 = J1 + m1 * lg_1 * lg_1 + (m2 + m3) * l_1 * l_1;
    double a2 = J2 + m2 * lg_2 * lg_2 + m3 * l_2 * l_2;
    double a3 = (m2 * lg_2 + m3 * l_2) * l_1;
    double a4 = J3 + m3 * lg_3 * lg_3;
    double a5 = m3 * l_1 * lg_3;
    double a6 = m3 * l_2 * lg_3;

    double b1 = (m1 * lg_1 + m2 * l_1 + m3 * l_1) * g;
    double b2 = (m2 * lg_2 + m3 * l_2) * g;
    double b3 = (m3 * lg_3) * g;

    // 初始化矩阵 A 和 B
    A.resize(3, std::vector<double>(3));
    B.resize(3, std::vector<double>(3));

    // 计算矩阵 A 和 B
    A[0][0] = a1 + a2 + a4 + 2 * a5 * cos(theta[1] + theta[2]) + 2 * a3 * cos(theta[1]) + 2 * a6 * cos(theta[2]);
    A[0][1] = a2 + a4 + a3 * cos(theta[1]) + a5 * cos(theta[1] + theta[2]) + 2 * a6 * cos(theta[2]);
    A[0][2] = a4 + a5 * cos(theta[1] + theta[2]) + a6 * cos(theta[2]);
    A[1][0] = A[0][1];
    A[1][1] = a2 + a4 + 2 * a6 * cos(theta[2]);
    A[1][2] = a4 + a6 * cos(theta[2]);
    A[2][0] = A[0][2];
    A[2][1] = A[1][2];
    A[2][2] = a4;

    B[0][0] = -b1 * sin(theta[0]) - b2 * sin(theta[0] + theta[1]) - b3 * sin(theta[0] + theta[1] + theta[2]);
    B[1][0] = -b2 * sin(theta[0] + theta[1]) - b3 * sin(theta[0] + theta[1] + theta[2]);
    B[2][0] = -b3 * sin(theta[0] + theta[1] + theta[2]);
}
