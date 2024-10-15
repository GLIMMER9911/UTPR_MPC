#ifndef UTPR_LPV_HPP_
#define UTPR_LPV_HPP_

#include <vector>
#include <cmath>

void utpr_lpv(const std::vector<double>& x, double T, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B);

#endif
