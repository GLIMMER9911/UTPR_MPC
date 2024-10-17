#include "utpr_mpc.hpp"


UtprMpc::UtprMpc() : nx(6), nu(2), N(10), T(0.1), iter_max(4) {
	this->count_ = 0;

	initializeParameters();
}

void UtprMpc::initializeParameters() {
	A = Eigen::MatrixXd::Zero(nx, nx);
	B = Eigen::MatrixXd::Zero(nx, nu);
	Q = Eigen::MatrixXd::Identity(nx, nx);
	R = Eigen::MatrixXd::Identity(nu, nu);
	Phi = Eigen::MatrixXd::Zero(N * nx, nx);
	Gamma = Eigen::MatrixXd::Zero(N * nx, N * nu);
	u_bar = Eigen::VectorXd::Constant(nu, 1.0);
	x0 = Eigen::VectorXd::Zero(nx, 1);
	ref = Eigen::VectorXd::Zero(nx, 1);
	rho = Eigen::VectorXd::Zero(N * nx);
	H = Eigen::MatrixXd::Zero(N * nu, N * nu);
	g = Eigen::MatrixXd::Zero(N * nu, 1);

	torque.resize(2);

	DEBUG_PRINT("A matrix: " << A.rows() << " " << A.cols());
	DEBUG_PRINT("B matrix: " << B.rows() << " " << B.cols());
	DEBUG_PRINT("Q matrix: " << Q.rows() << " " << Q.cols());
	DEBUG_PRINT("R matrix: " << R.rows() << " " << R.cols());
	DEBUG_PRINT("Phi matrix: " << Phi.rows() << " " << Phi.cols());
	DEBUG_PRINT("Gamma matrix: " << Gamma.rows() << " " << Gamma.cols());

	//std::cout << "u_bar matrix: " << u_bar.rows() << " " << u_bar.cols() << " " << u_bar << std::endl;

}

void UtprMpc::setInitialState(const std::array<double, 6>& state) {
	x0 = Eigen::VectorXd::Map(state.data(), state.size());

	DEBUG_PRINT("------------- setInitialState -------------");
	DEBUG_PRINT("x0: " << x0.rows() << " " << x0.cols());
	DEBUG_PRINT(x0);

}

void UtprMpc::setQ(const Eigen::MatrixXd& Q) {
	this->Q = Q;

	DEBUG_PRINT("------------- setQ -------------");
	DEBUG_PRINT("Q: " << Q.rows() << " " << Q.cols());
	DEBUG_PRINT(Q);
}

void UtprMpc::setR(const Eigen::MatrixXd& R) {
	this->R = R;

	DEBUG_PRINT("------------- setR -------------");
	DEBUG_PRINT("R: " << R.rows() << " " << R.cols());
	DEBUG_PRINT(this->R);
}

void UtprMpc::setT(const double T) {
	this->T = T;

	DEBUG_PRINT("------------- setT -------------");
	DEBUG_PRINT("T: " << this->T);
}

void UtprMpc::setN(const int N) {
	this->N = N;
	this->Phi = Eigen::MatrixXd::Zero(N * nx, nx);
	this->Gamma = Eigen::MatrixXd::Zero(N * nx, N * nu);

	DEBUG_PRINT("------------- setN -------------");
	DEBUG_PRINT("N: " << this->N);
	DEBUG_PRINT("Phi matrix: " << Phi.rows() << " " << Phi.cols());
	DEBUG_PRINT("Gamma matrix: " << Gamma.rows() << " " << Gamma.cols());
}

void UtprMpc::setUBar(const Eigen::VectorXd& u_bar) {
	this->u_bar = u_bar;

	DEBUG_PRINT("------------- setUBar -------------");
	DEBUG_PRINT("u_bar matrix: " << u_bar.rows() << " " << u_bar.cols());
	DEBUG_PRINT(this->u_bar);
}

void UtprMpc::setref(const Eigen::VectorXd& ref) {
	this->ref = ref;

	DEBUG_PRINT("------------- setUBar -------------");
	DEBUG_PRINT("u_bar matrix: " << ref.rows() << " " << ref.cols());
	DEBUG_PRINT(this->ref);
}

void UtprMpc::setIter_max(int iter_max) {
	this->iter_max = iter_max;

	DEBUG_PRINT("------------- setIter_max -------------");
	DEBUG_PRINT("iter_max: " << this->iter_max);
}

auto UtprMpc::compute(const std::array<double, 6>& state) -> std::vector<double> {

	Eigen::VectorXd x = Eigen::VectorXd::Zero(nx, 1);
	if (this->count_ == 0) {
		x = x0;
		Eigen::MatrixXd rho_mtr = kroneckerProduct(Eigen::VectorXd::Ones(N + 1, 1), x).eval();
		this->rho = Eigen::Map<Eigen::VectorXd>(rho_mtr.data(), rho_mtr.size());
	}
	else {
		x = Eigen::VectorXd::Map(state.data(), state.size());
	}

	auto start_time = std::chrono::high_resolution_clock::now();
	solveMpc(x);

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_time - start_time;
	//std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

	this->count_++;
	return torque;
}

void UtprMpc::InitialMPCParameters() {
	this->Q_hat = kroneckerProduct(Eigen::MatrixXd::Identity(N, N), Q).eval();
	this->R_hat = kroneckerProduct(Eigen::MatrixXd::Identity(N, N), R).eval();
	this->Acons = kroneckerProduct(Eigen::MatrixXd::Identity(N, N), Eigen::MatrixXd::Identity(nu, nu)).eval();
	
	MatrixXd eye_nu = MatrixXd::Identity(nu, nu);
	MatrixXd minus_eye_nu = -MatrixXd::Identity(nu, nu);
	MatrixXd eye_combined(2 * nu, nu);
	eye_combined << eye_nu, minus_eye_nu;

	this->Acons2 = kroneckerProduct(Eigen::MatrixXd::Identity(N, N), eye_combined).eval();
	this->ucons = kroneckerProduct(Eigen::MatrixXd::Ones(N, 1), u_bar).eval();
	this->lcons = -this->ucons;
	this->bcons = kroneckerProduct(Eigen::VectorXd::Ones(2 * N, 1), u_bar).eval();

	//DEBUG_PRINT("------------- InitialMPCParameters -------------");
	//DEBUG_PRINT("Acons: " << this->Acons.rows() << " " << this->Acons.cols());

	//std::cout << "Q_hat matrix: " << Q_hat.rows() << " " << Q_hat.cols() << std::endl;
	//std::cout << "R_hat matrix: " << R_hat.rows() << " " << R_hat.cols() << std::endl;
	//std::cout << "Acons matrix: " << Acons.rows() << " " << Acons.cols() << std::endl;
	//std::cout << "lcons matrix: " << lcons.rows() << " " << lcons.cols() << std::endl;
}

void UtprMpc::buildModelMatrices(const Eigen::VectorXd& x) {
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
	//DEBUG_PRINT("A matrix: " << A);
	//DEBUG_PRINT("B matrix: " << B);

}

void UtprMpc::solveMpc(const Eigen::VectorXd& x) {

	Eigen::setNbThreads(4);

	// 设置QP问题
	const int n = N * nu; // 变量数目
	const int nCE = 0, nCI = 2 * n;

	//DEBUG_PRINT("x : " << x);

	double cost_mpc = 0;

	VectorXd x_mpc = x;
	VectorXd x_hat = ref - x;

	Eigen::MatrixXd rho_conv = MatrixXd::Zero(rho.size(), iter_max + 1);
	rho_conv.col(0) = rho;

	MatrixXd result = MatrixXd::Zero(N * nu, 1);

	//this->fprintRho(this->count_, this->rho);

	for (int j = 0; j < iter_max; ++j) {
		for (int i = 0; i < N; ++i) {
			int x_range_start = i * nx;
			int x_range_end = x_range_start + nx;
			int u_range_start = i * nu;
			int u_range_end = u_range_start + nu;
			VectorXd x_slice = rho.segment(x_range_start, nx);
			buildModelMatrices(x_slice);

			if (i == 0) {
				Phi.block(x_range_start, 0, nx, nx) = A;
			}
			else {
				Phi.block(x_range_start, 0, nx, nx) = A * Phi.block(x_range_start - nx, 0, nx, nx);
			}
			Gamma.block(x_range_start, u_range_start, nx, nu) = B;
			for (int l = 0; l < i; ++l) {
				int z_range_start = l * nu;
				int z_range_end = z_range_start + nu;
				Gamma.block(x_range_start, z_range_start, nx, nu) = A * Gamma.block(x_range_start - nx, z_range_start, nx, nu);
			}
		}

		double max_value = 1e20;
		Gamma = Gamma.cwiseMin(max_value);
		Gamma = Gamma.cwiseMax(-max_value);

		Phi = Phi.cwiseMin(max_value);
		Phi = Phi.cwiseMax(-max_value);

		double regularization = 1e-10;
		Gamma = Gamma / (1 + regularization * Gamma.norm());

		Phi = Phi / (1 + regularization * Phi.norm());

		//auto start_time = std::chrono::high_resolution_clock::now();

		//Eigen::setNbThreads(4);
		Eigen::SparseMatrix<double> Gamma_sparse = Gamma.sparseView();
		Eigen::SparseMatrix<double> Phi_sparse = Phi.sparseView();
		Eigen::SparseMatrix<double> Q_hat_sparse = Q_hat.sparseView();

		g.noalias() = Gamma_sparse.transpose() * Q_hat_sparse * Phi_sparse * x_mpc;
		H = Gamma_sparse.transpose() * Q_hat_sparse * Gamma_sparse + R_hat;
		H = 0.5 * (H + H.transpose());

		//auto end_time = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed = end_time - start_time;
		//std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

		//DEBUG_PRINT("Phi matrix row col: " << Phi.rows() << " " << Phi.cols());
		//DEBUG_PRINT("Phi matrix: " << Phi);
		//DEBUG_PRINT("Gamma matrix row col: " << Gamma.rows() << " " << Gamma.cols());
		//DEBUG_PRINT("Gamma matrix: \n" << Gamma);

		//DEBUG_PRINT("===== current iter_max : " << j);
		// OSQP solver
		// n is the size of output, m is the size of input, p is P, q is Q, a is A, l and u
		//DEBUG_PRINT("===== current iter_max : " << j);
		//result = qpSolver(H, g, this->Acons, this->lcons, this->ucons, n, 0);
		//DEBUG_PRINT("===== result : " << result(0) << " " << result(1) );



		//this->fprintResult(j, result);

		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> H_row_major(H.data(), H.rows(), H.cols());

		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> g_row_major(g.data(), g.rows(), g.cols());

		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Acons2_row(Acons2.data(), Acons2.rows(), Acons2.cols());

		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> bcons_row(bcons.data(), bcons.rows(), bcons.cols());
		// double* aris_b = bcons_row.data();

		std::vector<double> aris_CE(nCE * n);
		std::vector<double> aris_ce(nCE);

		std::vector<double> res(n), mem(n * n * 2 + 8 * (nCE + nCI) + 3 * n);
		auto r = aris::dynamic::s_quadprog(n, nCE, nCI, H_row_major.data(), g_row_major.data(), aris_CE.data(), aris_ce.data(), Acons2_row.data(), bcons_row.data(), res.data(), mem.data());

		for (int i = 0; i < res.size(); ++i) {
			result(i, 0) = res[i];
		}

		rho.head(nx) = x_mpc.col(0);
		rho.segment(nx, N * nx) = Phi * x_mpc.col(0) + Gamma * result;
		rho_conv.col(j + 1) = rho;

	}

	this->torque[0] = result(0);
	this->torque[1] = result(1);

	if (this->count_ % 1 == 0) {
		std::cout << " ---------------- cout " << this->count_ << " ---------------- " << std::endl;
		std::cout << "torque: " << result(0) << " " << result(1) << std::endl;
	}

	int size_rho = this->rho.size();
	if (size_rho > nx) {
		auto el_size = size_rho - nx;
		this->rho.segment(0, el_size) = this->rho.segment(nx, el_size);
		this->rho.tail(nx).setZero();
	}

	//this->fprintRho(this->count_, this->rho);
	//this->fprintMatrix(this->count_, H, g);

}

void UtprMpc::fprintMatrix(int index, const Eigen::MatrixXd& H, const Eigen::MatrixXd& g)
{
	std::ofstream file0("matrix.txt", std::ios::app); // 打开文件以进行写入
	if (file0.is_open()) {
		file0 << "................ " << index << " ................" << std::endl;
		file0 << "H" << std::endl;
		file0 << H << std::endl; // 使用Eigen的输出运算符
		file0 << "g" << std::endl;
		file0 << g << std::endl;
		file0 << "R_hat" << std::endl;
		file0 << R_hat << std::endl;
		file0 << "Q_hat" << std::endl;
		file0 << Q_hat << std::endl;
		file0 << "Gamma" << std::endl;
		file0 << Gamma << std::endl; // 使用Eigen的输出运算符
		file0 << "Phi" << std::endl;
		file0 << Phi << std::endl; // 使用Eigen的输出运算符
			
		file0.close(); // 关闭文件
		//std::cout << "Matrix saved to matrix.txt" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing" << std::endl;
	}
}

void UtprMpc::fprintRho(int index, const Eigen::MatrixXd& rho)
{
	std::ofstream file("rho_matrix.txt", std::ios::app); // 打开文件以进行写入
	if (file.is_open()) {
		file << "................ " << index << " ................" << std::endl;
		file << rho << std::endl; // 使用Eigen的输出运算符

		file.close(); // 关闭文件
		//std::cout << "Matrix saved to matrix.txt" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing" << std::endl;
	}
}

void UtprMpc::fprintResult(int index, const Eigen::MatrixXd& result)
{
	std::ofstream file2("result.txt", std::ios::app); // 打开文件以进行写入
	if (file2.is_open()) {
		file2 << "................ " << index << " ................" << std::endl;
		file2 << result << std::endl; // 使用Eigen的输出运算符
		file2.close(); // 关闭文件
		//std::cout << "Matrix saved to matrix.txt" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing" << std::endl;
	}
}

double UtprMpc::si(double x) {
	return (x == 0) ? 1.0 : std::sin(x) / x;
}