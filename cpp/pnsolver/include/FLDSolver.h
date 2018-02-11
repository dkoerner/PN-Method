#pragma once
#include <Eigen/Sparse>
#include <PNVolume.h>





struct FLDSolver
{

	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve(PNVolume &volume, double tol, int maxIterations);
	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve2(PNVolume &volume, double tol, int maxIterations);

};
