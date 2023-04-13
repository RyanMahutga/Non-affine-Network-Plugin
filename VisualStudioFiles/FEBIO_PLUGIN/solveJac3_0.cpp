/* solveJac3_0.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves for the network stresses given a deformation.

CALLED FROM: jacobianSolver.cpp

CALLS ON: networkSolver.cpp

INPUTS:

OUTPUTS: None

NOTE:

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 03-04-2022
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "math.h"
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

const double R = 8.314; // universal gas constant
const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void solveJac3_0(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, Vector3d& err, double P, 
	double fib_vol_m)
{
	// solving the network
	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, guess, stab_node);

	// reading the !network! stress
	net_stress = netResults_n.net_stress;

	// parameters for pressure term
	double c_star = 150, // M - bathing solution concentration
		T = 310, // temperature K
		k = 200; //k gives cfcd give 60mEq/L at 30% actin

	double J = F.determinant(); // volume of network

	double phim = fib_vol_m / (x_scale*x_scale*J); // actin volume fraction

	double c_fcd = k * phim; // fixed charge density for actin

	double pressure = R * T * (sqrt(c_fcd*c_fcd + 4 * c_star*c_star) - 2 * c_star); // calculating pressure
	netResults_n.pressure = pressure; // storing pressure

	// calculating !total! stress
	net_stress(0, 0) = net_stress(0, 0) - pressure;
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	// residual for zero applied stress
	err(0) = (appliedStress(0)-net_stress(0, 0));
	err(1) = (appliedStress(0)-net_stress(1, 1));
	err(2) = (appliedStress(0)-net_stress(2, 2));
}