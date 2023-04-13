/* solveJac.cpp -------------------------------------------------------------------------------

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

void solveJac(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, VectorXd& err, double P, 
	double fib_vol_m)
{
	// solve network
	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, guess, stab_node);

	net_stress = netResults_n.net_stress; // pull network stress

	double c_star = 150, // M - bathing solution concentration
		T = 310, // temperature K
		k = 200, //k gives cfcd give 60mEq/L at 30% actin
		J = F.determinant(), // volume of network
		phim = fib_vol_m / (x_scale * x_scale * J); // fiber fraction of actin

	double c_fcd = k * phim; // fixed charge density of actin

	double pressure = R * T * (sqrt(c_fcd * c_fcd + 4 * c_star * c_star) - 2 * c_star);// hydrostatic pressure
	netResults_n.pressure = pressure; // store pressure

	// total stress
	net_stress(0, 0) = net_stress(0, 0) - pressure; 
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	// error between applied stress and net stress
	err(0) = (appliedStress(0) - net_stress(0, 0));
	err(1) = (appliedStress(1) - net_stress(0, 1));
	err(2) = (appliedStress(2) - net_stress(0, 2));
	err(3) = (appliedStress(3) - net_stress(1, 1));
	err(4) = (appliedStress(4) - net_stress(1, 2));
	err(5) = (appliedStress(5) - net_stress(2, 2));
}