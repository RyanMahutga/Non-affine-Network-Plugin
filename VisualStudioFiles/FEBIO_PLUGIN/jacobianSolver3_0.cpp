/* jacobianSolver3_0.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves for the network deformation gradient that gives the stresses defined in appliedStress
(which should be zero, or very low) using a damped Newton loop with a numerical Jacobian.

CALLED FROM: Main

CALLS ON: solveJac3_0.cpp

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
#include "time.h"
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>

using namespace std;
using namespace Eigen;

void jacobianSolver3_0(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P, int type_flag)
{
	network net0 = network0; // undeformed network
	network net = network_n; // deformed network

	double  dlam = 1e-5, // dlam for numerical jacobian
		mag_err = 1, err_thresh = 10,  stored_err = 1000; // error magnitude, threshhol, and previous iteration stored error

	Matrix3d dF, F_o; // change in F, original F
		F_o = F; // original F
	MatrixXd J(3,3); // Jacobian initialize
	Vector3d epsilon, errors, // change applied to F, errors (residual of stress)
		init_err, err; // initial error (prior to numerical jacobian), err error in numerical jacobian 

	// setting error threshhold for convergence
	err_thresh = abs(1e-3 * appliedStress.maxCoeff());
	if (err_thresh < 1)
	{
		err_thresh = 1;
	}

	srand(time(NULL)); // seed random

	for (int n = 0; n < 100; n++)
	{
		// clear Jacobian
		J.setZero();

		if (n == 1) // set F to original F
		{
			F(0, 0) =  F_o(0, 0);
			F(1, 1) =  F_o(1, 1);
			F(2, 2) =  F_o(2, 2);
		}

		// initial solve
		solveJac3_0(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
			guess, stab_node, appliedStress, net_stress, init_err, P, fib_vol_m);

		// error magnitude
		mag_err = sqrt(net_stress(0, 0) * net_stress(0, 0) + net_stress(1, 1) * net_stress(1, 1) + net_stress(2, 2) * net_stress(2, 2));

		if (mag_err < err_thresh) // convergence check
		{
			netResults_n.F = F;
			//cout << "Converged! Error: " << mag_err << endl << endl;
			break;
		}

			for (int i = 0; i < 3; i++) // numerical Jacobian calculation
			{
				// setting dF - F moved by dlam in each dimension
				dF = F; 

				if (i == 0)
				{
					dF(0, 0) = F(0, 0) + (dlam);
				}
				else if (i == 1)
				{
					dF(1, 1) = F(1, 1) + (dlam);
				}
				else if (i == 2)
				{
					dF(2, 2) = F(2, 2) + (dlam);
				}

				solveJac3_0(num_nodes, num_fibers, x_scale, dF, network0,
					network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress,
					net_stress, err, P, fib_vol_m);

				// J is change in error/dlam
				J(0, i) = (err(0) - init_err(0)) / (dlam); // 00
				J(1, i) = (err(1) - init_err(1)) / (dlam); // 01
				J(2, i) = (err(2) - init_err(2)) / (dlam); // 02
			}

			// residuals are errors
			errors(0) = -init_err(0);
			errors(1) = -init_err(1);
			errors(2) = -init_err(2);

			epsilon = J.bdcSvd(ComputeThinU | ComputeThinV).solve(errors); // epsilon is change to be applied to F0

			stored_err = mag_err; // storing error

		// setting upper and lower bounds for how far the F can change on each iteration
		for (int j = 0; j < 3; j++)
		{
			if (epsilon(j) > 0.01)
			{
				epsilon(j) = 0.01;
			}
			else if (epsilon(j) < -0.01)
			{
				epsilon(j) = -0.01;
			}
		}

		cout << F << endl << errors << endl<<netResults_n.pressure << endl<< endl;

		// damping the Newton loop
		if (n < 79) // initially
		{
			if (n < 10) // damped 0.1 initially
			{
				F(0, 0) += 0.5 * epsilon(0);
				F(1, 1) += 0.5 * epsilon(1);
				F(2, 2) += 0.5 * epsilon(2);

			}
			else // damped 0.5 maximum
			{
				F(0, 0) += epsilon(0);
				F(1, 1) += epsilon(1);
				F(2, 2) += epsilon(2);
			}
		}
		else // for many iterations higher damping to help with oscillations
		{
			if (n < 89)
			{
				F(0, 0) += 0.05 * epsilon(0);
				F(1, 1) += 0.05 * epsilon(1);
				F(2, 2) += 0.05 * epsilon(2);
			}
			else
			{
				F(0, 0) += 0.01 * epsilon(0);
				F(1, 1) += 0.01 * epsilon(1);
				F(2, 2) += 0.01 * epsilon(2);
			}
		}

		netResults_n.F = F; // setting new F in netResults
	}
}