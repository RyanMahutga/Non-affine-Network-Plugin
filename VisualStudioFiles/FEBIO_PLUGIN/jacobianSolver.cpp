/* jacobianSolver.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves for the network deformation gradient that gives the stresses defined in appliedStress 
using a damped Newton loop with a numerical Jacobian.  

CALLED FROM: Main

CALLS ON: solveJac.cpp

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

void jacobianSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess,
	VectorXi stab_node, VectorXd appliedStress, Matrix3d& net_stress, double& pressure, Vector3d x_guess,
	double fib_vol_m, double phi_m, double P, int type_flag)
{
	network net0 = network0; // network undeformed
	network net = network_n; // network deformed

	double dlam = 1e-5, // dlam for numerical jacobian 
		mag_err = 1, err_thresh = 1e-3, stored_err = 10, // magnitude of error, error threshhold, stored previous iteration error
		scalef = 1.0; // newton damping coefficient

	Matrix3d dF, F_o; // dF - F for numerical jacobian, F_o is original F
	MatrixXd J(9, 9); // jacobian
	VectorXd init_err(6), err(6); // stress terms (6 from Cauchy)
	VectorXd epsilon(9), errors(9); // epsilon (change to F), errors (residual of stress)

	// set convergence threshold
	err_thresh = abs(1e-4 * appliedStress.maxCoeff());
	if (err_thresh < 1)
	{
		err_thresh = 1;
	}

	if (type_flag == 1)
	{
		dlam = 5e-5;
	}

	F_o = F; // original F

	srand(time(NULL)); // seed random

	for (int n = 0; n < 150; n++)
	{
		// clear vars
		J.setZero();

		// initial solve
		solveJac(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
			guess, stab_node, appliedStress, net_stress, init_err, P, fib_vol_m);

		mag_err = init_err.norm(); // error 

		stored_err = 1.1 * mag_err;

		if (mag_err < err_thresh) // check convergence
		{
			netResults_n.F = F;
			//cout << "Converged! Error: " << mag_err << endl << endl;
			break;
		}
		for (int i = 0; i < 9; i++) // loop numerical Jacobian
				{
					dF = F; // reset F

					// add dlam to F terms
					if (i==0)
					{ 
						dF(0, 0) = F(0, 0) + (dlam);
					}
					else if (i == 1)
					{
						dF(0, 1) = F(0, 1) + (dlam);
					}
					else if (i==2)
					{ 
						dF(0, 2) = F(0, 2) + (dlam);
					}
					else if (i==3)
					{ 
						dF(1, 0) = F(1, 0) + (dlam);
					}
					else if (i==4)
					{ 
						dF(1, 1) = F(1, 1) + (dlam);
					}
					else if (i == 5)
					{
						dF(1, 2) = F(1, 2) + (dlam);
					}
					else if (i == 6)
					{
						dF(2, 0) = F(2, 0) + (dlam);
					}
					else if (i == 7)
					{
						dF(2, 1) = F(2, 1) + (dlam);
					}
					else if (i == 8)
					{
						dF(2, 2) = F(2, 2) + (dlam);
					}

					solveJac(num_nodes, num_fibers, x_scale, dF, network0,
						network_n, fiberResults_n, netResults_n, guess, stab_node, appliedStress,
						net_stress, err, P, fib_vol_m);

					// numerical jacobian
					J(0, i) = (err(0) - init_err(0)) / (dlam); // 00
					J(1, i) = (err(1) - init_err(1)) / (dlam); // 01
					J(2, i) = (err(2) - init_err(2)) / (dlam); // 02
					J(3, i) = (err(1) - init_err(1)) / (dlam); //10
					J(4, i) = (err(3) - init_err(3)) / (dlam); //11
					J(5, i) = (err(4) - init_err(4)) / (dlam);//12
					J(6, i) = (err(2) - init_err(2)) / (dlam);//20
					J(7, i) = (err(4) - init_err(4)) / (dlam);//21
					J(8, i) = (err(5) - init_err(5)) / (dlam);//22

		}
				// stress residual
				errors(0) = -init_err(0);
				errors(1) = -init_err(1);
				errors(2) = -init_err(2);
				errors(3) = -init_err(1);
				errors(4) = -init_err(3);
				errors(5) = -init_err(4);
				errors(6) = -init_err(2);
				errors(7) = -init_err(4);
				errors(8) = -init_err(5);

				epsilon = J.bdcSvd(ComputeThinU | ComputeThinV).solve(errors); // calculatin changes for F

				stored_err = mag_err; // storing error

			// set boundary on change in F
			for (int j=0; j<9; j++) // loop epsilon values
			{ 
				if (j == 0 || j == 4 || j == 8) // on the diagonal of F
				{
					if (epsilon(j) > 0.1)
					{
						epsilon(j) = 0.1;
					}
					else if (epsilon(j) < -0.1)
					{
						epsilon(j) = -0.1;
					}
				}
				else // shear terms
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
			}

			// damping the Newton loop
			if (n < 10)
			{
				scalef = 0.25;
			}
			else if (n < 20)
			{
				scalef = 0.5;
			}
			else
			{
				scalef = 1.0;
			}
	
			F(0, 0) += scalef * epsilon(0);
			F(0, 1) += scalef * epsilon(1);
			F(0, 2) += scalef * epsilon(2);
			F(1, 0) += scalef * epsilon(3);
			F(1, 1) += scalef * epsilon(4);
			F(1, 2) += scalef * epsilon(5);
			F(2, 0) += scalef * epsilon(6);
			F(2, 1) += scalef * epsilon(7);
			F(2, 2) += scalef * epsilon(8);

			// check lower bounds
			if (F(0, 0) < 0.2)
			{
				F(0, 0) = 0.2;
			}
			if (F(1, 1) < 0.2)
			{
				F(1, 1) = 0.2;
			}
			if (F(2, 2) < 0.2)
			{
				F(2, 2) = 0.2;
			}
	}
	netResults_n.F = F;
}