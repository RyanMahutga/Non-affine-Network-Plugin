/* elasticity.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves the force balance on nodes in a periodic fiber RVE subjected to a user defined stretch.
The solution can use any constitutive equations for fibers (defined in fiberConstEqn.cpp), the network equilibrium
is found through a damped Netwon's Method relying on a Conjugate Gradient solver in the Eigen Library.

CALLED FROM: solveJac.cpp & solveJac3_0.cpp

CALLS ON: calcForcesPeriodic.cpp, calcJacobian.cpp, componentStress2.cpp

INPUTS:

OUTPUTS: None

NOTE: This code modifies 'nodes', bnd_nodes','fibers', 'net_stress', 'fib_forces', 'fib_stress', and 'rve_stretch'

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 03-04-2022
*/

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>
#include <ctime>
#include "dataStructureDefinitions.h"

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void elasticity(Matrix3d F, network network0, fiberResults fiberResults_n,  Matrix3d net_stress, double fib_vol, double (&c)[3][3][3][3])
{
	double C1[3][3][3][3], C2[3][3][3][3], CH[3][3][3][3], C_total[3][3][3][3], C_sum[3][3];
	// Zero C_IJKL
	for (int I = 0; I < 3; I++)
	{
		for (int J = 0; J < 3; J++)
		{
			C_sum[I][J] = 0.0;
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					C1[I][J][K][L] = 0.0;
					C2[I][J][K][L] = 0.0;
					CH[I][J][K][L] = 0.0;
					C_total[I][J][K][L] = 0.0;
					c[I][J][K][L] = 0.0;
				}
			}
		}
	}

	double scale = network0.fiber_vol_fract,
		detF = F.determinant(); // volume of network
	/*
	double c_star = 150.0, // M - bathing solution concentration
		T = 310.0, // temperature K
		k_fcd = 0.0, //k gives cfcd give 60mEq/L at 30% actin
		phim = fib_vol / (scale * scale * detF); // fiber fraction of actin

	const double R = 8.314; // universal gas constant
	double c_fcd = k_fcd * phim; // fixed charge density of actin

	double pressure = R * T * (sqrt(c_fcd * c_fcd + 4 * c_star * c_star) - 2 * c_star);// hydrostatic pressure
	*/
	int num_fibs = network0.num_fibers;

	Matrix3d F_inv = F.inverse(), F_invT = F_inv.transpose(), P;
	MatrixXd fibers = network0.fibers;
	MatrixXi fibtype = network0.fib_type;
	MatrixXd nodes = network0.nodes;
	MatrixXd init_lens = network0.init_lens;
	MatrixXd fib_rads = network0.fib_rads;
	MatrixXd fib_areas = network0.fib_areas;

	VectorXd fib_forces = fiberResults_n.fib_forces;

	P = detF * net_stress * F_inv.transpose(); // calculate PK1 stress from Cauchy

	int node1, node2;
	Vector3d realNode2, vect, vectF, unit, unitF, node_force1, node_force2;
	double vect_len,k_fib,s2 = 1/(scale*scale), fib_force, dFdlam, lambda;
	// Populate C_IJKL
	for (int e = 0; e < num_fibs; e++) // loop fibers
	{
		// calculate fiber direction vector
		node1 = (int)(fibers(e, 0) + 0.001); // first (static) node of fiber 
		node2 = (int)(fibers(e, 1) + 0.001); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,2:4))

		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes(node2, 0) + fibers(e, 2);
		realNode2(1) = nodes(node2, 1) + fibers(e, 3);
		realNode2(2) = nodes(node2, 2) + fibers(e, 4);

		vect(0) = realNode2(0) - nodes(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		vectF = F * vect;

		unitF = vectF / vectF.norm();

		fiberConstEqn(vectF, node_force1, node_force2, lambda, init_lens(e),
			fibtype(e), fib_areas(e), fib_force, dFdlam);

		//Force = F*
		// fiber stiffness
		if (fibtype(e) == 1)
		{
			k_fib = 1e6 * fib_areas(e);
		}
		else if (fibtype(e) == 2)
		{
			k_fib = 1e6 * fib_areas(e);;
		}
		else if (fibtype(e) == 3)
		{
			k_fib = 1e6 * fib_areas(e);
		}
		if (fib_forces(e) < 0.0)
		{
			k_fib = k_fib / 1000.0;
		}
		//cout << "Fiber " << e << ". Force: " << fib_forces(e) << endl;

		for (int I = 0; I < 3; I++)
		{
			for (int L = 0; L < 3; L++)
			{
				for (int j = 0; j < 3; j++)
				{
					C_sum[I][L] = C_sum[I][L] + k_fib * unit(L) * init_lens(e) * unit(I)
						- F_inv(L, j) * fib_forces(e) * unitF(j) * init_lens(e) * unit(I);
				}
			}
		}
	}

		/*
		for (int I = 0; I < 3; I++)
		{
			for (int J = 0; J < 3; J++)
			{
				for (int K = 0; K < 3; K++)
				{
					for (int L = 0; L < 3; L++)
					{
						// intermediate summations of quantities
						for (int Q = 0; Q < 3; Q++)
						{
							C1[I][J][K][L] = C1[I][J][K][L] + s2 * k_fib * unit(L) * init_lens(e) * unit(I) * F_inv(J, Q) * F_invT(K, Q);
							for (int w = 0; w < 3; w++)
							{
								C2[I][J][K][L] = C2[I][J][K][L] + F_inv(L,w) * P(w,I) * F_inv(J, Q) * F_invT(K,Q);
							}
						}
					}
				}
			}
		}*/

	for (int I = 0; I < 3; I++)
	{
		for (int L = 0; L < 3; L++)
		{
		for (int J = 0; J < 3; J++)
		{
			for (int K = 0; K < 3; K++)
			{
					for (int p = 0; p < 3; p++)
					{
						C_total[I][J][K][L] = s2 * F_inv(K, p) * F_inv(J, p) * C_sum[I][L]; // -CH[I][J][K][L]; // full material elasticity tensor
					}
				}
			}
		}
	}

	// horrible ugly loop to get spatial elasticity (probably an easier way, 
	// but there are too many dimensions for me to think through)
	for (int I = 0; I < 3; I++)
	{
		for (int J = 0; J < 3; J++)
		{
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							for (int k = 0; k < 3; k++)
							{
								for (int l = 0; l < 3; l++)
								{
									c[i][j][k][l] = c[i][j][k][l] + F(i, I) * F(j, J) * F(k, K) * F(l, L) * C_total[I][J][K][L]  ;
								}
							}
						}
					}
				}
			}
		}
	}


}