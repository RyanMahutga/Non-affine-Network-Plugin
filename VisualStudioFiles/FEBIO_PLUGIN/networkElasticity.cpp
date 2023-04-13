/* elasticity.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves the elasticity in a periodic network from fibers and nodes in an RVE subjected to a defined stretch.
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

void networkElasticity(Matrix3d F, Matrix3d U_inv, Matrix3d U, Matrix3d R, network network0, fiberResults fiberResults_n,  Matrix3d net_stress, double fib_vol, double (&c)[3][3][3][3])
{
	double C1[3][3][3][3], C2[3][3][3][3], CH[3][3][3][3], dU_dC[3][3][3][3], C_total[3][3][3][3], C_sum[3][3];
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

	// calculating dU_dC ------------------------------------------------------------------------------------
	
	Matrix3d Id; Id << 1, 0, 0, 0, 1, 0, 0, 0, 1; // identitity

	// approximated dU_dC
	for (int M = 0; M < 3; M++)
	{
		for (int W = 0; W < 3; W++)
		{
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					dU_dC[M][W][K][L] = 1 / 4 * (U_inv(M, W) * Id(L, K) + U_inv(M, L) * Id(K, W));
				}
			}
		}
	}
	
	/*
	// Mandel form of dC_dU - technically more correct, but in practice the above approximation is more relaible
	MatrixXd C_s(6, 6); C_s << 2 * U(0, 0), 0, 0, 0, sqrt(2) * U(0, 2), sqrt(2) * U(0, 1),
		0, 2 * U(1, 1), 0, sqrt(2) * U(1, 2), 0, sqrt(2) * U(1, 0),
		0, 0, 2 * U(2, 2), sqrt(2) * U(2, 1), sqrt(2) * U(2, 0), 0,
		0, 0, 2 * sqrt(2) * U(1, 2), 2 * U(1, 1), 2 * U(1, 0), 0,
		2 * sqrt(2) * U(2, 0), 0, 0, 0, 2 * U(2, 2), 2 * U(2, 1),
		0, 2 * sqrt(2) * U(0, 1), 0, 2 * U(0, 2), 0, 2 * U(0, 0) ;

	//inverse of the mandel form
	MatrixXd C_inv(6, 6); C_inv = C_s.inverse();

	// extracting all 36 unique values (remember C ijkl = C jikl = C ijlk = C klij
	// for storage in the 81 spaces (3x3x3x3)

	for (int i = 0; i < 3; i++)
	{
		dU_dC[i][i][0][0] = C_inv(i, 0);
		dU_dC[i][i][1][1] = C_inv(i, 1);
		dU_dC[i][i][2][2] = C_inv(i, 2);
		dU_dC[i][i][1][2] = C_inv(i, 3) / sqrt(2);
		dU_dC[i][i][2][0] = C_inv(i, 4) / sqrt(2);
		dU_dC[i][i][0][1] = C_inv(i, 5) / sqrt(2);
		// by symmetry 
		dU_dC[i][i][2][1] = C_inv(i, 3) / sqrt(2);
		dU_dC[i][i][0][2] = C_inv(i, 4) / sqrt(2);
		dU_dC[i][i][1][0] = C_inv(i, 5) / sqrt(2);
	}

	// 12mn and 21mn 
	dU_dC[1][2][0][0] = C_inv(3, 0) / sqrt(2);
	dU_dC[1][2][1][1] = C_inv(3, 1) / sqrt(2);
	dU_dC[1][2][2][2] = C_inv(3, 2) / sqrt(2);
	dU_dC[1][2][1][2] = C_inv(3, 3) / 2;
	dU_dC[1][2][2][0] = C_inv(3, 4) / 2;
	dU_dC[1][2][0][1] = C_inv(3, 5) / 2;
	// symmetric part
	dU_dC[1][2][2][1] = C_inv(3, 3) / 2;
	dU_dC[1][2][0][2] = C_inv(3, 4) / 2;
	dU_dC[1][2][1][0] = C_inv(3, 5) / 2;

	// by symmetry from above
	dU_dC[2][1][0][0] = C_inv(3, 0) / sqrt(2);
	dU_dC[2][1][1][1] = C_inv(3, 1) / sqrt(2);
	dU_dC[2][1][2][2] = C_inv(3, 2) / sqrt(2);
	dU_dC[2][1][1][2] = C_inv(3, 3) / 2;
	dU_dC[2][1][2][0] = C_inv(3, 4) / 2;
	dU_dC[2][1][0][1] = C_inv(3, 5) / 2;
	// symmetric part
	dU_dC[2][1][2][1] = C_inv(3, 3) / 2;
	dU_dC[2][1][0][2] = C_inv(3, 4) / 2;
	dU_dC[2][1][1][0] = C_inv(3, 5) / 2;

	//20mn and 02mn
	dU_dC[2][0][0][0] = C_inv(4, 0) / sqrt(2);
	dU_dC[2][0][1][1] = C_inv(4, 1) / sqrt(2);
	dU_dC[2][0][2][2] = C_inv(4, 2) / sqrt(2);
	dU_dC[2][0][1][2] = C_inv(4, 3) / 2;
	dU_dC[2][0][2][0] = C_inv(4, 4) / 2;
	dU_dC[2][0][0][1] = C_inv(4, 5) / 2;
	// symmetric part
	dU_dC[2][0][2][1] = C_inv(4, 3) / 2;
	dU_dC[2][0][0][2] = C_inv(4, 4) / 2;
	dU_dC[2][0][1][0] = C_inv(4, 5) / 2;

	// by symmetry from above
	dU_dC[0][2][0][0] = C_inv(4, 0) / sqrt(2);
	dU_dC[0][2][1][1] = C_inv(4, 1) / sqrt(2);
	dU_dC[0][2][2][2] = C_inv(4, 2) / sqrt(2);
	dU_dC[0][2][1][2] = C_inv(4, 3) / 2;
	dU_dC[0][2][2][0] = C_inv(4, 4) / 2;
	dU_dC[0][2][0][1] = C_inv(4, 5) / 2;
	// symmetric part
	dU_dC[0][2][2][1] = C_inv(4, 3) / 2;
	dU_dC[0][2][0][2] = C_inv(4, 4) / 2;
	dU_dC[0][2][1][0] = C_inv(4, 5) / 2;

	// 01mn and 10mn
	dU_dC[0][1][0][0] = C_inv(5, 0) / sqrt(2);
	dU_dC[0][1][1][1] = C_inv(5, 1) / sqrt(2);
	dU_dC[0][1][2][2] = C_inv(5, 2) / sqrt(2);
	dU_dC[0][1][1][2] = C_inv(5, 3) / 2;
	dU_dC[0][1][2][0] = C_inv(5, 4) / 2;
	dU_dC[0][1][0][1] = C_inv(5, 5) / 2;
	// symmetric part
	dU_dC[0][1][2][1] = C_inv(5, 3) / 2;
	dU_dC[0][1][0][2] = C_inv(5, 4) / 2;
	dU_dC[0][1][1][0] = C_inv(5, 5) / 2;

	// by symmetry from above
	dU_dC[1][0][0][0] = C_inv(5, 0) / sqrt(2);
	dU_dC[1][0][1][1] = C_inv(5, 1) / sqrt(2);
	dU_dC[1][0][2][2] = C_inv(5, 2) / sqrt(2);
	dU_dC[1][0][1][2] = C_inv(5, 3) / 2;
	dU_dC[1][0][2][0] = C_inv(5, 4) / 2;
	dU_dC[1][0][0][1] = C_inv(5, 5) / 2;
	// symmetric part
	dU_dC[1][0][2][1] = C_inv(5, 3) / 2;
	dU_dC[1][0][0][2] = C_inv(5, 4) / 2;
	dU_dC[1][0][1][0] = C_inv(5, 5) / 2;
	*/

	// calculating network part ------------------------------------------------------------------------------

	double scale = network0.fiber_vol_fract,
		detF = F.determinant(); // volume of network
	/* 
	// calculating osmotic pressure term
	double c_star = 150.0, // M - bathing solution concentration
		T = 310.0, // temperature K
		k_fcd = 0.0, //k gives cfcd give 60mEq/L at 30% actin
		phim = fib_vol / (scale * scale * detF); // fiber fraction of actin

	const double R = 8.314; // universal gas constant
	double c_fcd = k_fcd * phim; // fixed charge density of actin

	double pressure = R * T * (sqrt(c_fcd * c_fcd + 4 * c_star * c_star) - 2 * c_star);// hydrostatic pressure
	*/

	Matrix3d F_inv = F.inverse(), F_invT = F_inv.transpose(), P;

	// initialize fiber terms for solution of elasticity
	int num_fibs = network0.num_fibers; // number of fibers
	MatrixXd fibers = network0.fibers; // fibers connectivity
	MatrixXi fibtype = network0.fib_type; // fiber type
	MatrixXd nodes = network0.nodes; // nodal coordinated
	MatrixXd init_lens = network0.init_lens; // fiber undeformed lengths
	MatrixXd fib_rads = network0.fib_rads; // fiber radii
	MatrixXd fib_areas = network0.fib_areas; // fiber cross-sectionala areas

	VectorXd fib_forces = fiberResults_n.fib_forces, fib_stiffness = fiberResults_n.fib_stiffness; // fiber force and stiffness

	int node1, node2; // node number integers
	Vector3d realNode2, vect, vectF, unit, unitF, node_force1, node_force2; // vectors for nodal coordinate, fiber vectors, fiber unit vectors, and node forces
	double vect_len, k_fib, s2 = 1 / (scale * scale), // 1/scale^2 for proper units due to scaling of init_lens/V_0
		fib_force, dFdlam, lambda;

	VectorXd fib_stretch = fiberResults_n.fib_stretch; // fiber stretch

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

		unit = 1 / vect_len * vect; // fiber unit vector (undeformed)

		vectF = F * vect; unitF = vectF / vectF.norm(); // deformed unit vector

		lambda = vectF.norm() / init_lens(e); // fiber stretch

		// fiber stiffness
		k_fib = fib_stiffness(e);

		// sum over all fibers (indices are I & K with dummy (sum) index M)
		for (int I = 0; I < 3; I++) // index I 
		{
			for (int K = 0; K < 3; K++) // index K
			{
				for (int M = 0; M < 3; M++) // dummy sum index M
				{
					C_sum[I][K] = C_sum[I][K] + (k_fib * (unit(K) - unit(M) * U_inv(M, K)) - F_inv(K, M) * fib_forces(e) * unitF(M)) 
						* init_lens(e)* unit(I);
				}
			}
		}
	}

	// calculating total material elasticity tensor ----------------------------------------------------------------

	for (int I = 0; I < 3; I++)
	{
		for (int J = 0; J < 3; J++)
		{
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					for (int M = 1; M < 3; M++)
					{
						for (int Q = 0; Q < 3; Q++)
						{
							C_total[I][J][K][L] = C_total[I][J][K][L] + s2 * C_sum[I][Q] * U_inv(J, M) * dU_dC[M][Q][K][L]; // -CH[I][J][K][L]; // full material elasticity tensor
						}
					}

				}
			}
		}
	}
	// horrible ugly loop to get spatial elasticity (probably an easier way, --------------------------------------------
	// but there are too many dimensions for me to think through) -------------------------------------------------------
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