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

void elasticity(network network0, int num_fibs, Matrix3d F, Matrix3d net_stress, double (&C)[3][3][3][3])
{
	// Zero C_IJKL
	for (int I = 0; I < 3; I++)
	{
		for (int J = 0; J < 3; J++)
		{
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					C[I][J][K][L] = 0.0;
				}
			}
		}
	}

	Matrix3d F_inv = F.inverse(), P;
	double J = F.determinant();

	MatrixXi fibers = network0.fibers;
	MatrixXi fibtype = network0.fib_type;
	MatrixXd nodes = network0.nodes;
	MatrixXd init_lens = network0.init_lens;

	P = J * net_stress * F_inv; // calculate PK1 stress from Cauchy

	int node1, node2;
	Vector3d realNode2, vect, unit;
	double vect_len,k_fib;
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

		// fiber stiffness
		if (fibtype(e) == 1)
		{ 
			k_fib = 1e6;
		}
		else if (fibtype(e) == 2)
		{ 
			k_fib = 1e6;
		}
		else if (fibtype(e) == 3)
		{ 
			k_fib = 1e6;
		}

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
							for (int w = 0; w < 3; w++)
							{
								C[I][J][K][L] = C[I][J][K][L] + ();
							}
						}
					}
				}
			}
		}
	}



}