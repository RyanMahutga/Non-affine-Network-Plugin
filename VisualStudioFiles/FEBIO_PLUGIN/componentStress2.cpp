
/* componentStress2.cpp ----------------------------------------------------------------------------

SUMMARY: This function calculates the volume averaged stress tensor from the boundary nodes of a
periodic network along with the stresses in constituents of type 1,2,3,0.

CALLED FROM: networkSolver.cpp

CALLS ON: None

INPUTS: nodes - Nx3 array of interior node coordinates(xyz)
fibers - Mx5 matrix of nodal coordinates(columns 1 and 2) and

OUTPUTS: None

NOTE:

Created by : Ryan Mahutga - Barocas Research Group - University of Minnesota
Date Modified : 03-04-22

*/

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"

using namespace std;
using namespace Eigen;

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, netResults& netResults_n, 
	fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale )
{
	MatrixXd comp_stress(18, 3); // component stress matrix (3x3 for 4 constituents in 4 rows)

	VectorXi fib_type = network_n.fib_type; // fiber type
	VectorXd fib_forces = fiberResults_n.fib_forces; // fiber forces
	VectorXd init_lens = network_n.init_lens; // fiber initial lengths
	Vector3d unit, vect, // fiber unit vector and fiber vector
		realNode2; // real node2 coordinates

	double V = F.determinant(); // network volume

	Matrix3d sum0, sum1, sum2, sum3, sum4,sum5, sum_tot; // stress sums
	sum0 = sum1 = sum2 = sum3 = sum4 = sum5 = sum_tot = Matrix3d::Zero();

	int node1, node2; // node numbers
	double vect_len = 0.0; // fiber vector length

	for (int f = 0; f < num_fibs; f++)
	{
		node1 = (int)(fibers_n(f, 0) + 0.003);
		node2 = (int)(fibers_n(f, 1) + 0.003);

		//Generating fiber vectors
		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes_n(node2, 0) + fibers_n(f, 2)*F(0, 0) + fibers_n(f, 3)*F(0, 1) + fibers_n(f, 4)*F(0, 2);
		realNode2(1) = nodes_n(node2, 1) + fibers_n(f, 2)*F(1, 0) + fibers_n(f, 3)*F(1, 1) + fibers_n(f, 4)*F(1, 2);
		realNode2(2) = nodes_n(node2, 2) + fibers_n(f, 2)*F(2, 0) + fibers_n(f, 3)*F(2, 1) + fibers_n(f, 4)*F(2, 2);

		vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		for (int i = 0; i < 3; i++) // loop over coordinate directions
		{
			for (int j = 0; j < 3; j++) // loop over coordinate directions
			{
				if (fib_type(f) == 1) // sum type 1 fibers
				{
					sum1(i, j) = sum1(i,j) + fib_forces(f)*vect_len*unit(i)*unit(j);
				}
				else if (fib_type(f) == 2) // sum type 2 fibers
				{
					sum2(i, j) = sum2(i, j) + fib_forces(f)*vect_len*unit(i)*unit(j);
				}
				else if (fib_type(f) == 3) // sum type 3 fibers
				{
					sum3(i, j) = sum3(i, j) + fib_forces(f)*vect_len*unit(i)*unit(j);
				}
				else if (fib_type(f) == 4) // sum type 3 fibers
				{
					sum4(i, j) = sum4(i, j) + fib_forces(f) * vect_len * unit(i) * unit(j);
				}
				else if (fib_type(f) == 5) // sum type 3 fibers
				{
					sum5(i, j) = sum5(i, j) + fib_forces(f) * vect_len * unit(i) * unit(j);
				}
				else if (fib_type(f) == 0) // sum type 0 fibers
				{
					sum0(i, j) = sum0(i, j) + fib_forces(f)*vect_len*unit(i)*unit(j);
				}
			}
		}
	}
		sum_tot = sum0 + sum1 + sum2 + sum3 + sum4 + sum5; // total stresses
		netResults_n.net_stress = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum_tot; // calculate and store network stress

		// component stress
		comp_stress.block(0, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum1; // type 1 fibers
		comp_stress.block(3, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum2; // type 2 fibers
		comp_stress.block(6, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum3; // type 3 fibers
		comp_stress.block(9, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum0; // type 0 fibers
		comp_stress.block(12, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum4; // type 4 fibers
		comp_stress.block(15, 0, 3, 3) = 1.0 / V * 1.0 / x_scale * 1.0 / x_scale * sum5; // type 5 fibers

		netResults_n.component_stress = comp_stress; // store component stresses		
}
