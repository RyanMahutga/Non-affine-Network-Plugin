/* calcJacobain.cpp calculates the nuumerical JAcobian for the periodic RVE

SUMMARY: This function calculates the analytical Jacobian for the RVE of fibers
defined by fiberConstEqn.cpp. It stores the Jacobian as a sparse matrix for 
analysis using the Sparse Solvers in Eigen.

CALLED FROM: networkSolver.cpp
CALLS ON: calcForcesPeriodic.cpp

INPUTS: nodes_n - Nx3 matrix of nodal coordinates
		nodes0 - Nx3 matrix of nodal coordinates for tethering
		fibers - Mx5 fiber matrix relating fiber end nodes to one another
		init_lens - Mx1 vector of fiber initial lengths
		fib_type - Mx1 vector of fiber types
		fib_areas - Mx1 vector of fiber cross-sectional areas
		rve_stretch - 3x1 vector of RVE prinicpal stretches
		fib_forces - Mx1 vector of fiber forces (initially empty)
		node_forces - Nx3 matrix of x,y,z components of forces on nodes (initially empty)
		num_fibers - M integer number of fibers
		num_nodes - N integer number of nodes
		bnd_nodes - Px3 coordinates of boundary nodes
		J - 3*N x 3*N Jacobian matrix of dF_i/dx_j where F_i is the directional nodal force and x_j is the directional coordinate 

OUTPUTS: None

NOTE: nodes_n, fibers, fib_forces, node_forces, bnd_nodes, and J are modified by this code.

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 04-05-2018

*/

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"
#include "omp.h"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"

using namespace std;
using namespace Eigen;

void calcJacobian2(MatrixXd& nodes_n, MatrixXd& fibers_n, network& network_n, Matrix3d F, VectorXd& node_forces,
	int num_fibers, int num_nodes,  SparseMatrix<double>& J, VectorXi stab_node)
{
	VectorXd init_lens = network_n.init_lens; // fiber deformed initial lengths
	VectorXd fib_areas = network_n.fib_areas; // fiber deformed areas
	VectorXi fib_type = network_n.fib_type; // fiber types

	int xyz_num = 3 * num_nodes; // DOFs

	MatrixXd Jac = MatrixXd::Zero(xyz_num, xyz_num); //clear Dense Jacobian

	int node1, node2; // nodes numbers

	Vector3d realNode2, // nodal positions of node2
		old_vect, unit_old, // finer vector and unit vector
		node_force1, node_force2, //forces on node1 and node2 
		dlamdx1, dlamdx2; // dstretch/dcoordinate

		double fib_force, // fiber force
			dFdlam, lambda, // d forcce/ d stretch, stretch
			vect_len, // vector length
			dndx11, dndx12,  dndx21, dndx22, // d normal / dx 
			j_temp11,  j_temp12, j_temp21, j_temp22, // temporary jacobian entries
			inv_vect_len3; // 1/L^3

		for (int n = 0; n < num_fibers; n++)
		{
			node1 = (int)(fibers_n(n, 0) + 0.3); // first (static) node of fiber 
			node2 = (int)(fibers_n(n, 1) + 0.3); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,2:4))

			// real nodal coordinates of nearest node along fiber from node1 to node2
			realNode2(0) = nodes_n(node2, 0) + fibers_n(n, 2)*F(0, 0) + fibers_n(n, 3)*F(0, 1) + fibers_n(n, 4)*F(0, 2);
			realNode2(1) = nodes_n(node2, 1) + fibers_n(n, 2)*F(1, 0) + fibers_n(n, 3)*F(1, 1) + fibers_n(n, 4)*F(1, 2);
			realNode2(2) = nodes_n(node2, 2) + fibers_n(n, 2)*F(2, 0) + fibers_n(n, 3)*F(2, 1) + fibers_n(n, 4)*F(2, 2);

			// vector from node1 to node2
			old_vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
			old_vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber y span
			old_vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber z span

			vect_len = old_vect.norm(); // fiber length

			inv_vect_len3 = 1 / (vect_len*vect_len*vect_len); // 1/L^3

			// fiber unit vector
			unit_old = old_vect/vect_len;

			// calling fiber constitutive equation to get fib_forces and dFdlam
			fiberConstEqn(old_vect, node_force1, node_force2, lambda, init_lens(n), 
				fib_type(n), fib_areas(n), fib_force, dFdlam);

			// dlam/dx values for node1 x
			dlamdx1 = -1.0 / (init_lens(n)*vect_len)*old_vect;

			// dlam/dx values for node2 x
			dlamdx2 = 1.0 / (init_lens(n)*vect_len)*old_vect;

			// residual node forces
			node_forces(3 * node1 + 0) = node_forces(3 * node1 + 0) + node_force1(0);
			node_forces(3 * node1 + 1) = node_forces(3 * node1 + 1) + node_force1(1);
			node_forces(3 * node1 + 2) = node_forces(3 * node1 + 2) + node_force1(2);

			node_forces(3 * node2 + 0) = node_forces(3 * node2 + 0) + node_force2(0);
			node_forces(3 * node2 + 1) = node_forces(3 * node2 + 1) + node_force2(1);
			node_forces(3 * node2 + 2) = node_forces(3 * node2 + 2) + node_force2(2);

			// calculating Jacobian values
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if (j == i)
					{
						dndx11 = -1.0 / vect_len + old_vect(i)*old_vect(i) *inv_vect_len3; // dn/dx for F node1 and x node1
						dndx12 = 1.0 / vect_len - old_vect(i)*old_vect(i) *inv_vect_len3; // dn/dx for F node1 and x node2
						dndx21 = dndx11; // dn/dx for F node2 and x node1
						dndx22 = dndx12; // dn/dx for F node2 and x node2
					}
					else
					{
						dndx11 = old_vect(i)*old_vect(j) * inv_vect_len3; // dn/dx for F node1 and x node1
						dndx12 = -1.0* old_vect(i)*old_vect(j) * inv_vect_len3; // dn/dx for F node1 and x node2
						dndx21 = dndx11; // dn/dx for F node2 and x node1
						dndx22 = dndx12; // dn/dx for F node2 and x node2
					}

					// Calculate dFi / dxj
					j_temp11 = dFdlam * dlamdx1(j) * unit_old(i) + fib_force * dndx11; // temporary jacoian value F node1 x node1
					j_temp12 = dFdlam * dlamdx2(j) * unit_old(i) + fib_force * dndx12; // temporary jacoian value F node1 x node2
					j_temp21 = j_temp12; //-dFdlam * dlamdx1(j) * unit_old(i) - fib_force * dndx21;  // temporary jacoian value F node2 x node1
					j_temp22 = -dFdlam * dlamdx2(j) * unit_old(i) - fib_force * dndx22;  // temporary jacoian value F node2 x node2

					// this is dFi / dxi
					Jac(3 * node1 + i, 3 * node1 + j) += j_temp11;
					// this is dFj / dxj
					Jac(3 * node2 + i, 3 * node2 + j) += j_temp22;
					// this is dFi / dxj
					Jac(3 * node1 + i, 3 * node2 + j) += j_temp12;
					// this is dFj / dxi
					Jac(3 * node2 + i, 3 * node1 + j) = Jac(3 * node1 + i, 3 * node2 + j);
				}
			}
		}

	// Constraining the stabilized node to their initial positions
	double max_J = 10*Jac.maxCoeff(); // maximum J
	for (int p = 0; p < 3; p++)
	{
		if (stab_node(p) > -1)
		{
			Jac(3 * stab_node(p) + 0, 3 * stab_node(p) + 0) = max_J;
			Jac(3 * stab_node(p) + 1, 3 * stab_node(p) + 1) = max_J;
			Jac(3 * stab_node(p) + 2, 3 * stab_node(p) + 2) = max_J;

			Jac(3 * stab_node(p) + 0, 3 * stab_node(p) + 1) = 0;
			Jac(3 * stab_node(p) + 0, 3 * stab_node(p) + 2) = 0;
			Jac(3 * stab_node(p) + 1, 3 * stab_node(p) + 0) = 0;
			Jac(3 * stab_node(p) + 1, 3 * stab_node(p) + 2) = 0;
			Jac(3 * stab_node(p) + 2, 3 * stab_node(p) + 0) = 0;
			Jac(3 * stab_node(p) + 2, 3 * stab_node(p) + 1) = 0;


			node_forces(3*stab_node(p) + 0) = 0;
			node_forces(3*stab_node(p)+ 1) = 0;
			node_forces(3*stab_node(p)+ 2) = 0;
		}
	}
	
	double drop_tol = 1e-18; // drop tolerance for sparse jacobian

	J = Jac.sparseView(drop_tol, 1); // setting to sparse matrix
}