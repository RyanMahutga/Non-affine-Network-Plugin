/* calcForcesPeriodic.cpp ----------------------------------------------------------

SUMMARY: This function calculates the fiber and nodal forces in the network given 
the current nodal positions. This function is used within the Newton loop to 
equilibrate the network, and to calculate the stress.

CALLED FROM: networkSolver.cpp 
CALLS ON: fiberConstEqn.cpp

INPUTS: nodes_n - Nx3 matrix of nodal coordinates
		nodes_mapped - Nx3 matrix of nodal coordinates mapped to undeformed
		fibers - Mx5 fiber matrix relating fiber end nodes to one another

OUTPUTS: None

NOTE: nodes_n, fibers, fib_forces, node_forces, and bnd_nodes are modified by this code.

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 03-04-2022

*/

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"
#include "dataStructureDefinitions.h"

using namespace std;
using namespace Eigen;

void calcForcesPeriodic(MatrixXd nodes_n, MatrixXd nodes_mapped, MatrixXd fibers_n, network& network_n,
	Matrix3d F, fiberResults& fiberResults_n, netResults& netResults_n, int num_fibers, int num_nodes)
{

	VectorXd init_lens = network_n.init_lens; // fiber legnth
	VectorXi fib_type = network_n.fib_type; // fiber type
	VectorXd fib_areas = network_n.fib_areas; // fiber cross-sectionala rea
	VectorXd fib_stretch(num_fibers), // fiber stretch 
		fib_stress(num_fibers), // fiber stress
		fib_forces(num_fibers), // fiber force
		fib_stiffness(num_fibers); // fiber tangent stiffness

	VectorXd node_forces(3*num_nodes); node_forces = VectorXd::Zero(3*num_nodes); // node forces

	int node1 = 0, // node1 number
		node2 = 0; // node2 number
	double vect_len = 0, // fiber vector length 
			lambda = 0, // fiber stretch
			fib_force = 0, // fiber force
			dFdlam = 0; // d fiber force/ d fiber stretch

	Vector3d node_force1(3); node_force1 = Vector3d::Zero(3); // node 1 forces
	Vector3d node_force2(3); node_force2 = Vector3d::Zero(3); // node 2 forces
	Vector3d vect(3); vect = Vector3d::Zero(3); // fiber vector
	Vector3d realNode2(3); realNode2 = Vector3d::Zero(3); // real node 2 position
	Vector3d unit(3); unit = Vector3d::Zero(3); // unit vector from node1 to node2

	for (int n = 0; n<num_fibers; n++) // loop through fibers to calculate forces
	{
		node1 = (int)(fibers_n(n, 0) + 0.001); // first (static) node of fiber 
		node2 = (int)(fibers_n(n, 1) + 0.001); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,2:4))

		// Shift node2 to real position (nearest connected neighbor) to node1
		realNode2(0) = nodes_n(node2, 0) + fibers_n(n, 2)*F(0, 0) + fibers_n(n, 3)*F(0, 1) + fibers_n(n, 4)*F(0, 2);
		realNode2(1) = nodes_n(node2, 1) + fibers_n(n, 2)*F(1, 0) + fibers_n(n, 3)*F(1, 1) + fibers_n(n, 4)*F(1, 2);
		realNode2(2) = nodes_n(node2, 2) + fibers_n(n, 2)*F(2, 0) + fibers_n(n, 3)*F(2, 1) + fibers_n(n, 4)*F(2, 2);

		vect(0) = realNode2(0) - nodes_n(node1, 0); // fiber x span
		vect(1) = realNode2(1) - nodes_n(node1, 1); // fiber x span
		vect(2) = realNode2(2) - nodes_n(node1, 2); // fiber x span

		vect_len = vect.norm(); // fiber length

		unit = 1 / vect_len * vect; // fiber unit vector

		// fiber constitutive equation
		fiberConstEqn(vect, node_force1, node_force2, lambda, init_lens(n), fib_type(n), fib_areas(n), fib_force, dFdlam);

		fib_stretch(n) = lambda; // fiber stretch
		fib_forces(n) = fib_force; // store fiber force
		fib_stress(n) = fib_force / fib_areas(n); // fiber PK1 stress
		fib_stiffness(n) = dFdlam;

		// storing nodal forces
		node_forces(3 * node1 + 0) = node_forces(3 * node1 + 0) + node_force1(0);
		node_forces(3 * node1 + 1) = node_forces(3 * node1 + 1) + node_force1(1);
		node_forces(3 * node1 + 2) = node_forces(3 * node1 + 2) + node_force1(2);

		node_forces(3 * node2 + 0) = node_forces(3 * node2 + 0) + node_force2(0);
		node_forces(3 * node2 + 1) = node_forces(3 * node2 + 1) + node_force2(1);
		node_forces(3 * node2 + 2) = node_forces(3 * node2 + 2) + node_force2(2);

	}
	// store fiber results
	fiberResults_n.fib_forces = fib_forces; 
	fiberResults_n.fib_stiffness = fib_stiffness;
	fiberResults_n.fib_stretch = fib_stretch;
	fiberResults_n.fib_stress = fib_stress;
	netResults_n.node_forces = node_forces;
}