/* networkSolver.cpp -------------------------------------------------------------------------------

SUMMARY: This function solves the force balance on nodes in a periodic fiber RVE subjected to a user defined stretch.
The solution can use any constitutive equations for fibers (defined in fiberConstEqn.cpp), the network equilibrium 
is found through a damped Netwon's Method relying on a Conjugate Gradient solver in the Eigen Library.

CALLED FROM: solveJac.cpp & solveJac3_0.cpp

CALLS ON: calcForcesPeriodic.cpp, calcJacobian.cpp, componentStress2.cpp

INPUTS: 

OUTPUTS: None

NOTE: This code modifies 'nodes_n', bnd_nodes','fibers_n', 'net_stress', 'fib_forces', 'fib_stress', and 'rve_stretch'

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

void networkSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess, 
	VectorXi stab_node, double netSolve)
{
	// check for valid stretches (looks for jacobian approaching zero or negative)
	if (F(0, 0) < 0.001 || F(1, 1) < 0.001 || F(2, 2) < 0.001)
	{
		cout << "Deformation Flattening Network! " << endl << endl;
	}

	int iter_max = 300; // maximum Newton Iterations 

	VectorXd residual(3 * num_nodes); residual = VectorXd::Zero(3 * num_nodes); // residual nodal forces
	VectorXd node_forces(3 * num_nodes); node_forces = VectorXd::Zero(3 * num_nodes); // nodal forces 

	MatrixXd nodes_n = network_n.nodes;// deformed nodes
	MatrixXd nodes = network0.nodes; // undeformed nodes

	MatrixXd fibers = network0.fibers; // undeformed fibers
	MatrixXd fibers_n = network_n.fibers; // deformed fibers

	MatrixXd nodes_n0 = network_n.nodes;// initial deformed nodes
	MatrixXd nodes0 = network0.nodes; //initial undeformed nodes

	MatrixXd fibers0 = network0.fibers; // initial undeformed fibers
	MatrixXd fibers_n0 = network_n.fibers; // initial deformed fibers

	int flag = 1, count = 0;

	double err = 1e-6, // error
		err_tol = 1e-15, // error tolerance
		damp_fact = 1.0; // Newton damping coefficient

	MatrixXd vect(num_fibers, 3); vect = MatrixXd::Zero(num_fibers, 3); // fiber vector from node1 to real position of node2
	VectorXd lambda(num_fibers); lambda = VectorXd::Zero(num_fibers); // fiber stretch
	VectorXd dFdlam(num_fibers); dFdlam = VectorXd::Zero(num_fibers); // derivative of force with respect to lambda (analytically calculated)
	VectorXd epsilon(3 * num_nodes); epsilon = VectorXd::Zero(3 * num_nodes);

	MatrixXd nodes_mapped(num_nodes, 3); nodes_mapped = nodes; // nodal positions mapped back to undeformed unit cube (update every iteration)

	Matrix3d F_inv; // inverse of def grad

	F_inv = F.inverse(); // inverse deformation mapping

	SparseMatrix<double> J(3 * num_nodes, 3 * num_nodes); //J = MatrixXd::Zero(3 * num_nodes, 3 * num_nodes); // initialize Jacobian

	for (int k = 0; k < num_nodes; k++) // initial guesses for deformed nodal coordinates
	{
			nodes_n(k, 0) = F(0, 0) * (nodes(k, 0)) + F(0, 1) * (nodes(k, 1)) + F(0, 2) * (nodes(k, 2));
			nodes_n(k, 1) = F(1, 0) * (nodes(k, 0)) + F(1, 1) * (nodes(k, 1)) + F(1, 2) * (nodes(k, 2));
			nodes_n(k, 2) = F(2, 0) * (nodes(k, 0)) + F(2, 1) * (nodes(k, 1)) + F(2, 2) * (nodes(k, 2));

		// setting mapped nodes
		nodes_mapped(k, 0) = F_inv(0, 0) * nodes_n(k, 0) + F_inv(0, 1) * nodes_n(k, 1) + F_inv(0, 2) * nodes_n(k, 2);
		nodes_mapped(k, 1) = F_inv(1, 0) * nodes_n(k, 0) + F_inv(1, 1) * nodes_n(k, 1) + F_inv(1, 2) * nodes_n(k, 2);
		nodes_mapped(k, 2) = F_inv(2, 0) * nodes_n(k, 0) + F_inv(2, 1) * nodes_n(k, 1) + F_inv(2, 2) * nodes_n(k, 2);
	}

	// updating fiber crossings
	updateCrossing(nodes_mapped, nodes_n, fibers_n, num_nodes, num_fibers, F); // update crossing conditions

	network_n.nodes = nodes_n; // reset nodes
	network_n.fibers = fibers_n; // reset fibers

	network0.nodes = nodes_mapped; // reset nodes0
	network0.fibers = fibers_n; // reset fibers0 

	if (netSolve > 0.99)
	{
		for (int iter = 0; iter < iter_max; iter++) // Newton's Method
		{
			J.setZero(); // set all values of Jacobian to zero
			node_forces.setZero(); // set node forces to zero

			// Jacobain and Residual Calculation
			calcJacobian2(nodes_n, fibers_n, network_n, F, node_forces, num_fibers, num_nodes, J, stab_node); // Jacobain calculation muust be done at least once to store Jacobian

			residual = -1.0 * node_forces; // residuals

			err = residual.norm(); // error in the sum of nodal forces

			if (err < err_tol)
			{
				break; // break out of for loop if residual error drops below tolerance
			}
			else if (iter == iter_max - 1)
			{
				cout << endl << "Network did not Converge!" << endl ;
			}

			// clearing epsilon
			epsilon = VectorXd::Zero(3 * num_nodes); // zero out epsilon (shift vector)

			// solving using Conjugate Gradient Solver (this should be in the ballpark of the matlab solution)
			ConjugateGradient<SparseMatrix<double>, Lower | Upper, DiagonalPreconditioner<double>>  solver;
			//DGMRES<SparseMatrix<double>, Lower | Upper, DiagonalPreconditioner IncompleteCholesky<double>>  solver;
			//BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>>  solver;
			solver.setTolerance(0.001 * err);
			solver.setMaxIterations(3 * num_nodes);
			solver.compute(J);
			epsilon = solver.solve(residual);

			// damping Newton Loop
			if (iter < 5 || iter > 251)
			{
				damp_fact = 0.25;
			}
			else if (iter < 10 || iter > 201)
			{
				damp_fact = 0.5;
			}
			else
			{
				damp_fact = 1.0;
			}

			// updating nodal positions in deformed and mapped
			for (int k = 0; k < num_nodes; k++)
			{
				nodes_n(k, 0) = nodes_n(k, 0) + damp_fact * epsilon(3 * k + 0);
				nodes_n(k, 1) = nodes_n(k, 1) + damp_fact * epsilon(3 * k + 1);
				nodes_n(k, 2) = nodes_n(k, 2) + damp_fact * epsilon(3 * k + 2);
			}

			// store nodes and fibers
			network_n.nodes = nodes_n;
			network_n.fibers = fibers_n;
		}

		for (int k = 0; k < num_nodes; k++)
		{
			// writing in the nodes mapped back to the undeformed unit cube
			nodes_mapped(k, 0) = F_inv(0, 0) * nodes_n(k, 0) + F_inv(0, 1) * nodes_n(k, 1) + F_inv(0, 2) * nodes_n(k, 2);
			nodes_mapped(k, 1) = F_inv(1, 0) * nodes_n(k, 0) + F_inv(1, 1) * nodes_n(k, 1) + F_inv(1, 2) * nodes_n(k, 2);
			nodes_mapped(k, 2) = F_inv(2, 0) * nodes_n(k, 0) + F_inv(2, 1) * nodes_n(k, 1) + F_inv(2, 2) * nodes_n(k, 2);
		}

		network0.nodes = nodes_mapped;
		network0.fibers = fibers_n;

		// update crossings
		updateCrossing(nodes_mapped, nodes_n, fibers_n, num_nodes, num_fibers, F); // update crossing conditions
	}

	// Final Force Calculation
	calcForcesPeriodic(nodes_n, nodes_mapped, fibers_n, network_n, F,
		fiberResults_n, netResults_n, num_fibers, num_nodes);

	// Final stress calculation
	componentStress(nodes_n, fibers_n, network_n, netResults_n, fiberResults_n, F, num_fibers, x_scale);

	netResults_n.Jacobian = MatrixXd(J); // store jacobian
}
