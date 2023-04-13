/* calcForcePeriodic.h is the header file for solving the stretch control problem for periodic RVEs

Created by : Ryan Mahutga - Alford Lab - University of Minnesota
Date Modified : 04-10-23
*/

#pragma once
#include "Eigen/Eigen"
#include "dataStructureDefinitions.h"

using namespace Eigen;

void stressCalc(Matrix3d F, int netnum, double netSave, double netSolve, network& network0, fiberResults& fiberResults_n, Matrix3d& net_stress);

void networkSolver(int num_nodes, int num_fibers, double x_scale, Matrix3d F, network& network0,
	network& network_n, fiberResults& fiberResults_n, netResults& netResults_n, bool guess, VectorXi stab_node, double netSolve);
// solve network

void updateCrossing(MatrixXd& nodes_mapped, MatrixXd& nodes_n, MatrixXd& fibers_n, int num_nodes, int num_fibers, Matrix3d F);
// updating crossing conditions if fibers leave RVE boundary

void calcForcesPeriodic(MatrixXd nodes_n, MatrixXd nodes_mapped, MatrixXd fibers_n, network& network_n,
	Matrix3d F, fiberResults& fiberResults_n, netResults& netResults_n, int num_fibers, int num_nodes); // calc forces

void fiberConstEqn(Vector3d vect1, Vector3d& node_force1, Vector3d& node_force2, double& lambda, double init_len,
	int fibtype, double fib_area, double& fib_force, double& dFdlam); // fiber Constitutive models

void calcJacobian2(MatrixXd& nodes_n, MatrixXd& fibers_n, network& network_n, Matrix3d F, VectorXd& node_forces,
	int num_fibers, int num_nodes, SparseMatrix<double>& J, VectorXi stab_node); // Jacobain and Residual calculation

void componentStress(MatrixXd nodes_n, MatrixXd fibers_n, network& network_n, netResults& netResults_n, 
	fiberResults& fiberResults_n, Matrix3d F, int num_fibs, double x_scale); 
// component stress calculation

void networkElasticity(Matrix3d F, Matrix3d U_inv, Matrix3d U, Matrix3d R, network network0, fiberResults fiberResults_n, Matrix3d net_stress, double fib_vol, double(&c)[3][3][3][3]);

void appendResults(Matrix3d F, MatrixXd component_stress, Matrix3d net_stress, VectorXd fib_stretch,
	double pressure, int time);

void overwriteResults(Matrix3d F, MatrixXd component_stress, Matrix3d net_stress, VectorXd fib_stretch,
	VectorXd fib_stress, int time);
