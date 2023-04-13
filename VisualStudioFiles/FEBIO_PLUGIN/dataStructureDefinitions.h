/* Data storage structures for the network solution

Ryan Mahutga - Univerity of Minnesota - 04/10/2023

Network:
	doubles:
	num_fibers - number fo fibers in network, N
	num_bnd - number of boundary nodes in network, B
	num_nodes - number of interior nodes in network, M
	fiber_vol_fract - fiber volume fraction in RVE
	
	Nx1 Vectors:
	fib_type - fiber type definition 
	fib_rads - fiber radii (in m)
	fib_areas- fiber cross sectional areas (in m^2)
	init_lens - fiber undeformed/ intitial lengths (in computational units)

	Nx5 Matrix:
	fibers - fiber definition cols 1 and 2 represent end nodes, cols 3,4,5 represent boundary crossings in x,y,z respectively 
			(ex. fibers=[0 5, 1, -1 0] connects nodes 0 and 5 with a crossings from node 0 through the positive x-face and a negative y-face)

	Mx3 Matrix
	nodes - X,Y,Z positions of each interior node (in computational units)

	Bx3 Matrix:
	bnd_nodes - X<Y<Z coordinated of each boundary node ordered as they appear in the fiber definition (in computational units)

fiberResults:
	Nx1 Vectors:
	fib_stress - fiber PK1 stress (or fiber force/fiber area) (in Pa)
	fib_stretch - fiber stretch (current length / undeformed length)
	fib_forces - fiber forces (in N)

netResults:
	double:
	pressure - pressure in the network (Pa)

	3*Mx1 Vector:
	node forces - linear vector of all node x,y,z forces

	3x3 Matrix:
	net_stress - full 3D network stress tensor
	F - network deformation gradient tensor

	3M x 3M Matrix
	Jacobian - network jacobian representing df_i/dx_j for all nodes		

	Last Update: 03/04/22
*/

#pragma once
#include "Eigen/Eigen"

using namespace Eigen;

struct network // network attributes
{
	int num_fibers;
	int num_bnd;
	int num_nodes;
	double fiber_vol_fract;

	VectorXi fib_type;
	VectorXd fib_rads;
	VectorXd fib_areas;
	VectorXd init_lens;

	MatrixXd fibers;
	MatrixXd nodes;
	MatrixXd nodes_mapped;
	MatrixXd bnd_nodes;
};

struct fiberResults // fiber data
{
	VectorXd fib_stress;
	VectorXd fib_stretch;
	VectorXd fib_forces;
	VectorXd fib_stiffness;
};

struct netResults // network results
{
	double pressure;
	int net_flag;
	VectorXd node_forces;
	Matrix3d net_stress;
	MatrixXd component_stress;
	Matrix3d F;
	MatrixXd Jacobian;
};