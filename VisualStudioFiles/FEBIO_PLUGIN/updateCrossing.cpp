/* updateCrossing.cpp updates the fiber vector according to how nodes enter and exit the RVE

SUMMARY: This function updates the crossing conditions in 'fibers' based on 'nodes_n'. It finds
the edges of the RVE with a margin of error, and determines if a node on the end of a fiber is 
outside the boundary. If it is, the fiber crossing is updated and the nodes is reflected back
into the RVE.

CALLED FROM: calcForcesPeriodic.cpp & calcBndNodes.cpp

CALLS ON: None

INPUTS: nodes_n - Nx3 matrix of nodal coordinates
		fibers - Mx5 matrix of fiber connectivity
		num_fibers - M number of fibers
		rve_stretch - 3x1 vector of primary RVE stretches

OUTPUTS: None

NOTE: This function does modify the values for nodes_n and fibers in main.

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota
Contact: mahut005@umn.edu

Last Update: 04-02-2018
*/

#include <iostream>
#include "Eigen/Eigen"
#include <complex>
#include <cmath>
#include "calcForcesPeriodic.h"

const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

void updateCrossing(MatrixXd& nodes_mapped, MatrixXd& nodes_n, MatrixXd& fibers_n, const int num_nodes, 
const int num_fibers, Matrix3d F)
{
	// defining local values
	double marg = 1e-16; // margin of error for finding bounds
	int node1 = 0, // static node in RVE on fiber from node1 to node2
		node2 = 0, // dynamic fiber node (i.e. node that moves to be the nearest neighbor to node1 along the fiber
		q = 0, p = 0,k = 0, h = 0; // crossing counters

	double current_node1, // current nodal coordinates of static node
		   current_node2; // current nodal coordinates of dynamic node

	// defining local data structures
	MatrixXd nodes2 = nodes_mapped; // unchanging node positions

	// finding boundaries (always 0.5 and -0.5 for mapped networks)
	double up_bnd = 0.5, // upper boundaries defined with margin of error
		   lo_bnd = -0.5; // lower bound in all cases

	// looping fibers to update fiber vector crossings, shift nodes, and find fiber lengths and direction vectors
	for (int n = 0; n < num_fibers; n++)
	{
		node1 = (int)(fibers_n(n, 0) + 0.0003); // first (static) node of fiber 
		node2 = (int)(fibers_n(n, 1) + 0.0003); // second (dynamic) node of fiber (i.e. node that is moved according to fibers(m,3:5))

		for (int i = 0; i < 3; i++) // 0=X, 1=Y, 2=Z
		{
			current_node1 = nodes2(node1, i); // static node

			current_node2 = nodes2(node2, i) ; // dynamic node

			q = 0; p = 0; h = 0; k = 0; // setting counters to zero

			while (current_node1 > up_bnd || current_node2 > up_bnd || // out of upper bound
				current_node1 < lo_bnd || current_node2 < lo_bnd) // out of lower bound
			{
				if (current_node1 > up_bnd) // seeing if node1 is in the RVE
				{
					fibers_n(n, 2 + i) = (int)(fibers_n(n,2 + i) - 1.2);
					current_node1 = current_node1 - 1.0;
					q++;
					nodes_mapped(node1, i) = nodes2(node1, i) - q*(1.0);
				}
				else if (current_node1 < lo_bnd)
				{
					fibers_n(n, 2 + i) = (int)(fibers_n(n, 2 + i) + 1.2);
					p++;
					current_node1 = current_node1 + 1.0;
					nodes_mapped(node1, i) = nodes2(node1, i) + p*(1.0);
				}
				if (current_node2 > up_bnd) // seeing if node2 is in the RVE
				{
					fibers_n(n, 2 + i) = (int)(fibers_n(n, 2 + i) + 1.2);
					h++;
					current_node2 = current_node2 - 1.0;
					nodes_mapped(node2, i) = nodes2(node2, i) - h*(1.0);
				}
				else if (current_node2 < lo_bnd)
				{
					fibers_n(n, 2 + i) = (int)(fibers_n(n, 2 + i) - 1.2);
					k++;
					current_node2 = current_node2 + 1.0;
					nodes_mapped(node2, i) = nodes2(node2, i) + k*(1.0);
				}
			}
		}
	}

	// update deformed nodal positions
	for (int n = 0; n < num_nodes; n++)
	{
		nodes_n(n, 0) = F(0, 0)*(nodes_mapped(n, 0)) + F(0, 1)*(nodes_mapped(n, 1)) + F(0, 2)*(nodes_mapped(n, 2));
		nodes_n(n, 1) = F(1, 0)*(nodes_mapped(n, 0)) + F(1, 1)*(nodes_mapped(n, 1)) + F(1, 2)*(nodes_mapped(n, 2));
		nodes_n(n, 2) = F(2, 0)*(nodes_mapped(n, 0)) + F(2, 1)*(nodes_mapped(n, 1)) + F(2, 2)*(nodes_mapped(n, 2));
	}

}


