
#define _CRT_SECURE_NO_WARNINGS

#include "FENetwork.h"
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

using namespace std;
using namespace Eigen;

// stress calculation for network
void stressCalc(Matrix3d F, int netNum, double netSave, double netSolve, network& network0, fiberResults& fiberResults_n, Matrix3d& net_stress)
{
	// reading in data from PeriodicNetwork.txt
	ifstream inFile, inFile1, inFile2, inFile3;

	// Load stress free network ---------------------------------------------------------------
	double x_scale; // RVE scale
	VectorXi stab_node(3); stab_node = VectorXi::Zero(3); // stabilization node

	string s = "PeriodicNetwork" + to_string(netNum) + ".txt";

	inFile.open(s);
	if (inFile.fail()) { // if can't open file
		cout << "Unable to Open Network File!" << endl;
		//cin.get();
		exit(1); // terminate with error
	}

	int num_fibers = 0;
	int num_nodes = 0;
	// reading in first line (number of fibers, number of nodes)
	inFile >> num_fibers >> num_nodes >> x_scale >> stab_node(0) >> stab_node(1) >> stab_node(2);

	network0.fiber_vol_fract = x_scale;
	network0.num_fibers = num_fibers;
	network0.num_nodes = num_nodes;

	// initialize values
	MatrixXd fibers = MatrixXd::Zero(num_fibers, 5);
	MatrixXd nodes = MatrixXd::Zero(num_nodes, 3);
	VectorXd fib_rads = VectorXd(num_fibers);
	VectorXd init_lens = VectorXd(num_fibers);
	VectorXi fib_type = VectorXi::Zero(num_fibers);

	int count = -1; // counter for fibers
	int count2 = 0; // counter for nodes
	while (!inFile.eof()) // read file until end
	{
		if (count < num_fibers && count > -1)
		{
			inFile >> fibers(count, 0) >> fibers(count, 1) >> fibers(count, 2) >> fibers(count, 3)
				>> fibers(count, 4) >> fib_type(count) >> init_lens(count) >> fib_rads(count);
		}
		else if (count >= num_fibers && count < num_fibers + num_nodes)
		{
			inFile >> nodes(count2, 0) >> nodes(count2, 1) >> nodes(count2, 2);
			count2++;
		}
		else if (count > -1)
		{
			//cout << "Problem Reading in Data!" << endl << endl;
			//cout << "Check to ensure there are no blank lines in input file." << endl;
			//cin.get();
			exit(1);
		}
		count++;
	}

	inFile.close();

	network0.nodes = nodes;
	network0.fibers = fibers;
	network0.fib_type = fib_type;
	network0.init_lens = init_lens;
	network0.fib_rads = fib_rads;
	VectorXd fib_areas = PI * fib_rads.cwiseProduct(fib_rads); // fiber areas
	network0.fib_areas = fib_areas;

	MatrixXd nodes_n = nodes;
	MatrixXd fibers_n = fibers;

	network network_n = network0; // deformed/modified network structure definition

	double VV = F.determinant();
	MatrixXd F_inv(3, 3); F_inv = F.inverse(); // inverse deformation mapping

	netResults netResults_n; // network results
	// solve network
	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, 1, stab_node, netSolve);

	MatrixXd component_stress = netResults_n.component_stress; // component stresses
	net_stress = netResults_n.net_stress; // pull network stress
	VectorXd fib_stretch = fiberResults_n.fib_stretch; // fiber stretch
	VectorXd fib_stress = fiberResults_n.fib_stress; // fiber stress

	if (netSave>0.99)
	{
		overwriteResults(F, component_stress, net_stress, fib_stretch, fib_stress, netNum);
	}
}