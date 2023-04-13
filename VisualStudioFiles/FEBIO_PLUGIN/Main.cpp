/* Main.cpp -------------------------------------------------------------------------------------------
NOTE: To compile this code you will need to have the Eigen library. (https://eigen.tuxfamily.org/)

DETAILED DESCRIPTION: 

INPUTS: 

OUTPUTS:


Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Contact: mahut005@umn.edu

Last Modifed: 03/04/2022
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>
#include <cstdio>
#include <ctime>
#include "dataStructureDefinitions.h"

const double R = 8.314; // universal gas constant
const double PI = 3.141592653589793238462643383279502884;

using namespace std;
using namespace Eigen;

#define EIGEN_USE_MKL_ALL

// Main function
int main()
{
	// reading in data from PeriodicNetwork.txt
	ifstream inFile, inFile1, inFile2, inFile3;

	// Load stress free network ---------------------------------------------------------------
	network network0; // underformed network structure definition
	double x_scale; // RVE scale
	VectorXi stab_node(3); stab_node = VectorXi::Zero(3); // stabalization node

	inFile.open("PeriodicNetwork_new0.txt");
	if (inFile.fail()) { // if can't open file
		//cout << "Unable to Open Network File!" << endl;
		//cin.get();
		exit(1); // terminate with error
	}

	int num_fibers = 0;
	int num_nodes = 0;
	// reading in first line (number of fibers, number of nodes)
	inFile >> num_fibers >> num_nodes >> x_scale>> stab_node(0) >> stab_node(1) >> stab_node(2);

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

	// ----------------------------------------------------------------------------

	// Loading in F0 ------------------------------------------------------------
	Matrix3d F(3, 3); F = Matrix3d::Zero(3, 3); // network deformation gradient
	int flag, flag0;

	inFile1.open("F0.txt");
	if (inFile1.fail()) { // if can't open file
		//cout << "No F0 File Available." << endl;
		F << 1.0, 0,0,0, 1.0, 0, 0 ,0 , 1.0;
		flag0 = 0;
	}
	else
	{
		flag0 = 1;

		while (!inFile1.eof()) // read file until end
		{
			inFile1 >> F(0, 0) >> F(0, 1) >> F(0, 2) >> F(1, 0) >> F(1, 1) 
				>> F(1, 2) >> F(2, 0) >> F(2, 1) >> F(2, 2);
			
		}
	}
	inFile1.close();

	double VV = F.determinant();
	MatrixXd F_inv(3, 3); F_inv = F.inverse(); // inverse deformation mapping
	// ------------------------------------------------------------------------------------

	// Loading in finite element deformation ------------------------------------------
	Matrix3d F_max;
	inFile3.open("FFE.txt");
	if (inFile3.fail()) { // if can't open file
		//cout << "Single Fit Deformation File does not exist." << endl << endl;
		flag = 0;
	}
	else
	{
		flag = 1;
		while (!inFile3.eof()) // read file until end
		{
			inFile3 >> F_max(0, 0) >> F_max(0, 1) >> F_max(0, 2) >> F_max(1, 0) >> F_max(1, 1)
				>> F_max(1, 2) >> F_max(2, 0) >> F_max(2, 1) >> F_max(2, 2);
		}
	}
	// -----------------------------------------------------------------------------------------------------

	// Load in testing limits ------------------------------------------------------------
	double lam_eq = 1.1, lam_uni=1.3 , lam_shear=0.1, lam_tri=1.05;

	if (flag == 0) // if FFE doesn't exist
	{
		inFile2.open("lam_lims.txt");
		if (inFile2.fail()) { // if can't open file
			//cout << "Unable to Open Limits File!" << endl;
			//cin.get();
			exit(1); // terminate with error
		}

		inFile2 >> lam_eq >> lam_uni >> lam_shear >> lam_tri;

		inFile2.close();
	}
	//----------------------------------------------------------------------------------

// fiber volume calculations --------------------------------------------------------
	double  phi_e = 0, phi_m = 0, phi_c = 0, // fiber fractions
		fib_vol_m = 0, fib_vol_e = 0, fib_vol_c = 0; // fiber volumes

	for (int f = 0; f < num_fibers; f++)
	{
		if (fib_type(f) == 1)
		{
			fib_vol_m = fib_vol_m + fib_areas(f) * init_lens(f); // fib vol actin
		}
		else if (fib_type(f) == 2)
		{
			fib_vol_e = fib_vol_e + fib_areas(f) * init_lens(f); // fib vol elastin
		}
		else
		{
			fib_vol_c = fib_vol_c + fib_areas(f) * init_lens(f); // fib vol collagen
		}
	}

	phi_m = fib_vol_m / (x_scale * x_scale * VV);
	phi_e = fib_vol_e / (x_scale * x_scale * VV);
	phi_c = fib_vol_c / (x_scale * x_scale * VV);

	double fib_vol = fib_vol_m + fib_vol_c + fib_vol_e; // total fiber volume
// ----------------------------------------------------------------------------
	network old_net = network0;

	Matrix3d F_c, F0, F_prev, net_stress, net_stress0; // deformations, net stress

	Vector3d x_guess, init_err; 
	init_err << 0, 0, 0; 
	x_guess << F(0, 0), F(1, 1), F(2, 2); // solution guess

	VectorXd zeroStress(6), error(6); zeroStress << 0, 0, 0, 0, 0, 0; error = zeroStress;  // zero stress BC

	double pressure = 0, pressure0 = 0; // flow rate and pressure terms

	fiberResults fiberResults_n; // fiber esults
	netResults netResults_n; // network results

	VectorXd fib_stretch(num_fibers); fib_stretch = VectorXd::Zero(num_fibers); // fiber stretches

	MatrixXd comp(12, 3); comp = MatrixXd::Zero(12, 3); // component stresses
	VectorXd Om(18), Om0(18); Om = VectorXd::Zero(18); // orientation

	// calculate zero stress
	if (flag0==0)
	{
		jacobianSolver3_0(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n, 0,
		stab_node, zeroStress, net_stress, pressure, x_guess, fib_vol, phi_m, 0 , 0);
	}
	else
	{
		solveJac3_0(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
		0, stab_node, zeroStress, net_stress, init_err, 0, fib_vol);

		netResults_n.F = F;
	}
	
	F = netResults_n.F; // deformation zero stress
	net_stress = netResults_n.net_stress; // stress at zero stress
	pressure = netResults_n.pressure; // pressure at zero stress
	net_stress0 = net_stress; // stress0
	pressure0 = pressure; // pressure0

	// fiber faction
	double V = F.determinant(), J = V, fiber_vol_fract = fib_areas.dot(init_lens) / (x_scale * x_scale * V);
	network0.fiber_vol_fract = fiber_vol_fract;
	
	// reset data for each step
	MatrixXd nodes00 = network0.nodes;
	MatrixXd fibers00 = network0.fibers;
	MatrixXd nodes0n = network_n.nodes;
	MatrixXd fibers0n = network_n.fibers;

	// setting deformations
	F0 = netResults_n.F;
	F_c = netResults_n.F;
	F_prev = netResults_n.F;

	// simulation parameters
	double n_pts = 19.0, lam, dlam, ratio;
	int test; // test number

	if (flag == 0) // if FFE doesn't exist
	{
		for (int te = 0; te < 10; te++)
		{
			test = te;
			if (test == 0 || test == 1 || test == 2) //equibiaxial
			{
				lam = lam_eq;
				dlam = (lam - 1.0) / n_pts;
			}
			else if (test == 3 || test == 4 || test == 5) // uniaxial
			{
				lam = lam_uni;
				dlam = (lam - 1.0) / n_pts;
			}
			else if (test == 6 || test == 7 || test == 8) // shear
			{
				lam = lam_shear;
				dlam = lam / n_pts;
			}
			else
			{
				lam = lam_tri;
				dlam = (lam - 1.0) / n_pts;
			}

			// reset F_c and F for each test
			F_c = F0;
			F = F0;

			// reset fibers and nodes for each test
			network0.fibers = fibers00;
			network_n.fibers = fibers0n;
			network0.nodes = nodes00;
			network_n.nodes = nodes0n;

			for (int p = 0; p < n_pts; p++)
			{
				// overwrite first iteration
				if (p == 0)
				{
					overwriteResults(F0, comp, net_stress0, fib_stretch, pressure0, test+1);
				}

				F_prev = F_c; // set F_prev to F_c

				if (test == 0) // xy biaxial
				{
					F_c(0, 0) += dlam * F0(0, 0);
					F_c(1, 1) += dlam * F0(1, 1);
					F_c(2, 2) = V / (F_c(0, 0) * F_c(1, 1));
				}
				else if (test == 1) //xz biaxial
				{
					F_c(0, 0) += dlam * F0(0, 0);
					F_c(2, 2) += dlam * F0(2, 2);
					F_c(1, 1) = V / (F_c(0, 0) * F_c(2, 2));
				}
				else if (test == 2) // yz biaxial
				{
					F_c(1, 1) += dlam * F0(1, 1);;
					F_c(2, 2) += dlam * F0(2, 2);;
					F_c(0, 0) = V / (F_c(1, 1) * F_c(2, 2));
				}
				else if (test == 3) // x uniaxial
				{
					F_c(0, 0) += dlam * F0(0, 0);
					ratio = F0(2, 2) / F0(1, 1);
					F_c(1, 1) = sqrt(V / (ratio*F_c(0, 0)));
					F_c(2, 2) = (V / (F_c(0, 0)*F_c(1, 1)));
				}
				else if (test == 4) // y uniaxial
				{
					F_c(1, 1) += dlam * F0(1, 1);
					ratio = F0(2, 2) / F0(0, 0);
					F_c(0, 0) = sqrt(V / (ratio*F_c(1, 1)));
					F_c(2, 2) = (V / (F_c(0, 0) * F_c(1, 1)));
				}
				else if (test == 5) // z uniaxial
				{
					F_c(2, 2) += dlam * F0(2, 2);
					ratio = F0(1,1) / F0(0, 0);
					F_c(0, 0) = sqrt(V / (ratio*F_c(2, 2)));
					F_c(1,1) = (V / (F_c(0, 0) * F_c(2, 2)));
				}
				else if (test == 6)// xy shear
				{
					F_c(0, 1) += dlam * F0(1, 1);
					F_c(1, 0) += dlam * F0(0, 0);
					F_c(2, 2) = V / (F_c(0, 0) * F_c(1, 1) - F_c(0, 1) * F_c(1, 0));
				}
				else if (test == 7)// xz shear
				{
					F_c(0, 2) += dlam * F0(2, 2);
					F_c(2, 0) += dlam * F0(0, 0);
					F_c(1, 1) = V / (F_c(0, 0) * F_c(2, 2) - F_c(0, 2) * F_c(2, 0));
				}
				else if (test == 8)//yz shear
				{
					F_c(1, 2) += dlam * F0(2, 2);
					F_c(2, 1) += dlam * F0(1, 1);
					F_c(0, 0) = V / (F_c(2, 2) * F_c(1, 1) - F_c(2, 1) * F_c(1, 2));
				}
				else if (test == 10)// strip
				{
					F_c(0, 0) += dlam * F0(0, 0);
					F_c(1, 1) = F0(1, 1);
					F_c(2, 2) = F0(2, 2);
				}
				else if (test == 11)// strip
				{
					F_c(0, 0) = F0(0, 0);
					F_c(1, 1) += dlam * F0(1, 1);
					F_c(2, 2) = F0(2, 2);
				}
				else if (test == 12)//  strip
				{
					F_c(0, 0) = F0(0, 0);
					F_c(1, 1) = F0(1, 1);
					F_c(2, 2) += dlam * F0(2, 2);
				}
				else if (test == 9)// equi-triaxial
				{
					F_c(0, 0) += dlam * F0(0, 0);
					F_c(1, 1) += dlam * F0(1, 1);
					F_c(2, 2) += dlam * F0(2, 2);
				}

				F = F_c; // setting F

				// run simulation for F_c
				solveJac(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
					1, stab_node, zeroStress, net_stress, error, 0, fib_vol);

				// pressure
				pressure = netResults_n.pressure;

				// stresses
				comp = netResults_n.component_stress;
				net_stress = netResults_n.net_stress;
				// fiber stretches
				fib_stretch = fiberResults_n.fib_stretch;
				// appending results
				appendResults(F, comp, net_stress, fib_stretch, pressure, test+1);
			}
		}
	}
	else if (flag == 1) // for FFE deformation
	{
		for (int te = 0; te < 2; te++)
			{
			test = 96 + te;

			// reset fibers and nodes for each test
			network0.fibers = fibers00;
			network_n.fibers = fibers0n;
			network0.nodes = nodes00;
			network_n.nodes = nodes0n;

			// reset F_c and F for each test
			F_c = F0;
			F = F0;

			for (int p = 0; p < n_pts; p++)
			{
				// overwrite first iteration
				if (p == 0)
				{
					overwriteResults(F0, comp, net_stress0, fib_stretch, pressure0, test);
				}

				if (te == 0) // test going from F0 to FFE
				{
					dlam = (F_max(0, 0) - F0(0, 0)) / n_pts;
					F_c(0, 0) += dlam;

					dlam = (F_max(0, 1) - F0(0, 1)) / n_pts;
					F_c(0, 1) += dlam;

					dlam = (F_max(0, 2) - F0(0, 2)) / n_pts;
					F_c(0, 2) += dlam;

					dlam = (F_max(1, 0) - F0(1, 0)) / n_pts;
					F_c(1, 0) += dlam;

					dlam = (F_max(1, 1) - F0(1, 1)) / n_pts;
					F_c(1, 1) += dlam;

					dlam = (F_max(1, 2) - F0(1, 2)) / n_pts;
					F_c(1, 2) += dlam;

					dlam = (F_max(2, 0) - F0(2, 0)) / n_pts;
					F_c(2, 0) += dlam;

					dlam = (F_max(2, 1) - F0(2, 1)) / n_pts;
					F_c(2, 1) += dlam;

					dlam = (F_max(2, 2) - F0(2, 2)) / n_pts;
					F_c(2, 2) += dlam;
				}
				else if (te == 1) // uniacxial stretch x
				{
					dlam = (F_max(0, 0) - F0(0, 0)) / n_pts;
					F_c(0, 0) += dlam;
				}
				else if (te == 2) // uniaxial stretch y
				{
					dlam = (F_max(1, 1) - F0(1, 1)) / n_pts;
					F_c(1, 1) += dlam;
				}
				else if (te == 3) // uniaxial stretch z
				{
					dlam = (F_max(2, 2) - F0(2, 2)) / n_pts;
					F_c(2, 2) += dlam;
				}

				F = F_c; // set F

				// solve network
				solveJac(num_nodes, num_fibers, x_scale, F, network0, network_n, fiberResults_n, netResults_n,
					1, stab_node, zeroStress, net_stress, error, 0, fib_vol);

				// stress
				comp = netResults_n.component_stress;
				net_stress = netResults_n.net_stress;
				// fiber stretch
				fib_stretch = fiberResults_n.fib_stretch;
				// appending results
				appendResults(F, comp, net_stress, fib_stretch, pressure, test);
			}
		}
	}
	return 0;
}

