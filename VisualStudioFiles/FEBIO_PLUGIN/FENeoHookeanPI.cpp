
#include "FENetwork.h"
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

//-----------------------------------------------------------------------------
// These macros define the material parameter list for the FENetwork material.
//
// The BEGIN_FECORE_CLASS macro takes the material class and its base class as 
// parameters. 
//
// The ADD_PARAMETER macro defines the actual parameter. It takes three parameters:
// - the variable as defined in the class.
// - a range speficier which defines the valid range of the parameter
// - a string that defines the name of the parameter as it will appear in the input file.
//
// The END_PARAMETER_LIST macro just defines the end of the parameter list.
BEGIN_FECORE_CLASS(FENetwork, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENetwork::FENetwork(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
// This (optional) function is called during initialization. This can be done
// to do one-time initialization. Since for this material this is not required,
// this function could have been ommitted. 
// Make sure to always call the base class (usually first). 
bool FENetwork::Init()
{
	// Don't forget the base class initialization first.
	if (FEElasticMaterial::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
// This (optional) function is used to validate material parameters.
// It is recommended to provide a valid range during material parameter definition
// (using ADD_PARAMETER2), but any additional verification can be done here.
// For instance, here it is done to set the m_lam and m_mu parameters which depend on the 
// material parameters. 
// Make sure to always call the base class (usually first).
bool FENetwork::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;

	// The Stress and Tangent functions are written in terms of the Lame parameters
	// so we calculate these here.
	m_lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	m_mu  = 0.5*m_E/(1+m_v);

	return true;
}

//-----------------------------------------------------------------------------
// This function needs to return the spatial (i.e. Cauchy stress) at the material point
// which is passed as a parameter. The FEMaterialPoint class contains all the state 
// variables of the material (and its base classes).
mat3ds FENetwork::Stress(FEMaterialPoint& mp)
{
	// The FEMaterialPoint classes are stored in a linked list. The specific material
	// point data needed by this function can be accessed using the ExtractData member.
	// In this case, we want to FEElasticMaterialPoint data since it stores the deformation
	// information that is needed to evaluate the stress.
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// We'll need the deformation gradient and its determinant in this function.
	// Note that we don't take the determinant of F directly (using mat3d::det)
	// but instead use the m_J member variable of FEElasticMaterialPoint.
	mat3d &Fb = pt.m_F;

	Matrix3d F;

	F(0, 0) = Fb(0, 0);
	F(0, 1) = Fb(0, 1);
	F(0, 2) = Fb(0, 2);

	F(1, 0) = Fb(1,0);
	F(1, 1) = Fb(1,1);
	F(1, 2) = Fb(1,2);

	F(2, 0) = Fb(2,0);
	F(2, 1) = Fb(2,1);
	F(2, 2) = Fb(2,2);

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
	// --------

	fiberResults fiberResults_n; // fiber esults
	netResults netResults_n; // network results

	// solve network
	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, 1, stab_node);

	Matrix3d net_stress = netResults_n.net_stress; // pull network stress

	double c_star = 150, // M - bathing solution concentration
		T = 310, // temperature K
		k = 200, //k gives cfcd give 60mEq/L at 30% actin
		J = F.determinant(), // volume of network
		phim = fib_vol_m / (x_scale * x_scale * J); // fiber fraction of actin

	const double R = 8.314; // universal gas constant
	double c_fcd = k * phim; // fixed charge density of actin

	double pressure = R * T * (sqrt(c_fcd * c_fcd + 4 * c_star * c_star) - 2 * c_star);// hydrostatic pressure
	netResults_n.pressure = pressure; // store pressure

	// total stress
	net_stress(0, 0) = net_stress(0, 0) - pressure;
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	//cout << "Network Solve. " << endl << F << endl << net_stress << endl << endl;
	// This is the actual computation of the Cauchy stress. 
	double detF = pt.m_J;

	// The FEElasticMaterialPoint class defines several useful functions for 
	// evaluating strain measures, such as the left Cauchy-Green tensor.
	mat3ds b = pt.LeftCauchyGreen();

	// This creates the second-order identity tensor which we need
	// to evaluate the Cauchy stress.
	mat3dd I(1);

	// This is the actual computation of the Cauchy stress. 
	mat3ds s;// = (b - I) * (m_mu / detF) + I * (m_lam * log(detF) / detF);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			s(i, j) = net_stress(i, j);
		}
	}

	// The Cauchy stress is returned.
	return s;
}

//-----------------------------------------------------------------------------
// This function calculates the spatial elasticity tangent tensor. 
// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
// which is a fourth-order tensor with major and minor symmetries.
tens4ds FENetwork::Tangent(FEMaterialPoint& mp)
{
	// As in the Stress function, we need the data from the FEElasticMaterialPoint
	// class to calculate the tangent.
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// Get the deformation gradient and its determinant
	mat3d &Fb = pt.m_F;
	double detF = pt.m_J;

	Matrix3d F;

	F(0, 0) = Fb(0, 0);
	F(0, 1) = Fb(0, 1);
	F(0, 2) = Fb(0, 2);

	F(1, 0) = Fb(1, 0);
	F(1, 1) = Fb(1, 1);
	F(1, 2) = Fb(1, 2);

	F(2, 0) = Fb(2, 0);
	F(2, 1) = Fb(2, 1);
	F(2, 2) = Fb(2, 2);

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
	// --------

	fiberResults fiberResults_n; // fiber esults
	netResults netResults_n; // network results

	// solve network
	networkSolver(num_nodes, num_fibers, x_scale, F, network0,
		network_n, fiberResults_n, netResults_n, 1, stab_node);

	Matrix3d net_stress = netResults_n.net_stress; // pull network stress

	double c_star = 150, // M - bathing solution concentration
		T = 310, // temperature K
		k = 200, //k gives cfcd give 60mEq/L at 30% actin
		J = F.determinant(), // volume of network
		phim = fib_vol_m / (x_scale * x_scale * J); // fiber fraction of actin

	const double R = 8.314; // universal gas constant
	double c_fcd = k * phim; // fixed charge density of actin

	double pressure = R * T * (sqrt(c_fcd * c_fcd + 4 * c_star * c_star) - 2 * c_star);// hydrostatic pressure
	netResults_n.pressure = pressure; // store pressure

	// total stress
	net_stress(0, 0) = net_stress(0, 0) - pressure;
	net_stress(1, 1) = net_stress(1, 1) - pressure;
	net_stress(2, 2) = net_stress(2, 2) - pressure;

	double c_sp[3][3][3][3];
	elasticity(network0, num_fibers, F, net_stress, x_scale, fib_vol_m, c_sp);

	// Calculate the modified Lame parameters
	double lam1 = m_lam / detF;
	double mu1  = (m_mu - m_lam*log(detF)) / detF;

	// define identity tensor and some useful
	// dyadic products of the identity tensor.
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	// evaluate the elasticity tensor
	tens4ds c;//= I4*(2.0*mu1) + IxI*lam1;

	for (int I = 0; I < 3; I++)
	{
		for (int J = 0; J < 3; J++)
		{
			for (int K = 0; K < 3; K++)
			{
				for (int L = 0; L < 3; L++)
				{
					c(I, J, K, L) = c_sp[I][J][K][L];
				}
			}
		}
	}

	// return the elasticity tensor
	return c;
}
