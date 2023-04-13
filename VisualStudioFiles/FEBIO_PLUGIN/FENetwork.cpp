
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
	ADD_PARAMETER(netNumber, FE_RANGE_DONT_CARE(), "netNum"); // network number identifier N in PeriodicNetworkN.txt
	ADD_PARAMETER(netSave, FE_RANGE_DONT_CARE(), "netSave"); // if 1 save the network data
	ADD_PARAMETER(netSolve, FE_RANGE_DONT_CARE(), "netSolve"); // if 1 solve the network problem
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E"); // Young's Modulus
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v"); // Poisson's Ratio
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

	Matrix3d F;F << Fb(0, 0), Fb(0, 1), Fb(0, 2), 
					Fb(1,0), Fb(1,1), Fb(1,2),
					Fb(2,0), Fb(2,1), Fb(2,2);

	Matrix3d net_stress; net_stress << 0, 0, 0, 0, 0, 0, 0, 0, 0; // initialize net stress

	// initialize network and fiber results
	network network0; fiberResults fiberResults_n;
	
	// solve stress given F
	stressCalc(F, netNumber, netSave, netSolve, network0, fiberResults_n, net_stress);

	double detF = pt.m_J;
	mat3dd I(1);

	// The FEElasticMaterialPoint class defines several useful functions for 
// evaluating strain measures, such as the left Cauchy-Green tensor.
	mat3ds b = pt.LeftCauchyGreen();

	// This is the actual computation of the Cauchy stress. 
	mat3ds s, s_bulk = (b - I) * (m_mu / detF) + I * (m_lam * log(detF) / detF);

	// converts from eigen matrix to mat3ds in FEBio
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			s(i, j) = net_stress(i, j);
		}
	}

	s = s + s_bulk; // add network and neo-Hookean stress

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

	// deformation matrix in eigen matrix form
	Matrix3d F; F << Fb(0, 0), Fb(0, 1), Fb(0, 2),
		Fb(1, 0), Fb(1, 1), Fb(1, 2),
		Fb(2, 0), Fb(2, 1), Fb(2, 2);

	mat3d U_i = pt.RightStretch(); // right stretch tensor
	Matrix3d U; U << U_i(0, 0), U_i(0, 1), U_i(0, 2),
		U_i(1, 0), U_i(1, 1), U_i(1, 2),
		U_i(2, 0), U_i(2, 1), U_i(2, 2);

	Matrix3d U_inv = U.inverse(); // inverse of right stretch

	Matrix3d R = U_inv * F; // rotation matrix

	Matrix3d net_stress; net_stress << 0, 0, 0, 0, 0, 0, 0, 0, 0; // initialize net stress

	// initialize netwock and fiber results
	network network0; fiberResults fiberResults_n;

	// calculate stress
	stressCalc(F, netNumber, 0, netSolve, network0, fiberResults_n, net_stress);

	double fib_vol_m = network0.fib_areas.transpose() * network0.init_lens; 

	double c_sp[3][3][3][3]; // initialize elasticity array

	// calculate elasticity
	networkElasticity(F, U_inv, U, R, network0,fiberResults_n, net_stress, fib_vol_m, c_sp);

	// evaluate the elasticity tensor
		// Calculate the modified Lame parameters
	// Calculate the modified Lame parameters
	double lam1 = m_lam / detF;
	double mu1 = (m_mu - m_lam * log(detF)) / detF;

	// define identity tensor and some useful
	// dyadic products of the identity tensor.
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	// evaluate the elasticity tensor
	tens4ds c, c_bulk = I4 * (2.0 * mu1) + IxI * lam1;

	// translate array to tens4d in FEBio
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

	c = c + c_bulk; // sum network and neo-Hookean elasticity

	// return the elasticity tensor
	return c;
}
