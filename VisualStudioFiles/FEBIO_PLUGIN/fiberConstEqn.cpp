/* fiberConstEqn.cpp ------------------------------------------------------------------------------

SUMMARY: fiberConstEqn.cpp is a function with several fiber constitutive equations implemented.
To change the fiber constitutive equation simply change the fibertype of fibers in the network or
change the corresponging fiber const eqn in this function.
This function also calculates dF/dlambda for the constitutive equations for use in the analytical 
Jacobian calculation. 

CALLED FROM: calcForcesPeriodic.cpp, calcJacobian2.cpp
CALLS ON: None

INPUTS: vect1 - 3x1 vector along fibers
		node_force1 - 3x1 vector of forces acting on node1 of fiber
		node_force2 - 3x1 vector of forces acting on node2 of fiber
		init_len - double value of fiber initial length
		fibtype - int value of fiber type
		fib_area - double value of fiber cross-sectional area
		fib_force - double value of fiber force
		dFdlam - double value of dF/dlambda for the constitutive model

OUTPUTS: None

NOTE: node_force1, node_force2, fib_force, and dFdlam are modified by this code.

Author: Ryan Mahutga - Alford Lab - University of Minnesota

Last Update: 04-10-2023

*/

#include <iostream>
#include "Eigen/Eigen"
#include "calcForcesPeriodic.h"

using namespace Eigen;
using namespace std;

const double PI = 3.141592653589793238462643383279502884;

void fiberConstEqn(Vector3d vect1, Vector3d& node_force1, Vector3d& node_force2, double& lambda, double init_len, 
	int fibtype,double fib_area, double& fib_force, double& dFdlam)
{
	double GS = 0, // Green Strain
		fib_mod = 0, // fiber modulus
		fib_length = 0; // fiber length

	fib_length = vect1.norm(); // fiber length

	lambda = fib_length / init_len; // fiber stretch

	if (fibtype == 1 || fibtype == 12) // fiber force from helical model in Freed, A.D. & Doehring, T.C. (2005) J. Biomech. Eng. 127  
	{ // 3 is for collagen

			fib_mod = 120e6;

			double R0 = 9.8; // [nm] roughly radius of collagen microfibril
			double r0 = 1.6; // [nm] roughly a collagen triple helix
			double H0 = 67.4; //[nm] d-pattern banding in collagen molecules
		if (fibtype==12)
		{
			R0 = 8.8; // [nm] roughly radius of collagen microfibril
		}
		double L0 = sqrt((2 * PI*R0)*(2 * PI*R0) + H0 * H0);

		double lambda_bar = L0 / H0;

		double E_bar = fib_mod * H0*H0 / (L0*(H0 + (1 + 37 / (6 * PI*PI) + 2 * (L0*L0) / ((PI*r0)*(PI*r0)))*(L0 - H0)));

		double sigma = 0;

		if (lambda <= lambda_bar) // see above paper for more information
		{
			double H = lambda * H0;
			double R = sqrt(L0*L0 - H * H) / (2 * PI);
			double eta = (R*R + H * H) / (L0*H*(1 + 4 * R*R / (r0*r0) + 6 * (20 / 9 + R * R / (r0*r0))*R*R / (H*H)));

			double dHdlam = H0;
			double dRdH = -H / (2 * PI*sqrt((L0*L0 - H * H)));
			double dEtadH = -(H*H + R * R) / ((H*H) * L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +
				2 / (L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +
				12 * (R*R) * ((H*H) + (R*R))*((R*R) / (r0*r0) + 20 / 9) / ((H*H*H*H) * L0*((6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
			double dEtadR = 2 * R / (H*L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) -
				((H*H) + (R*R))*(12 * (R*R) / ((r0*r0) * (H*H)) + 12 * R*((R*R) / (r0*r0) + 20 / 9) / (H*H) + 8 * R / (r0*r0)) / (H*L0*((6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
			double dEtadlam = dEtadH * dHdlam + dEtadR * dRdH*dHdlam;

			if (lambda >= 1)
			{
				sigma = eta * E_bar*(lambda - 1);
				dFdlam = fib_area * (eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
			}
			else 
			{
				sigma = eta * E_bar*(lambda - 1);
				dFdlam = fib_area * (eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
			}
		}
		else
		{
			sigma = E_bar * (lambda_bar - 1) + fib_mod * (lambda / lambda_bar - 1);
			dFdlam = fib_mod * fib_area / lambda_bar;
		}

			fib_force = fib_area * sigma;
	}
	else if (fibtype == 2 || fibtype==0 || fibtype==5) // fiber force from F = k(GS) where GS = Green strain = 0.5(lambda^2-1) 
	{//(fibtype 2 for elastin, 0 for failed fibers)
		GS = lambda - 1.0;

		if (fibtype == 2)
		{
			if (lambda >= 1.0)
			{
				fib_mod = 1e6; //see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
				dFdlam = fib_mod * fib_area;
			}
			else
			{
				fib_mod = 0.1e6; // 1% of elastin stiffness for failed fibers
				dFdlam = fib_mod * fib_area;
			}
		}
		else if (fibtype == 5) // no tension/compression asymmetry
		{
			fib_mod = 1e6; //see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
			dFdlam = fib_mod * fib_area;
		}
		else if (fibtype == 0) // failed/removed
		{
			if (lambda >= 1.0) 
			{
				fib_mod = 1; 
				dFdlam = fib_mod * fib_area;
			}
			else
			{
				fib_mod = 1; 
				dFdlam = fib_mod * fib_area;
			}
		}
	
		fib_force = fib_mod * fib_area*GS;
	}
	else if (fibtype == 3) // actin
	{
		GS = lambda - 1.0;

		// passive contribution
		if (lambda >= 1.0) //  give 12pN force at lam = 1.05
		{
			fib_mod = 0.5e6; //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.
			dFdlam = fib_mod * fib_area;
		}
		else
		{
			fib_mod = 0.05e6; // 1% of elastin stiffness for compressed
			dFdlam = fib_mod * fib_area;
		}

		fib_force = fib_mod * fib_area * GS;

		//active contribution
		if (lambda > 0.85 && lambda < 1.85)
		{
			double S0 = 0.2e6, l_min = 0.85, l_max = 1.35, WSS_inf = 5.0; // estimate from Alford integrative biology //0.036e9; // 0.0011e9 see Elastic Proteins: biological roles and mechanical properties Gosline, J. et al. (2002) Phil. Trans. R. Soc. Lond.

			double fib_S = S0;

			double f_lam = 1.0 - (l_max - lambda)*(l_max - lambda) / ((l_max - l_min)*(l_max - l_min));

			double fib_force_active = fib_S * f_lam*fib_area;

			double dFdlamS = fib_area * 2.0 * fib_S * (l_max - lambda) / ((l_max - l_min)*(l_max - l_min));

			fib_force = fib_mod * fib_area * GS + fib_force_active;

			dFdlam = dFdlam + dFdlamS;
		}
		else 
		{
			fib_force = fib_mod * fib_area*GS;
		}
	}
	else if (fibtype == 5) // fiber force from F = EA*(exp(GS)-1) where GS = Green strain = 0.5(lambda^2-1)
	{

		double lambda_limit = 3.15;
		double B;

		fib_mod = 1e5;//2569; //  fit to Ruberti single fiber experiment
		B = 2;//77.2;

		GS = 0.5*(lambda*lambda - 1);  // Green Strain

		if (lambda >= lambda_limit)
		{
			GS = 0.5*(lambda_limit*lambda_limit - 1);  // Green Strain at limit
			double force_exp = fib_mod * fib_area*(exp(GS*B) - 1) / B;
			double slope_at_limit = fib_mod * fib_area*lambda_limit*exp(B*GS);
			fib_force = force_exp + slope_at_limit * (lambda - lambda_limit);
			dFdlam = slope_at_limit;
		}
		else if (lambda >= 1)
		{
			fib_force = fib_mod * fib_area*(exp(B*GS) - 1) / B;
			dFdlam = fib_mod * fib_area*lambda*exp(B*GS);
		}
		else if (lambda < 1)
		{
			fib_force = fib_mod * fib_area*(exp(B*GS) - 1) / B;
			dFdlam = fib_mod * fib_area*lambda*exp(B*GS);
		}
	}

	// angles for forces
	Vector3d cosine = vect1 / fib_length;

	// nodal forces due to stretched fibers
	node_force1 = fib_force*cosine;
	node_force2 = -1.0 * fib_force*cosine;

}