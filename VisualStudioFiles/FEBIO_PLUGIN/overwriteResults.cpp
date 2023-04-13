/* overwriteResults.cpp -------------------------------------------------------------------------------

SUMMARY: This function overwrites the stress-stretch data used to fit material properties.

CALLED FROM: Main

CALLS ON: None

INPUTS:

OUTPUTS: None

NOTE:

Author: Ryan Mahutga - Barocas Research Group - University of Minnesota

Last Update: 03-04-2022
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "math.h"
#include "Eigen/Eigen"
#include "Eigen/Sparse"
#include "calcForcesPeriodic.h"
#include <vector>
#include <cstdio>
#include <ctime>
#include "dataStructureDefinitions.h"

using namespace std;
using namespace Eigen;

void overwriteResults(Matrix3d F, MatrixXd component_stress, Matrix3d net_stress, VectorXd fib_stretch,
	VectorXd fib_stress, int time)
{
	string str;

	// saving stretch results
	ofstream myfile1;
	str = "F" + to_string(time) + ".txt";
	myfile1.open(str);
	myfile1 << F << endl;
	myfile1.close();

	// saving stress results
	ofstream myfile2;
	str = "stress" + to_string(time) + ".txt";
	myfile2.open(str);
	myfile2 << net_stress << endl;
	myfile2.close();

	// saving stress results
	ofstream myfile3;
	str = "comp_stress" + to_string(time) + ".txt";
	myfile3.open(str);
	myfile3 << component_stress << endl;
	myfile3.close();

	// saving fiber stretches
	ofstream myfile4;
	str = "fib_stretch" + to_string(time) + ".txt";
	myfile4.open(str);
	myfile4 << fib_stretch << endl;
	myfile4.close();

	// saving fiber stresses
	ofstream myfile5;
	str = "fib_stress" + to_string(time) + ".txt";
	myfile5.open(str);
	myfile5 << fib_stress << endl;
	myfile5.close();
}