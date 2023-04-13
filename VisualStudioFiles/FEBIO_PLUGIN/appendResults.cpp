/* appendResults.cpp -------------------------------------------------------------------------------

SUMMARY: This function appends results to various text files for each iteration in main

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

void appendResults(Matrix3d F, MatrixXd component_stress, Matrix3d net_stress, VectorXd fib_stretch,
	double pressure, int time)
{
	string str;

	// saving stretch results
	ofstream myfile1;
	str = "F" + to_string(time) + ".txt";
	myfile1.open(str, ofstream::app);
	myfile1 << F << endl;
	myfile1.close();

	// saving stress results
	ofstream myfile2;
	str = "stress" + to_string(time) + ".txt";
	myfile2.open(str, ofstream::app);
	myfile2 << net_stress << endl;
	myfile2.close();

	// saving stress results
	ofstream myfile3;
	str = "comp_stress" + to_string(time) + ".txt";
	myfile3.open(str, ofstream::app);
	myfile3 << component_stress << endl;
	myfile3.close();

	ofstream myfile4;
	str = "fib_stretch" + to_string(time) + ".txt";
	myfile4.open(str, ofstream::app);
	myfile4 << fib_stretch << endl;
	myfile4.close();

	ofstream myfile5;
	str = "pressure" + to_string(time) + ".txt";
	myfile5.open(str, ofstream::app);
	myfile5 << pressure << endl;
	myfile5.close();
}