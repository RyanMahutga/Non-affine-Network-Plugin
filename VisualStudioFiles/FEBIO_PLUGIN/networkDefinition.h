
#include "Eigen/Eigen"

using namespace Eigen;

struct network
{
	int num_fibers;
	int num_bnd;
	int num_nodes;

	VectorXd fib_type;
	VectorXd fib_rads;
	VectorXd fib_areas;
	VectorXd init_lens;

	MatrixXd fibers;
	MatrixXd nodes;
	MatrixXd bnd_nodes;
};