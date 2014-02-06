/* Definitions for types used.  In particular, typenames for matrix
 * objects in the eigen package v3 for matrices and linear algebra are
 * defined here.  See http://eigen.tuxfamily.org/ for details.
 */



#ifndef DEFS_HPP_
#define DEFS_HPP_

#include "Eigen/Dense"


using namespace Eigen;

#define complex std::complex<double>

//Definition of the matrix type using eigen
typedef Matrix<int, Dynamic, Dynamic> intMatrix;

typedef Matrix<double, Dynamic, Dynamic> DoubleMatrix;

typedef Matrix<complex, Dynamic, Dynamic> ComplexMatrix;

typedef Matrix<double, Dynamic, 1> DoubleVector;

#endif   /* DEFS_HPP_ */

