#ifndef _MAIN_H_
#define _MAIN_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
// #include <Eigen/SelfAdjointEigenSolver.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::RowVectorXd RowVec;

#include "base/Argument.h"
#include "base/IO.h"
#include "base/Logger.h"
#include "base/Utils.h"
#include "base/TypeConversion.h"

#include "ReadData.h"
#include "Procrustes.h"

#include "GitVersion.h"

Logger* logger = NULL;

void centeringMatrix(Mat* mat) {
  Mat& m = *mat;
  m = m.rowwise() - m.colwise().mean();
}

int main(int argc, char** argv){
  Mat x;
  Mat y;
  readIntoMatrix("test.procrustes.x", 0, 0, &x);
  readIntoMatrix("test.procrustes.y", 0, 0, &y);
  std::cout << x << "\n";  
  std::cout << y << "\n";

  centeringMatrix(&x);
  centeringMatrix(&y);
  
  Procrustes procrustes;
  procrustes.compute(x, y);

  std::cout << procrustes.getA() << "\n";
  std::cout << procrustes.getRho() << "\n";
  std::cout << procrustes.getD() << "\n";
  return 0;
}

#endif /* _MAIN_H_ */
