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

#include "PCA.h"
#include "ReadData.h"

#include "GitVersion.h"

Logger* logger = NULL;


int main(int argc, char** argv){
  Mat x;
  readIntoMatrix("test.pca", 0, 0, &x);
  PCA pca;
  pca.compute(x);
  std::cout << x << "\n";  
  std::cout << pca.getD() << "\n";
  std::cout << pca.getV() << "\n";

  return 0;
}

#endif /* _MAIN_H_ */
