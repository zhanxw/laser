#ifndef _PROCRUSTES_H_
#define _PROCRUSTES_H_

#include "Common.h"

class Procrustes{
 public:
  /**
   * Find a mapping f, such that minimize the distance of \sum_i (f(x)_i  - y_i)^2
   * NOTE: @param x and @param y needs normalized
   * @return 0 if success
   */
  int compute(const Mat& x, const Mat& y) {
    // check dimension and centered
    if (x.cols() != y.cols() ||
        x.rows() != y.rows()) {
      /* logger->error("Dimension in Procrustes analysis mismatches"); */
      return -1;
    }
    bool centered;
    xMean = x.colwise().sum().array().abs();
    xMean /= x.rows();
    centered =  (xMean.array() < 1e-6).all();
    if (! centered ) {
      /* logger->info("Matrix X is not centered"); */
      xCenter = x;
      xCenter.rowwise() -= xMean;
    } else {
      xCenter = x;
    }
    /* logger->info("Done centering matrix X"); */
      
    yMean = y.colwise().sum().array().abs();
    yMean /= y.rows();
    centered = (yMean.array() < 1e-6).all();
    if (! centered ) {
      /* logger->info("Matrix Y is not centered"); */
      yCenter = y;
      yCenter.rowwise() -= yMean;
    } else {
      yCenter = y;
    }
    /* logger->info("Done centering Matrix Y"); */
    
    C = yCenter.transpose() * xCenter;
    svd.compute(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // C = U * S * V^t
    U = svd.matrixU();
    S = svd.singularValues();
    V = svd.matrixV();

    traceOfXtX = (xCenter.transpose() * xCenter).trace();
    traceOfS = S.sum();
    traceOfYtY = (yCenter.transpose() * yCenter).trace();

    A = V * (U.transpose());
    rho =  traceOfS / traceOfXtX;
    b = yMean - rho * A.transpose() * xMean;

    d = traceOfYtY - traceOfS * traceOfS / traceOfXtX;
    D = 1 - traceOfS * traceOfS / traceOfXtX / traceOfYtY;
    t = sqrt(1-D);
    return 0;
  }
  /**
   * apply f(x) to @param x.
   */
  void transform(const Mat& x, Mat* y) {
    if (x.rows() == 1) {
      (*y) = (rho * A.transpose() * x.transpose()).transpose() + b;
    } else if (x.cols() == 1) {
      (*y) = rho * A.transpose() * x + b.transpose();
    }
  }
  // U * S * V_t = original matrix
  const Mat& getU() const {
    return U;
  }
  const Vec& getS() const {
    return S; //svd.singularValues();
  }
  const Mat& getV() const {
    return V; //svd.matrixV();
  }
  const Mat& getA() const {
    return A; 
  }
  const RowVec& getB() const{
    return b;
  }
  double getD() const {
    return D;
  };
  double getRho() const {
    return rho;
  };
  double getT() const {
    return t;
  }
  void clear() {
    t = 0;
  };
 private:
  Mat C;
  Mat A;
  double rho;
  double d;
  double D;
  double traceOfS;
  double traceOfXtX;
  double traceOfYtY;
  Eigen::JacobiSVD<Mat> svd;
  Mat U;
  Mat V;
  Vec S;
  double t;
  RowVec xMean;
  RowVec yMean;
  RowVec b;
  Mat xCenter; // centered matrix of X
  Mat yCenter; // centered matrix of Y
};


#endif /* _PROCRUSTES_H_ */
