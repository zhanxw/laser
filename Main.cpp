#ifndef _MAIN_H_
#define _MAIN_H_

#include <Eigen/Core>
// #include <Eigen/Eigenvalues>
#include <Eigen/SelfAdjointEigenSolver.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;

#include "base/Logger.h"

void normalizeMatrix(Mat* m) {
  (*m).colwise().normalize();
}
class PCA{
 public:
  /**
   * Decompose n by n matrix @param m
   * NOTE: @param needs to be normalized
   * @return 0 if success
   */
  int compute(const Mat& m) {
    es.compute(m);
    return (es.info == Eigen::Success ? 0 : -1);
  };
  const Mat& getV() const {
    return es.eigenvectors();
  }
  const Mat& getD() const {
    return es.eigenvalues().asDiagonal();
  }

  // mat1.leftCols<cols>()
 private:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
};

class Procrustes{
 public:
  /**
   * Find a mapping f, such that minimize the distance of \sum_i (f(x)_i  - y_i)^2
   * NOTE: @param x and @param y needs normalized
   * @return 0 if success
   */
  int compute(const Mat& x, const Mat& y) {
    C = y.tranpose() * x;
    svd.compute(C, Eigen::ComputeThinU | Eigen::ComputeThinV);

    A = getV() * (getU().transpose());
    traceOfXtX = (x.transpose() * x).asDiagonal().sum();
    traceOfS = svd.getS().sum();
    rho =  traceOfS / traceOfXtX;
    // NOTE: b = 0, since we alreadyed centered @param x and @param y
    traceOfYtY = (y.transpose() * y).asDiagonal().sum();
    d = traceOfYtY - traceOfS * traceOfS / traceOfXtX;
    D = 1 - traceOfS * traceOfS / traceOfXtX / traceOfYtY;
    return 0;
  }
  /**
   * apply f(x) to @param x.
   */
  void transform(const Mat& x, Mat* y) {
    (*y) = rho * A * x;
  }
  // U * S * V_t = original matrix
  const Mat& getU() const {
    return svd.matrixU();
  }
  const Eigen::SingularValuesType & getS() const {
    return svd.singularValues();
  }
  const Mat& getV() const {
    return svd.matrixV();
  }
 private:
  Mat C;
  Mat A;
  Eigen::JacobiSVD<Mat> svd; 
  double rho;
  double d;
  double D;
  double traceOfS;
  double traceOfXtX;
  double traceOfYtY;
};

class Binomial{
 public:
  Binomial() {
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
  };
  ~Binomial() {
    gsl_rng_free (r);
  }
  int rbinom( int n, double p) const {
    // gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)
    return (int) (gsl_ran_binomail(p, n));
  }
 private:
  const gsl_rng_type * T;
  gsl_rng * r;
};

int resampleGenotype(int geno, int coverage, double errorRate, Binomial& b) {
  switch(geno) {
    case 0:
      return b.rbinom(coverage, errorRate);
    case 1:
      return b.rbinom(coverage, .5);
    case 2:
      return b.rbinom(coverage, 1.0 - errorRate);
    default:
      logger->error("Wrong genotype!\n");
      return -1;
  }
}
void resampleGenotype(const Mat& refGeno, const Vec& depth, const Vec& refCount, double errorRate, Binomial& b, Mat *resampleGeno){
  if (refGeno.cols() != depth.size() ||
      refGeno.cols() != refCount.size()) {
    return -1;
  }

  int nonEmptySite = 0;
  for (int i = 0; i < depth.size(); ++i ) {
    if (depth(i) > 0) {
      ++ nonEmptySite;
    }
  }
  
  Mat& n = *resampleGeno;
  n.resize(refGeno.rows(), nonEmptySite);

  int idx = 0;
  for (int i = 0; i < depth.size(); ++i ) {
    if (depth(i) == 0){
      continue;
    }
    for (int s = 0; s < refGeno.rows(); ++s) {
      m(s, idx) = resampleGenotype(refGeno(s,i), depth(i), errorRate);
    }
    ++idx;
  }
    int idx = 0;
      refGeno(s, i)
  }
  
};
int main(int argc, char** argv){

  // load reference coord
  Mat refCoord;
  if (readRefCoord(FLAG_refCoordFile, & refCoord)) {
    logger->error("Cannot read ref cood file\n");
  };
  // load reference geno
  Mat refGeno;
  if (readRefGeno(FLAG_refGenoFile, & refGeno)) {
    logger->error("Cannot read ref geno file\n");
  };

  if (refCoord.rows() == refGeno.rows()){
    logger->error("Ref coord and ref geno does not match");
  }
  
  // for each sample
  Mat resampleGeno;
  Mat resampleCoord;
  Vec depth;
  Vec refCount;
  std::string line;
  std::vector<std::string> fd;
  int lineNo = 0;
  LineReader lr (FLAG_seqFile);
  while (lr.readByLine(&line)){
    lineNo ++;
    //check line range
    // XXX
    
    // extract geno
    stringTokenize(line, "\t ", &fd);
    // check dimension
    if (fd.size() - 6 != refGeno.cols()) {
      logger->error("Dimension does not match at line: %d\n", lineNo);
      continue;
    }

    // resample geno
    resampleGenotype(refGeno, depth, refCount, &resampleGeno);
    
    // normalize geno
    normalizeMatrix(&resampleGeno);
    // pca resample
    pca.compute(resampleGeno);
    
    // procurstes
    Mat resampledRef = pca.getV().topLeft(refGeno.rows(), PC);
    Mat resampledNew = pca.getV().bottomLeft(1, PC);
    Mat resampledNewInOrig;

    procrustes.compute(resampleRef, refCoord);
    procrustes.transform(resampledNew, &resampledNewInOrig);

    //output
    for (int i = 0; i < 6; ++i) {
      fprintf(stdout, "%s\t", fd[i].c_str());
    }
    fprintf(stdout, procrustes.getD());
    for (int i = 0; i < PC; ++i) {
      fprintf(stdout, "\t%lf", resampleNewInOrig(0, i));
    }
    fputchar("\n", stdout);
  };
  return 0;
}

#endif /* _MAIN_H_ */
