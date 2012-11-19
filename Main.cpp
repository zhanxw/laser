#ifndef _MAIN_H_
#define _MAIN_H_

#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
// #include <Eigen/SelfAdjointEigenSolver.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;

#include "base/Argument.h"
#include "base/IO.h"
#include "base/Logger.h"
#include "base/Utils.h"

#include "GitVersion.h"

Logger* logger = NULL;

void normalizeMatrix(Mat* m) {
  for (int i = 0; i < (*m).cols(); ++i)
    (*m).col(i).normalize();
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
    V = es.eigenvectors();
    D = es.eigenvalues();
    return (es.info() == Eigen::Success ? 0 : -1);
  };
  const Mat& getV() const {
    return V;
  }
  const Vec& getD() const {
    return D;
  }

  // mat1.leftCols<cols>()
 private:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  Mat V;
  Vec D;
};

class Procrustes{
 public:
  /**
   * Find a mapping f, such that minimize the distance of \sum_i (f(x)_i  - y_i)^2
   * NOTE: @param x and @param y needs normalized
   * @return 0 if success
   */
  int compute(const Mat& x, const Mat& y) {
    C = y.transpose() * x;
    svd.compute(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // C = U * S * V^t
    U = svd.matrixU();
    S = svd.singularValues();
    V = svd.matrixV();
    
    A = getV() * (getU().transpose());
    traceOfXtX = (x.transpose() * x).trace();
    traceOfS = svd.singularValues().sum();
    rho =  traceOfS / traceOfXtX;
    // NOTE: b = 0, since we alreadyed centered @param x and @param y
    traceOfYtY = (y.transpose() * y).trace();
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
    return U;
  }
  const Vec& getS() const {
    return S; //svd.singularValues();
  }
  const Mat& getV() const {
    return V; //svd.matrixV();
  }
  double getD() const {
    return D;
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
    return (int) (gsl_ran_binomial(r, p, n));
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

int resampleGenotype(const Mat& refGeno, const Vec& depth, const Vec& refCount, double errorRate, Binomial& b, Mat *resampleGeno){
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
  
  Mat& m = *resampleGeno;
  m.resize(refGeno.rows(), nonEmptySite);

  int idx = 0;
  for (int i = 0; i < depth.size(); ++i ) {
    if (depth(i) == 0){
      continue;
    }
    for (int s = 0; s < refGeno.rows(); ++s) {
      m(s, idx) = resampleGenotype(refGeno(s,i), depth(i), errorRate, b);
    }
    ++idx;
  }
  return 0;
}

/**
 * @return a string representing current time, without '\n' at the end
 */
std::string currentTime() {
  time_t t = time(NULL);
  std::string s (ctime(&t));
  s = s.substr(0, s.size() - 1);
  return s;
};

int readIntoMatrix(const std::string& fn, int firstNumRowToSkip, int firstNumColToSkip, Mat* out) {
  Mat& m = *out;
  double d;
  
  std::string line;
  std::vector<std::string> fd;
  int lastCol = -1;
  int lineNo = 0;
  LineReader lr (fn);
  while (lr.readLine(&line)){
    lineNo ++;
    if (lineNo <= firstNumRowToSkip) {
      // line 1 is just like this:
      // Variance Percentage: 6.02522 4.52461 1.63338 1.10073 
      continue;
    }
    if (lastCol < 0 )
      lastCol == fd.size();
    if (fd.size() != lastCol) {
      logger->error("Coord file line %d has problem", lineNo);
      continue;
    }

    for (size_t i = firstNumColToSkip; i < fd.size(); ++i ){
      d = atof(fd[i]);
      m << d;
    }
  }
  int ncol = lastCol - firstNumColToSkip; // skiping first two column
  int nrow = m.size() / ncol;
  if (m.size() != ncol * nrow) {
    logger->error("Dimension does not match when reading coord files");
    return -1;
  } else {
    m.conservativeResize(nrow, ncol);
  }
  return 0;
} // end int readIntoMatrix(const std::string& fn, int firstNumRowToSkip, int firstNumColToSkip, Mat* out) {

int readRefCoord(const std::string& fn, Mat* out) {
  const int numRowToSkip = 1;
  const int numColToSkip = 2;
  return readIntoMatrix(fn, numRowToSkip, numColToSkip, out);
}

int readRefGeno(const std::string& fn, Mat* out) {
  const int numRowToSkip = 0;
  const int numColToSkip = 6;
  return readIntoMatrix(fn, numRowToSkip, numColToSkip, out);
}

int main(int argc, char** argv){
  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, refCoordFile, "--refCoord", "specify reference coordinates")
      ADD_STRING_PARAMETER(pl, refGenoFile, "--refGeno", "specify reference genotype file")
      ADD_STRING_PARAMETER(pl, seqFile, "--seqFile", "specify .seq file generated from pile-ups")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_BOOL_PARAMETER(pl, help, "--help", "Print detailed help message")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);

  if (FLAG_help) {
    pl.Help();
    return 0;
  }
  pl.Status();
  if (FLAG_REMAIN_ARG.size() > 0){
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    abort();
  }

  REQUIRE_STRING_PARAMETER(FLAG_refCoordFile, "Please provide input file using: --refCoord");
  REQUIRE_STRING_PARAMETER(FLAG_refGenoFile, "Please provide input file using: --refGeno");
  REQUIRE_STRING_PARAMETER(FLAG_seqFile, "Please provide input file using: --seqFile");  
  Logger _logger( (FLAG_outPrefix + ".log").c_str());
  logger = &_logger;
  logger->infoToFile("Program Version");
  logger->infoToFile(gitVersion);
  logger->infoToFile("Parameters BEGIN");
  pl.WriteToFile(logger->getHandle());
  logger->infoToFile("Parameters END");

  time_t startTime = time(0);
  logger->info("Analysis started at: %s", currentTime().c_str());

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
  Binomial binomial;
  Mat resampleGeno;
  Mat resampleCoord;
  Vec depth;
  Vec refCount;
  PCA pca;
  Procrustes procrustes;

  const int FLAG_PC = 4;
  const double FLAG_errorRate = 0.01;
  
  std::string line;
  std::vector<std::string> fd;
  int lineNo = 0;
  LineReader lr (FLAG_seqFile);
  while (lr.readLine(&line)){
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
    resampleGenotype(refGeno, depth, refCount, FLAG_errorRate, binomial, &resampleGeno);
    
    // normalize geno
    normalizeMatrix(&resampleGeno);
    // pca resample
    pca.compute(resampleGeno);
    
    // procurstes
    Mat resampledRef = pca.getV().topLeftCorner(refGeno.rows(), FLAG_PC);
    Mat resampledNew = pca.getV().bottomLeftCorner(1, FLAG_PC);
    Mat resampledNewInOrig;

    procrustes.compute(resampledRef, refCoord);
    procrustes.transform(resampledNew, &resampledNewInOrig);

    //output
    for (int i = 0; i < 6; ++i) {
      fprintf(stdout, "%s\t", fd[i].c_str());
    }
    fprintf(stdout, "%lf", procrustes.getD());
    for (int i = 0; i < FLAG_PC; ++i) {
      fprintf(stdout, "\t%lf", resampledNewInOrig(0, i));
    }
    fputc('\n', stdout);
  };

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);
  return 0;
}

#endif /* _MAIN_H_ */
