#ifndef _MAIN_H_
#define _MAIN_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "base/Argument.h"
#include "base/IO.h"
#include "base/Logger.h"
#include "base/Utils.h"
#include "base/TypeConversion.h"

#include "GitVersion.h"
#include "Common.h"
#include "PCA.h"
#include "Procrustes.h"
#include "Binomial.h"

Logger* logger = NULL;

/**
 * Normalize matrix by column
 * For missing genotype, we will keep it as zero
 * For monomorphic site, we will set all elements to zero
 */
void normalizeMatrix(Mat* m) {
  Mat& geno = *m;
  int row = geno.rows();
  int col = geno.cols();
  RowVec count(col);
  count.setZero(col);
  RowVec sum(col);
  sum.setZero(col);
  RowVec sum2(col);
  sum.setZero(col);

  // std::cerr << sum.head(10) << "\n";
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < col; ++j){
      if ( geno(i,j) >= 0) { //not missing
        sum(j) += geno(i,j);
        sum2(j) += geno(i,j) * geno(i,j);
        count(j) ++;
      }
    }
  }

  // std::cerr << sum.head(10) << "\n";
  RowVec avg = sum.array() / count.array();
  RowVec sd = ((sum2.array() - sum.array()*sum.array()/count.array()) / (count.array() - 1)).sqrt();

  // std::cerr<< avg.head(10) << "\n";
  // std::cerr<< sd.head(10) << "\n";

  for (int j = 0; j < col; ++j){
    for (int i = 0; i < row; ++i){
      if (sd(j) < 1e-6) { // monomorphic
        geno.col(j).setZero();
        continue;
      }
      if ( geno(i,j) >= 0) { //missing
        geno(i,j) = ( geno(i,j) - avg(j)  ) /sd(j);
      } else {
        geno(i,j) = 0.0;
      }
    }
  }
}

int resampleGenotype(const int geno, const int coverage, const double errorRate, const Binomial& b) {
  if (coverage <= 0) {
    fprintf(stderr, "Coverage is less than or equal to zero!");
    assert(false);
  }
  switch(geno) {
    case 0:
      return b.rbinom(coverage, errorRate);
    case 1:
      return b.rbinom(coverage, .5);
    case 2:
      return b.rbinom(coverage, 1.0 - errorRate);
    case -9:
      return -9; // for missing reference genotype,
    default:
      logger->error("Wrong genotype!\n");
      exit(1);
  }
  fprintf(stderr, "Should not reach here!");
  return -1;
}

int resampleGenotype(const Mat& refGeno, const Vec& depth, const Vec& refCount, const double errorRate, const Binomial& b, Mat *resampleGeno){
  if (refGeno.cols() != depth.size() ||
      refGeno.cols() != refCount.size()) {
    logger->error("Cannot resample genotypes because of dimension mismatch");
    return -1;
  }

  int nonEmptySite = 0;
  int d;
  for (int i = 0; i < depth.size(); ++i ) {
    d = (int) depth(i);
    if ( d > 0) {
      ++ nonEmptySite;
    } else if (d < 0) {
      fprintf(stderr, "Some depths are less than zero");
    }
  }

  // // debug code
  // Mat count;
  // count.resize(3,11);
  // count.setZero();

  Mat& m = *resampleGeno;
  m.resize(refGeno.rows() + 1, nonEmptySite);

  int idx = 0;
  int g, rg;
  for (int i = 0; i < depth.size(); ++i ) {
    d = (int) depth(i);
    if ( d == 0){
      continue;
    }
    for (int s = 0; s < refGeno.rows(); ++s) {
      g = (int) refGeno(s,i);
      if (g >= 0) {
        rg = (int) resampleGenotype( g, d, errorRate, b);
        m(s, idx) = rg;
        // if (d > 10)
        //   count(g, int(10.0 * rg / depth(i))) ++;
        // fprintf(stderr, "%d\t%d\n"
      } else {
        m(s, idx) = -9; // missing will be kept at -9
      }
    }
    m(refGeno.rows(), idx) = (int) refCount(i);
    ++idx;
  }

  // std::cout << "Resample result:" << "\n";
  // std::cout << count;
  // std::cout << "\n";
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

  std::vector<double> geno;
  std::string line;
  std::vector<std::string> fd;
  int lastCol = -1;
  int lineNo = 0;
  LineReader lr (fn);
  while (lr.readLine(&line)){
    lineNo ++;
    if (lineNo <= firstNumRowToSkip) {
      continue;
    }
    stringNaturalTokenize(line, "\t ", &fd);
    if (lastCol < 0 )
      lastCol = fd.size();
    if ((int)fd.size() != lastCol) {
      logger->error("File [ %s ] line [ %d ] has problem", fn.c_str(), lineNo);
      return -1;
    }

    for (size_t i = firstNumColToSkip; i < fd.size(); ++i ){
      if (!str2double(fd[i],&d)) {
        logger->error("Wrong genotype [ %s ] in line [ %d ] and column [ %zu ]\n", fd[i].c_str(), lineNo, i);
        geno.push_back(0.0);
      } else {
        geno.push_back(d);
      }
    }
  }
  int ncol = lastCol - firstNumColToSkip; // skiping first two column
  int nrow = geno.size() / ncol;
  if ((int)geno.size() != ncol * nrow) {
    logger->error("Dimension does not match when reading file [ %s ]", fn.c_str());
    return -1;
  }
  m.resize(nrow, ncol);
  for (int i = 0; i < nrow; ++i){
    for (int j = 0; j < ncol; ++j) {
      m(i, j) = geno[i * ncol + j];
    }
  }
  return 0;
} // end int readIntoMatrix(const std::string& fn, int firstNumRowToSkip, int firstNumColToSkip, Mat* out) {

int readRefCoord(const std::string& fn, Mat* out) {
  // line 1 is just like this:
  // Variance Percentage: 6.02522 4.52461 1.63338 1.10073
  const int numRowToSkip = 1;
  // first two column are FamID and PersonID
  const int numColToSkip = 2;
  return readIntoMatrix(fn, numRowToSkip, numColToSkip, out);
}

int readRefGeno(const std::string& fn, Mat* out) {
  const int numRowToSkip = 0;
  const int numColToSkip = 6;
  return readIntoMatrix(fn, numRowToSkip, numColToSkip, out);
}

void dumpMatrix(const char* fn, const Mat& X) {
  std::ofstream fout(fn);
  fout << X;
  fout << "\n";
  fout.close();
}


void outputHeader(FileWriter& fout, int numPC) {
  fout.printf("info_1\tinfo_2\tL\taveC\tt");
  for (int i = 0; i < numPC; ++i) {
    fout.printf("\tPC%d", i + 1);
  }
  fout.write("\n");
};

/**
 * @param t: procrustes similarity score
 */
void output(FileWriter& fout, const std::string& fam, std::string& pid, int coveredSite, double meanCov, double t, Mat& row) {
  fout.printf("%s\t%s\t", fam.c_str(), pid.c_str());
  fout.printf("%d\t", coveredSite);
  fout.printf("%g\t", meanCov);
  fout.printf("%s\n", toString(row).c_str());
};

int main(int argc, char** argv){
#if 0
  Binomial b;
  for (int i = 0; i < 100; ++i ) {
    std::cout << i << "\t" << b.rbinom(10, .99) << "\n";
  }
#else
  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, refCoordFile, "--refCoord", "specify reference coordinates")
      ADD_STRING_PARAMETER(pl, refGenoFile, "--refGeno", "specify reference genotype file")
      ADD_STRING_PARAMETER(pl, seqFile, "--seqFile", "specify .seq file generated from pile-ups")
      ADD_STRING_PARAMETER(pl, outPrefix, "--out", "output prefix")
      ADD_INT_PARAMETER_WITH_DEFAULT(pl, PC, 4, "--pc", "number of PC to use (integer: >=1, default 4)")
      ADD_DOUBLE_PARAMETER_WITH_DEFAULT(pl, errorRate, 0.01, "--errorRate", "error rate (double: [0.0, 0.5), default: 0.01 )")
      ADD_PARAMETER_GROUP(pl, "Auxilliary Functions")
      ADD_BOOL_PARAMETER(pl, debug, "--debug", "specify whether to output intermediate files")
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
    logger->error("Cannot read ref cood file");
  }else{
    logger->info("Finish read ref cood file including [ %d ] samples and [ %d ] column", refCoord.rows(), refCoord.cols());
  }
  // load reference geno
  Mat refGeno;
  if (readRefGeno(FLAG_refGenoFile, & refGeno)) {
    logger->error("Cannot read ref geno file");
  } else {
    logger->info("Finish read ref geno file including [ %d ] samples and [ %d ] column ( [ %d ] genotypes )", refGeno.rows(), refGeno.cols(), (refGeno.cols() - 6)/2);
  };

  if (refCoord.rows() != refGeno.rows()){
    logger->error("Ref coord and ref geno does not match");
    exit(1);
  }

  // if (false){
  //   logger->info("Begin PCA on original reference");
  //   PCA pca;
  //   Mat r = refGeno;
  //   normalizeMatrix(&r);
  //   pca.compute(r * r.transpose());
  //   std::ofstream ofs1( (FLAG_outPrefix + ".refGeno.V").c_str());
  //   ofs1 << pca.getV() ;
  //   ofs1.close();
  //   std::ofstream ofs2( (FLAG_outPrefix + ".refGeno.D").c_str());
  //   ofs2 << pca.getD() ;
  //   ofs2.close();
  //   logger->info("Finished PCA on original reference");
  // }
  // const int FLAG_PC = 4;
  // const double FLAG_errorRate = 0.01;

  FileWriter fout( (FLAG_outPrefix + ".out").c_str());
  outputHeader(fout, FLAG_PC);

  // for each sample
  Binomial binomial;
  Mat resampleGeno;
  Vec depth;
  Vec refCount;
  PCA pca;
  Procrustes procrustes;

  int nonEmptySite = 0;
  double meanCov = 0;
  double similarityScore = 0;
  Mat resampledNewInOrig;
  resampledNewInOrig.resize(1, FLAG_PC);



  int ret = 0;
  std::string line;
  std::vector<std::string> fd;
  int lineNo = 0;
  LineReader lr (FLAG_seqFile);
  while (lr.readLine(&line)){
    lineNo ++;
    // clean output variable
    nonEmptySite = 0;
    meanCov = 0;
    similarityScore = 0;
    resampledNewInOrig.setZero();

    //check line range
    // XXXX

    // extract geno
    stringTokenize(line, "\t ", &fd);
    // check dimension
    if ((int)fd.size() != refGeno.cols() * 2 + 6) { // 6 is the first 6 columns
      logger->error("Skip: Dimension does not match at line: %d\n", lineNo);
      continue;
    }

    // resample geno
    depth.resize(refGeno.cols());
    refCount.resize(refGeno.cols());
    for (int i = 0; i < refGeno.cols(); ++i) {
      depth[i] = atoi(fd[i*2 + 6]);
      refCount[i] = atoi(fd[i*2 + 1 + 6]);
    }
    if (resampleGenotype(refGeno, depth, refCount, FLAG_errorRate, binomial, &resampleGeno)) {
      logger->error("Resample step failed!");
      output(fout, fd[0], fd[1], nonEmptySite, meanCov, similarityScore, resampledNewInOrig);
      continue;
    } else {
      logger->info("Resampled genotype is [ %d x %d ]", resampleGeno.rows(), resampleGeno.cols());
    }
    nonEmptySite = resampleGeno.cols();
    meanCov = resampleGeno.bottomRows(1).sum() / nonEmptySite;
    
    // normalize geno and remove monomorphic
    logger->info("Center and scale matrix [ %d by %d ]", (resampleGeno).rows(), (resampleGeno).cols());
    normalizeMatrix(&resampleGeno);

    // pca resample
    logger->info("Perform PCA decomposition");

    ret = pca.compute(resampleGeno * resampleGeno.transpose());
    if (ret) {
      logger->error("PCA decomposition failed, and we saved the results.");
      output(fout, fd[0], fd[1], nonEmptySite, meanCov, similarityScore, resampledNewInOrig);
      continue;

      if (FLAG_debug) {
        Mat M = resampleGeno * resampleGeno.transpose();
        std::ofstream ofs( (FLAG_outPrefix + "." + fd[0] + ".M").c_str() );
        ofs << M;
        ofs.close();

        ofs.open( (FLAG_outPrefix + "." + fd[0] + ".resampledGeno").c_str() );
        ofs << resampleGeno;
        ofs.close();
      }
    }
    logger->info("Top PC eigenvalues: %s", toString(pca.getD()).c_str());

    // procurstes
    Mat resampledRef = pca.getV().topLeftCorner(refGeno.rows(), FLAG_PC);
    Mat resampledNew = pca.getV().bottomLeftCorner(1, FLAG_PC);
    logger->info("Resampled coord before Procrustes: %s", toString(resampledNew).c_str());

    ret = procrustes.compute(resampledRef, refCoord);
    if (ret) {
      logger->error("Procrustes analysis failed, and we saved the results.");
      output(fout, fd[0], fd[1], nonEmptySite, meanCov, similarityScore, resampledNewInOrig);
      continue;

      if (FLAG_debug) {
        std::ofstream ofs( (FLAG_outPrefix + "." + fd[0] + ".resampleCoord").c_str() );
        ofs << pca.getV().leftCols(FLAG_PC);
        ofs.close();
        ofs.open( (FLAG_outPrefix + "." + fd[0] + ".origCoord").c_str() );
        ofs << refCoord;
        ofs.close();
      }
    }
    similarityScore = procrustes.getT();
    
    procrustes.transform(resampledNew, &resampledNewInOrig);
    logger->info("Resampled coord after procrustes: %s", toString(resampledNewInOrig).c_str());

    //output
    output(fout, fd[0], fd[1], nonEmptySite, meanCov, similarityScore, resampledNewInOrig);
  };

  time_t endTime = time(0);
  logger->info("Analysis ends at: %s", currentTime().c_str());
  int elapsedSecond = (int) (endTime - startTime);
  logger->info("Analysis took %d seconds", elapsedSecond);

  fout.close();
#endif

  return 0;
}

#endif /* _MAIN_H_ */
