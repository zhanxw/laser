#include "ReadData.h"

#include <vector>
#include <string>

#include "IO.h"
#include "Logger.h"
#include "base/Utils.h"
#include "base/TypeConversion.h"

extern Logger* logger;

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
