#ifndef _READDATA_H_
#define _READDATA_H_

#include "Common.h"

int readIntoMatrix(const std::string& fn, int firstNumRowToSkip, int firstNumColToSkip, Mat* out);
  
#endif /* _READDATA_H_ */
