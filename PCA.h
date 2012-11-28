#ifndef _PCA_H_
#define _PCA_H_

class PCA{
 public:
  /**
   * Decompose n by n matrix @param m
   * NOTE: each column of V = lambda^{1/2} * v and D is in decreasing order
   * NOTE: @param needs to be normalized
   * @return 0 if success
   */
  int compute(const Mat& m) {
    es.compute(m);
    D = es.eigenvalues().reverse(); // make eigen values from biggest to smallest
    V.resize(es.eigenvectors().rows(), es.eigenvectors().cols());
    for (int i = 0; i < V.cols() ; ++i) {
      V.col(i) = es.eigenvectors().col(V.cols() - 1 - i) * sqrt(D(i));
      // V.col(i) = es.eigenvectors().col(V.cols() - 1 - i) ;
    }
    return (es.info() == Eigen::Success ? 0 : -1);
  };
  const Mat& getV() const {
    return V;
  }
  const Vec& getD() const {
    return D;
  }

 private:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  Mat V;
  Vec D;
};


#endif /* _PCA_H_ */
