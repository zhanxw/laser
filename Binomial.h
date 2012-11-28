#ifndef _BINOMIAL_H_
#define _BINOMIAL_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


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
  int rbinom(int n, double p) const {
    if (p < 0.0 || p > 1.0) {
      fprintf(stderr, "Wrong p");
      exit(1);
    }
    // gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)
    return (int) (gsl_ran_binomial(r, p, n));
  }
 private:
  const gsl_rng_type * T;
  gsl_rng * r;
};


#endif /* _BINOMIAL_H_ */
