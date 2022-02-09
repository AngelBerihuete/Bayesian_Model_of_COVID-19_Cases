// The following functions are defined in the paper
//
// A Bayesian Model of COVID-19 Cases Based on the Gompertz Curve
// 
// doi.org/10.3390/math9030228
//

functions {
  real lambdat_lpdf (vector y, real a, real M0, real c, real capital_t) {
    real logprob1;
    vector [num_elements (y)] logprob2;
    real b;
    b = log(a) - log(M0);
    logprob1  = a * exp(- b * exp(-c * capital_t)) - a * exp(-b); // Gompertz curve (eq. 2)
    for (i in 1:num_elements(y)){
      logprob2[i] = log(a) + log(b) + log(c) - b * exp(-c * y[i]) - c * y[i];
    }
    return sum(logprob2) - logprob1; // The log-likelihood (eq. 5)   
  }

  real IG_lpdf (real x, real mu, real lambda){
    real prob;
    prob = (lambda/(2*pi()*(x^3)))^0.5*exp(-lambda*(x - mu)^2/(2*mu^2*x));
    return log(prob); // The prior distribution of "a" as Inverse Gaussian (eq. 6)
    }
}

data {
  int<lower=0> N; // Total cases at time T (last date observed).
  vector[N] y;  // Arrival times assumed independent.
  real<lower=0> P; // Inhabitants per 100000 people. 
  real<lower=0> capital_t;  // T last date observed.
}

transformed data{
real upper_a;
upper_a = P*100000; // Inhabitants.
}

parameters {
  real<lower=501, upper=upper_a> a; // Upper asymptote of infections (set lower because "b" must be well-defined).
  real<lower=1, upper=500> M0; // Inital cases at time 0 in a specific region.
  real<lower=0.01, upper=0.2> c; // Coefficient of the exponential decay rate of the relative growth. 
}

model {
  target += IG_lpdf(a / P | 399.947, 525.21);   // log-prior
  target += lambdat_lpdf(y | a, M0, c, capital_t);  // log-likelihood
}

generated quantities {
  real b;
  real tmax;
  b = log(a) - log(M0); // The displacement along the time.
  tmax = log(b)/c; // The point of the expected maximum rate of increase. 
}
