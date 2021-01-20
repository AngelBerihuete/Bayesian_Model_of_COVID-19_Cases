functions {
  real lambdat_lpdf (vector y, real a, real M0, real c, real capital_t) {
    real logprob1;
    vector [num_elements (y)] logprob2;
    real b;
    b = log(a) - log(M0);
    logprob1  = a * exp(- b * exp(-c * capital_t)) - a * exp(-b);
    for (i in 1:num_elements(y)){
      logprob2[i] = log(a) + log(b) + log(c) - b * exp(-c * y[i]) - c * y[i];
    }
    return sum(logprob2) - logprob1;
  }

  real IG_lpdf (real x, real mu, real lambda){
    real prob;
    prob = (lambda/(2*pi()*(x^3)))^0.5*exp(-lambda*(x - mu)^2/(2*mu^2*x));
    return log(prob);
    }
}

data {
  int<lower=0> N;
  vector[N] y;  // arrival times assumed independent
  real<lower=0> P;
  real<lower=0> capital_t;  //
}

transformed data{
real upper_a;
upper_a = P*100000;
}

parameters {
  real<lower=501, upper=upper_a> a; // set lower because "b" must be well-defined
  real<lower=1, upper=500> M0;
  real<lower=0.01, upper=0.2> c;
}

model {
  // log-prior (parameters are independent!!)
  target += IG_lpdf(a / P | 399.947, 525.21);
  // log-likelihood
  target += lambdat_lpdf(y | a, M0, c, capital_t); // defined above
}

generated quantities {
  real b;
  real tmax;
  b = log(a) - log(M0);
  tmax = log(b)/c;
}
