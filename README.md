# Bayesian Model of COVID-19 Cases

From the paper "A Bayesian Model of COVID-19 Cases Based on the Gompertz Curve"

We perform Bayesian inference for a non-homogeneous Poisson process with an
intensity function based on the Gompertz curve.

Files include in this repo are

- demo.R : a demo to deal with data an code.
- casos_diagnostico_provincia.csv : collected data from the Spanish National
  Network for Epidemiological Surveillance ([RENAVE](https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos), by its Spanish initials)
- COVID-19_cases.stan : a Bayesian model based on the Gompertz curve.
