This repository will contain all the code for the numerical simulation presented in the paper "Utility of Compartmental Models to Test the Competing Hypotheses of Pathogen Evolution and Human Intervention."
Instructions to Understand the MATLAB Code File: 

1. Main_MonteCarlo.m file contains the simulation code for the model fitting with COVID-19 infected data across the seven epidemic waves fitted using Monte Carlo simulations (10,000) for each parameter.
2. nonparametricDistfitter.m contains the code for the probability distribution figure generated for the key epidemiological parameters for the model.
3. virulencetransmission.m conations the code for the figures (boxplot, bar) that show the relationship between different key parameters of the model, to capture the effects of host-pathogen interplay during the pandemic.
4. compare.m contains the statistical summary of the model parameters, percentage change of parameters between the consecutive waves, and cliff delta calculation.
5. contartrat.m contains the code for finding the extreme R0 values (anomaly) to use as a threshold in the sensitivity analysis.
6. The R codes kbb1.R to kbb6.R contain the code for the Bayesian estimation that we used to calculate the posterior probability of the viral evolution hypotheses.
