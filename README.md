# Optimal-Treatment-Regimes-QAL
The code in this repository can be used to reproduce the results found in "Constructing optimal treatment length strategies to maximize quality-adjusted lifetimes".

# Table 1

Reproducing each row of Table 1 (Simulation Scenario 2) can be done as follows:

1. Open the folder 'SimExample' and then open the file 'SimExample.R'
2. Select model and sample size parameters based on which row you want to reproduce (ie. if trying to reproduce the first row you would set `est = 'IPW'`, `nuisance = 'logit'`, `n = 250`)
3. Set `K` and `interval` parameters based on which segment of the table you want to reproduce. If trying to reproduce the segment where (in the paper) K = 6, set `K = 60`, `upp = 26` and `interval = 10`. If trying to reproduce the segment where K = 25, set `K = 100`, `upp = 36` and `interval = 4`.
4. Set Table 1 specific parameters `scenario = 2` and `eta_opt = c(1,-1,-1)`.
5. Run SimExample.R 500 times with seeds 5001, ..., 5500.
6. Save each output with the name 'QAL{seed}.rda'. (ie. 'QAL5001.rda', 'QAL5002.rda',...)
7. Set your working directory to be the directory containing the output Rda files.
8. Open cp3Github_Table.R.
9. At bottom of file change parameters `K` and `upp` to match the those set in step 3 above.
10. Run cp3Github_Table.R and get results.
