# an example of the sim
source("cp3Github_support.R")

seed=1
n=250
K=60
interval=10
upp=26
eta_opt=c(1,-1,-1)
scenario=1
nuisance='logit'
est='IPW'

result = sim_gene(n=250, K=60, interval=10, upp=26, eta_opt=c(1,-1,-1), scenario=1, nuisance=nuisance, est=est, 
                     seed=seed)