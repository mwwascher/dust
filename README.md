DustMCMC2.jl contains the code used to run the MCMC algorithm described in "A mechanistic modeling and estimation framework for
environmental pathogen surveillance". The algorithm is divided into two functions: DustMCMC_cal() and DustMCMC_est(). DustMCMC_cal() is used to estimate beta given observed dust from a known number n of infectious individuals shedding. DustMCMC_est() is used to estimated an unknown number of infectious individuals n given 

DustMCMC2_cal(obs, n, nreps, filename) has the following arguments
obs - observed copies of viral RNA in 50mg of dust
n - known number of infectious individuals shedding
filename - name of the file to which the posterior samples pf beta should be saved
nreps - number of MCMC samples

DustMCMC2_est(filename, obs, nmax, nreps, filename_out) has the following arguments
filename - name of the file in which the posterior estimates of beta to be used for estimating n are stored
obs - observed copies of viral RNA in 50mg of dust
nmax - defines the prior on n to be Uniform(1,nmax)
nreps - number of MCMC samples 
filename_out - name of the file to which the posterior samples of n should be saved

The values used for calibration and estimation in the simulated data are contained in the file dust_values.txt. The values used for calibration and estimation in the real data example are contained in table 3 of the paper.
