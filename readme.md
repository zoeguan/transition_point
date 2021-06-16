
# Simulation Scripts
- transition_point_fns.R: functions for calculating transition point/interval
- sim_fns.R: functions for running simulations
- edat_orig.R: generate cleaned list of curatedMetagenomicData datasets
- tuning_\*.R: tune hyperparameters for simulation scenarios
- transition_point_\*.R: produce results for simulation scenarios
- plot_\*.R: produce plots for simulation scenarios

# Instructions for Simulations
- run edat_orig.R to generate edat_orig.RData
- run tuning_\*.R to generate hyperparams_\*.RData 
- run transition_point_\*.R: to generate transition_point_\*.RData
	- the scripts are parallelized
- run plot_\*.R to generate Figure 1 and Appendix Figures A.1-A.9

# Data Application Scripts
- cholesterol_single.Rnw: produce results for data application scenario 1
- cholesterol_multi.Rnw: produce results for data application scenario 2

# Instructions for Data Application
- compile cholesterol_single.Rnw to generate cholesterol_single.pdf, which contains results for scenario 1 and Figure 2
- compile cholesterol_multi.Rnw to generate cholesterol_multi.pdf, which contains results for scenario 2 and Figure 3