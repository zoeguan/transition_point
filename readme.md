
# Scripts
- transition_point_fns.R: functions for calculating transition point/interval
- sim_fns.R: functions for running simulations
- edat_orig.R: generate cleaned list of curatedOvarianData datasets
- transition_point_thm1_tuning.R, transition_point_thm3_tuning.R: tune hyperparameters for simulation scenarios
- transition_point_thm1.R, transition_point_thm2.R, transition_point_thm3.R, transition_point_thm4.R: produce results for simulation scenarios 
- plot_thm1.R, plot_thm2.R, plot_thm3.R, plot_thm4.R: produce Figures 1-3 and 6-8
- cholesterol_single.Rnw: produce results for data application scenario 1


# Instructions for Simulations
- run edat_orig.R to generate edat_orig.RData
- run transition_point_thm1_tuning.R and transition_point_thm3_tuning.R to generate hyperparams_thm1.RData and hyperparams_thm3.RData
- run transition_point_thm1.R, transition_point_thm2.R, transition_point_thm3.R, transition_point_thm4.R to generate transition_point_thm1.RData, transition_point_thm2.RData, transition_point_thm3.RData,transition_point_thm4.RData
	- the scripts are parallelized
	- transition_point_thm3.R, transition_point_thm4.R are more time consuming and were split into parallel jobs, one for each value of sigma; the results were then combined using the script combine_thm3_thm4_results.R, which generates transition_point_thm3.RData, transition_point_thm4.RData
- run plot_thm1.R, plot_thm2.R, plot_thm3.R, plot_thm4.R to generate Figures 1-3 and 6-8


# Instructions for Data Application
- compile cholesterol_single.Rnw to generate cholesterol_single.pdf, which contains results for scenario 1 and Figure 4
- compile cholesterol_multi.Rnw to generate cholesterol_multi.pdf, which contains results for scenario 2 and Figure 5