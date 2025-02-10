Steps of analysis: 


# Preprocessing 







# Analysis 

- Run group_behavioural_plots.Rmd, which will source group_behavioural_plots.R

This runs the group-level behavioural analyses as per (Zhang & Glascher, 2020) using a series of mixed models and t-tests where appropriate. no questionnaires are included, all analyses are ran at the group-level. 

- Run questionnaire_behavioural_plots.Rmd, which will source questionnaire_behavioural_plots.R

This runs the subsequent analyses determining the influence of questionnaire scores on the same analyses ran at the group-level in group_behavioural_plots.Rmd. Specifically, we run a series of mixed-models, select the best one, and then apply that model to each scale/sub-scale, correcting for multiple comparisons. 

We plot moderation and simple-slopes graphs, which are not new analyses per se, but visualisations of the underlying mixed-effects model results.

- Run cca.Rmd, which will source cca.R

This runs the CCA analysis using two sets of CVs 

- a questionnaire CV containing the sub-scales from the different questionnaire measures
- a behavioural CV containing task measures assessing pro-social behaviour.
- a computational CV containing model parameters assessing pro-social behaviour