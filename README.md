# whisking_model

A depository for the paper: Active inference, whisking and serotonin. 

This paper is a theoretical study of observations reported in (https://www.biorxiv.org/content/10.1101/762534v3) where rodents with a serotonin knocked-out gene show
a modulation of a whisking strategy when exploring an object just based on the tactile information. Here, we show that we can reproduce this aberrant levels of serotonin through
habitual or sensory precision in the active inference framework.

The main analysis can be found in the whiskers_main.m. The model is in generate_mdp_serotonin.m. The code for showing the whisking behaviour is then in MDP_whsisking_plot.m 
and MDP_whisking_plot_comparison.m. Lastly, the code LFP_serotonin_tactile.m represents the figures of LFP data.

To run these codes, it is necessary to download an SPM file from https://www.fil.ion.ucl.ac.uk/spm/
