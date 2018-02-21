# Translucent-windows-2018-code
Code to generate and analyze data for D'Andrea et al. Ecol. Lett. 2018

NoiseModel.r contains functions to create the competition kernel and run the Lotka-Volterra competition models presented in D'Andrea et al. Ecology Letters 2018. 

KmeansGap.r contains the fucntion used to measure clustering in the data. It receives a data frame with traits and abundances of each species, and returns the gap index and number of clusters in the species assemblage, along with a z-score (i.e. standardized effect size of the gap index) and p-value.

Any questions about the code can be directed to the author, Rafael D'Andrea (rdandrea@illinois.edu).
