The main program of simulation is the file named simu_main_all.m,
which includes the all codes that produce the results of simulation part in the paper.

The main program simu_main_all.m includes six parts:
%% Part 1:  Estimate B and H for the six examples and compare with LFM.
%% Part 2:  Model selection: select the number of factors q for the six examples.
%% Part 3: Evaluate the Running time for the large-sample and high dimensional simulated data.
%% Part 4: Validate the effeciency gain by one-step update.
%% Part 5: As a reviewer suggests, we compare the two methods to 
% exert the identifiability conditions in the first step of the algorithm  by Example 1 and  Example 5
%% Part 6: Generate data used for low rank model in R programming language.

The main program of real data analysis is included in the files named real_gene.m and real_arrythmia.m.