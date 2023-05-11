---
title: "MBBE"
author: "Mark Sale"
date: "11 May 2023"
output: html_document
--- 

# Model Based Bioequivalence (MBBE)
 
Model Based Bioequivalence
Currently hard coded to be saved in u:\fda\mbbe\modelaveraging (source.dir), and run in c:\fda\mbbe\modelaveraging (home.dir). Path to data file (U:\fda\mbbe\Modelaveraging\data_seq.csv) is harded code in the .mod files
Algorithm: 

 1. Read the nmodels (5 currently from the source directory (u:\fda\mbbe\modelaveraging\model1 - u:\fda\mbbe\modelaveraging\model5).\
 2. Parse the data file name (U:\fda\mbbe\Modelaveraging\data_seq.csv). \
 3. Read the data file name, count the number of subjects.\
 4. Generate samp.size bootstrap samples (currently 100), and write to data_rep??.csv where ?? is the sample number (currently 1-100).\
 5. Edit the model file to read the data_rep??.csv file, copy new model file to b?.mod in each model?/? folder where ? is the model number (currently 1-5).\
 6. Run bootstrap.\
 7. Read the xml file from each bootstrap sample and construct a data frame with nmodels + 1 columns. The columns will contain the BIC for each model, each sample..\ Included in bootstrap sample evaluation will be identifability, as defined in "Model selection and averaging of nonlinear mixed-effect models for robust phase III dose selection". Delta.parms is the absolute fractional difference between parameters between pre and post SADDLE_REST criteria for identifability. The final column will contain an integer for which model is "best" (lowest BIC) for each sample.\
 8. Simulate samp.size Monte Carlo samples using the "best" model and the parameter estimates for that sample/that model. Print out PK data from $TABLE, including treatment and sequence.\
 9. The Monte Carlo simulation can be done using the original NONMEM estimation model, or (more likely), a user provided data set that reflects the study design of the virtual study. For example, the data set may include multiple studies, sparse data, other popultions, other doses etc. Regardless of the actual study designs, the virtual study designs can reflect a standard BE study, e.g, a replicate cross over design.\
 10. Run NCA on the resulting data (using PKNCA).\
 11. Do TOST on NCA for each Monte Carlo simulation (not done yet).\
 12. Calculate power (not done yet).\
 
Command to run MBBE is:\
run.mbbe.json(Args.json)\

Where Args.json is the path to the json file with the arguments. An example of the contents of an Args.json file is below:\
{\
"home.dir": "c:/fda/mbbe/",\
"model.source": "u:/fda/mbbe/mbbe/",\
"nmodels":        2,\
"ngroups":        4,\
"samp.size":        8,\
"Reference.groups": [        1,        2 ],\
"Test.groups": [        3,        4 ],\
"numParallel":       16,\
"crash_value":   999999,\
"nmfe.path": "c:/nm744/util/nmfe74.bat",\
"delta.parms":      0.1,\
"use_check_identifiable": true,\
"NCA.end.time":       72,\
"rndseed":        1,\
"use.simulation.data": true,\
"simulation.data.path": "U:/fda/mbbe/mbbe/data_sim.csv" \
}\

  
