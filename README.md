# mbbe
 Model Based Bioequivalence
Currently hard coded to be saved in u:\fda\mbbe\modelaveraging (source.dir), and run in c:\fda\mbbe\modelaveraging (home.dir). Path to data file (U:\fda\mbbe\Modelaveraging\data_seq.csv) is harded code in the .mod files
Algorithm:
 1. read the nmodels (5 currently from the source directory (u:\fda\mbbe\modelaveraging\model1 - u:\fda\mbbe\modelaveraging\model5).
 2. Parse the data file name (U:\fda\mbbe\Modelaveraging\data_seq.csv). 
 3. Read the data file name, count the number of subjects
 4. Generate samp.size bootstrap samples (currently 100), and write to data_rep??.csv where ?? is the sample number (currently 1-100)
 5. Edit the model file to read the data_rep??.csv file, copy new model file to b?.mod in each model?/? folder where ? is the model number (currently 1-5)
 6. Run bootstrap
 7. Read the xml file from each bootstrap sample and construct a data frame with nmodels + 1 columns. The columns will contain the BIC for each model, each sample. The final column will contain an integer for which model is "best" (lowest BIC) for each sample.
 8. Simulate samp.size Monte Carlo samples using the "best" model and the parameter estimates for that sample/that model. Print out PK data from $TABLE, including treatment and sequence
 9. Run NCA on the resulting data (using PKNCA)
 10. Do TOST on NCA for each Monte Carlo simulation (not done yet)
 11. Calculate power (not done yet).
