# eICU featurization
Collection of functions to automate feature extraction
Everything to be updated by March 2021

### MARCC and general package installation instructions:
MARCC modules R/3.6.1 and gcc/5.5 should be loaded
ml R/3.6.1
ml gcc/5.5

Then in terminal (not in R) run the following from the eICU_featurization directory (cd into dir)
Rscript download_packages.R

It should download all packages and a little extra to run the files in MARCC. 

MARCC_test_script.R contains some of the test cases I've run. Check it out. This file is then sbatched and submitted to SLURM using the submit.sh bash file. It currently requests 15 minutes of 1 core on the Shared node.

### list of extraction functions
<ol>
  <li>extract_patient - extracts features from patient and hosptial table for patients of interest - last updated: 12-11-2020</li>
  <li>extract_lab - extracts lab features from lab table for patients of interest - version 0 (WIP) updated: 2-11-2021</li>
  <li>extract_diagnosis - coming soon</li>
  <li>extract_GCS - coming soon</li>
  <li>extract_GCS_timeseries - coming soon</li>
  <li>extract_medication - coming soon</li>
  <li>extract_PTS_time_domain - extracts time domain PTS features from PTS data. updated 5-19-2021</li>
  <li>extract_PTS_FFT_wavelet - coming soon</li>
  <li>extract_PTS_information_theory - coming soon</li>
  <li>extract_respiratory - coming soon</li>
  <li>extract_sofa - Jason - coming soon</li>
</ol>

### list of supervised learning functions and Others

<ol>
  <li>model_guide.R - contains automated supervised machine learning prototyping code - can be run with build_prototype - V0.2</li>
  <li>Performance_metric_plotting.R - contains CSV summary of model performance and plots of each performance metric </li>
  <li>required_custom_functions.R contains snippets of useful functions and is loaded into functions automatically upon code execution </li>
  <li>feature_ranking.R Contains feature ranking of build_prototype created RData trained models - generates both CSV rankings and plots for each of the 3 algorithms integrated into build_prototype as of March 2021 </li>
   <li>forward_selection.R Contains code to run forward selection on features to select the most optimal combination of features to use for model training. Refer to test_script_with_FS.R or forward_selection.R for example code - MAY 2021</li>
</ol>
