# eICU featurization
Collection of functions to automate feature extraction
Everything to be updated by March 2021

### MARCC and general package installation instructions:
MARCC modules R/3.6.1 and gcc/5.5 should be loaded
ml R/3.6.1
ml gcc/5.5

Then in terminal (not in R) run the following from the eICU_featurization directory (cd into dir)
Rscript download_packages.R

It should download all packages and a littl extra to run the files in MARCC. 

MARCC_test_script.R contains some of the test cases I've run. Check it out. This file is then sbatched and submitted to SLURM using the submit.sh bash file. It currently requests 15 minutes of 1 core on the Shared node.

### list of extraction functions
<ol>
  <li>extract_patient - extracts features from patient and hosptial table for patients of interest - last updated: 12-11-2020</li>
  <li>extract_lab - extracts lab features from lab table for patients of interest - version 0 (WIP) updated: 2-11-2021</li>
  <li>extract_diagnosis - coming soon</li>
  <li>extract_GCS - coming soon</li>
  <li>extract_GCS_timeseries - coming soon</li>
  <li>extract_medication - coming soon</li>
  <li>extract_PTS_time_domain - coming soon</li>
  <li>extract_PTS_FFT_wavelet - coming soon</li>
  <li>extract_PTS_information_theory - coming soon</li>
  <li>extract_respiratory - coming soon</li>
  <li>extract_sofa - Jason - coming soon</li>
</ol>

### list of supervised learning functions and Others

<ol>
  <li>model_guide.R - contains supervised machine learning prototyping code - can be run with build_prototype - V0.2</li>
  <li>Performance_metric_plotting.R - contains plotting code and others - automatically runs within build_prototype</li>
  <li>required_custom_functions.R and useful_functions.R are being combined and consolidated as of 2/17/2021 </li>
</ol>
