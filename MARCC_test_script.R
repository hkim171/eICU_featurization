### MARCC TEST CODE

#MY DIR - for local machine
#code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
#code_dir <- paste0(code_folder_dir, "/eICU_featurization")


#MARCC DIR
code_folder_dir <- "/home-2/hkim171@jhu.edu/scratch/"
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

source(paste0(code_dir, "/model_guide.R"))

build_prototype(code_dir = code_dir, 
                save_dir = code_folder_dir,
                
                merge_identifier = "patientunitstayid",
                feature_dir = paste0(code_dir, "/test_data/test_feature_space.csv"),
                label_dir = paste0(code_dir, "/test_data/test_label.csv"),
                outcomes = c("Expired", "Alive"),
                experiment_name = "test",
                Already_trained = FALSE
)
