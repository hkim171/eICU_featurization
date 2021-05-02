### MARCC TEST CODE

#LOCAL DIR - for local machine
code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

# #MARCC DIR
# code_folder_dir <- "/home-2/hkim171@jhu.edu/scratch/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")

#SAFE DIR - for safe desktop crunchr env 
code_folder_dir <- "~/workspace/Temporary/hkim171/scratch/"
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

#SAFE only items - uncomment if using SAFE desktop crunchr environment
install.packages("Cairo")
options(bitmapType='cairo')

#Build ML prototype example code using example CSV files

source(paste0(code_dir, "/model_guide.R"))



build_prototype(code_dir = code_dir,
                save_dir = code_folder_dir,
                num_outer_loop = 2,
                cross_validation_k = 3,
                merge_identifier = "patientunitstayid",
                feature_dir = paste0(code_dir, "/test_data/test_feature_space.csv"),
                label_dir = paste0(code_dir, "/test_data/test_label.csv"),
                outcomes = c("Expired", "Alive"),
                experiment_name = "test_4_7_21",
                Already_trained = TRUE
)

# build_prototype(code_dir = code_dir, 
#                 save_dir = code_folder_dir,
#                 num_outer_loop = 5,
#                 cross_validation_k = 5,
#                 merge_identifier = "patientunitstayid",
#                 feature_dir = paste0(code_dir, "/test_data/test_feature_space.csv"),
#                 label_dir = paste0(code_dir, "/test_data/test_label.csv"),
#                 outcomes = c("Expired", "Alive"),
#                 experiment_name = "test",
#                 Already_trained = TRUE
# )

source(paste0(code_dir, "/Performance_metric_plotting.R"))
# performance_metric_plotting(experiment_name = "test_4_7_21", saved_file_location = paste0(code_folder_dir, "/", "test_4_7_21"))


#Random forest feature ranking example
experiment_folder_dir <- code_folder_dir
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

source(paste0(code_dir, "/feature_ranking.R"))

random_forest_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = "test_4_7_21",
                   num_outer_loop = 5,
                   how_many_top_features = 50,
                   use_trained_rf = TRUE)

random_forest_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = "test_4_7_21",
                   num_outer_loop = 5,
                   how_many_top_features = 50,
                   use_trained_rf = FALSE)

XG_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = "test_4_7_21",
                   num_outer_loop = 5,
                   how_many_top_features = 50,
                   shap = TRUE)

XG_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = "test_4_7_21",
                   num_outer_loop = 5,
                   how_many_top_features = 50,
                   shap = FALSE)

GLM_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = "test_4_7_21",
                   num_outer_loop = 5,
                   how_many_top_features = 50)
