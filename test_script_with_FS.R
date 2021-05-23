### TEST CODE

#LOCAL DIR - for local machine
code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

# #MARCC DIR
# code_folder_dir <- "/home-2/hkim171@jhu.edu/scratch/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")

#SAFE DIR - for safe desktop crunchr env 
# code_folder_dir <- "~/workspace/Temporary/hkim171/scratch/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")

#SAFE only items - uncomment if using SAFE desktop crunchr environment
# install.packages("Cairo")
# options(bitmapType='cairo')

#Build ML prototype example code using example CSV files

#parameters:
feature_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_feature_space.csv"
label_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_label.csv"
outcomes =  c("Expired", "Alive")
experiment_name = "FS_test_5_23_21"
merge_identifier = "patientunitstayid"
num_outer_loop = 2

save_dir <- code_folder_dir

source(paste0(code_dir, "/model_guide.R"))
source(paste0(code_dir, "/forward_selection.R"))
source(paste0(code_dir, "/Performance_metric_plotting.R"))
source(paste0(code_dir, "/feature_ranking.R"))

#forward selection - creates a new folder with new features 
forward_selection(code_dir  = code_dir ,
                  feature_dir = feature_dir,
                  label_dir = label_dir,
                  outcomes = outcomes,
                  merge_identifier = merge_identifier,
                  relative_like_threshold = 1,
                  save_dir = save_dir,
                  experiment_name = experiment_name
)

#train model with the new FS feature space. 

fs_feature_dir = paste0(save_dir, "/", experiment_name, "_FS_results/", "FS_feature_space.csv")

build_prototype(code_dir = code_dir,
                save_dir = code_folder_dir,
                num_outer_loop = num_outer_loop,
                cross_validation_k = 3,
                merge_identifier = merge_identifier,
                feature_dir = fs_feature_dir,
                label_dir = paste0(code_dir, "/test_data/test_label.csv"),
                outcomes = outcomes,
                experiment_name = experiment_name,
                Already_trained = FALSE
)

# build_prototype(code_dir = code_dir, 
#                 save_dir = code_folder_dir,
#                 num_outer_loop = num_outer_loop,
#                 cross_validation_k = 3,
#                 merge_identifier = merge_identifier,
#                 feature_dir = fs_feature_dir,
#                 label_dir = paste0(code_dir, "/test_data/test_label.csv"),
#                 outcomes = outcomes,
#                 experiment_name = experiment_name,
#                 Already_trained = TRUE
# )


#Random forest feature ranking example
experiment_folder_dir <- save_dir
code_dir <- paste0(code_folder_dir, "/eICU_featurization")

random_forest_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = experiment_name,
                   num_outer_loop = num_outer_loop,
                   how_many_top_features = 50,
                   use_trained_rf = TRUE)

random_forest_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = experiment_name,
                   num_outer_loop = num_outer_loop,
                   how_many_top_features = 50,
                   use_trained_rf = FALSE)

XG_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = experiment_name,
                   num_outer_loop = num_outer_loop,
                   how_many_top_features = 50,
                   shap = TRUE)

XG_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = experiment_name,
                   num_outer_loop = num_outer_loop,
                   how_many_top_features = 50,
                   shap = FALSE)

GLM_rank(experiment_folder_dir <- experiment_folder_dir,
                   code_dir <- code_dir,
                   experiment_name = experiment_name,
                   num_outer_loop = num_outer_loop,
                   how_many_top_features = 50)
