#han
#last updated 3/16/20

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) 
library(BBmisc)
library(Hmisc)
library(caret)
library(dbscan)
library(ConsensusClusterPlus)
library(Rtsne)
library(plotly)
library(mclust)
library(ggfortify)
library(orca)
library(cluster)
library(kernlab)
library(grid)
library(gridExtra)
library(clusterSim)
require(RColorBrewer)
require(scales)
require(plotly)
require(htmlwidgets)

source("~/useful_functions.R")
functionlist <- list.files("~/eICU/eICU_feature_extract/")
for (functions in functionlist) {
  source(paste0("~/eICU/eICU_feature_extract/", functions))
}

save_dir <- "~/Clustering/4450_features/VIF_analysis_clustering_7_30_GCS//"
setwd(save_dir)

dir.create(save_dir)
dir.create(paste0(save_dir, "/figures/"))

############## SMD PLOT FUNCTION ######################
SMD_plot <- function(df, varnames, category_name, list_combination, save_directory) {
  combination_list <- list_combination
  category <- category_name
  var_names <- varnames
  df_clusters <- df
  
  for (combinations in combination_list) {
    cc <- df_clusters[which(df_clusters$cluster == combinations[1] | df_clusters$cluster == combinations[2]), ]
    cc_cluster <- cc %>% dplyr::select(cluster)
    cc <- normalize(cc[, -1], method = "standardize", range = c(-1, 1), margin = 1L, on.constant = "quiet")
    
    
    c1 <- cc[which(cc_cluster$cluster == combinations[1]), ]
    c2 <- cc[which(cc_cluster$cluster == combinations[2]), ]
    
    cmeans <- as.data.frame(cbind(colMeans(c1), colMeans(c2)))
    colnames(cmeans) <- c("c1", "c2")
    cmeans$c1_minus_c2 <- cmeans$c1 - cmeans$c2
    cmeans$c2_minus_c1 <- cmeans$c2 - cmeans$c1
    
    #plot just lab values
    var_cmeans <- cmeans[which(rownames(cmeans) %in% var_names), ]
    var_cmeans <- rownames_to_column(var_cmeans)
    # colnames(var_cmeans) <- c("vars", "c1", "c2")
    colnames(var_cmeans) <- c("vars", "c1", "c2", "c1_minus_c2", "c2_minus_c1")
    
    var_cmeans <- var_cmeans[order(var_cmeans$c1, decreasing = F), ]
    var_cmeans$vars <- factor(var_cmeans$vars, levels = unique(var_cmeans$vars))
    
    max_lim <- ceiling(max(var_cmeans[,c(2:5)]))
    
    if (nrow(var_cmeans) > 400) {
      text_size = 3.5
    }else if (nrow(var_cmeans) > 200) {
      text_size = 7
    } else if (nrow(var_cmeans) > 80) {
      text_size = 15
    } else {
      text_size = 20
    }
    
    cluster_name1 <- paste0("cluster ", combinations[1])
    cluster_name2 <- paste0("cluster ", combinations[2])
  
    
    g <- ggplot(var_cmeans) + 
      geom_point(aes(x = vars, y = c1, col = cluster_name1, group = 1), size = 1, alpha = 0.4) + 
      geom_point(aes(x = vars, y = c2, col = cluster_name2, group = 2), size = 1, alpha = 0.4) +
      geom_line(aes(x = vars, y = c1_minus_c2, col = cluster_name1, group = 1), size = 1) +
      geom_line(aes(x = vars, y = c2_minus_c1, col = cluster_name2, group = 2), size = 1) +
      geom_line(aes(x = vars, y = c1, col = cluster_name1, group = 1), size = 1, alpha = 0.5) +
      geom_line(aes(x = vars, y = c2, col = cluster_name2, group = 2), size = 1, alpha = 0.5) +
      geom_hline(yintercept = 0, alpha = 0.3, size = 1) + 
      ylim((-1 * max_lim) ,max_lim) +
      ggtitle(paste0(category,": Cluster ", combinations[1], " vs ", combinations[2], " SMD \n",test_alg, " ", test_distance, " k = ", i)) +
      ylab("Standardized variable values") +
      xlab(paste0(category, " variables")) +
      theme(plot.background = element_rect(fill = "white"), 
            panel.background = element_rect(fill = "white"), 
            panel.grid.major = element_line(size = 0.3, linetype = "dashed", color = "grey"), 
            panel.grid.minor = element_line(size = 0.2, linetype = "solid", color = "grey"), 
            legend.title=element_blank()) +
      theme(plot.title = element_text(size = 30),
            axis.text.y = element_text(size = text_size), 
            axis.text.x = element_text(size = 20), 
            legend.text = element_text(size = 25), 
            axis.title.x = element_text(size = 25),
            axis.title.y = element_text(size = 25)) +
      coord_flip() 
    ggsave(plot = g, filename = paste0(save_dir, "figures/", i, "_", test_name,"_SMD_", category, "_", combinations[1], "_vs_", combinations[2],".png"), width = 15, height = 20)
    
  }
}

#load labels
loc_label <- readRDS("~/Clustering/4450_discharge_loc_label.Rds")
loc_label$label[which(loc_label$label == "Favorable")] <- "Good"
loc_label$label[which(loc_label$label == "Unfavorable")] <- "Bad"

GCS_label <- readRDS("~/Clustering/4450_mGCS_label.Rds")
GCS_label$label[which(GCS_label$label == "Favorable")] <- "Good"
GCS_label$label[which(GCS_label$label == "Unfavorable")] <- "Bad"

surv_label <- readRDS("~/Clustering/4450_survival_label.Rds")
surv_label$label[which(surv_label$label == "Alive")] <- "Good"
surv_label$label[which(surv_label$label == "Expired")] <- "Bad"

#Clustering analyis begins here: 

# df <- readRDS("~/TBI/dec_19/24_all_features_imputed_3_23_20.Rds")
df <- readRDS("~/Clustering/4450_features/5_14_4450_PTS_LAB_DEMO_eICU.Rds")
pids <- df %>% dplyr::select(patientunitstayid)
df <- df[,-1]

VIFfeatures <- as.character(read.csv("~/Clustering/4450_features/VIF_analysis_mGCS_vifs_only_isotonic_youden/VIF_147_features_VIF_pts_demo_diag_2020-05-20.csv")[,1])

# VIFfeatures <- VIFfeatures[-grep("motor|eye|verbal", VIFfeatures, ignore.case = T)]

df <- df[, which(colnames(df) %in% VIFfeatures)]
write.csv(colnames(df), "Used_features", row.names = F)

# df <- remove_binary(df)

# df <- df %>% dplyr::select(-hospitalid, -wardid)

# binary <- which_binary(df, name_boolean = T)
# df <- binary_to_factor(df)


# autoplot(prcomp(df), variance_percentage = T)
# 
# tsne.df <- Rtsne(X = df)
# plot(tsne.df$Y)
# 
# autoplot(prcomp(df), variance_percentage = T, scale = 0)
# 
# tsne.df <- Rtsne(X = df)
# plot(tsne.df$Y)

#some clustering and imaging with cluster assignment agreement. 

#consensus clustering
number_clusters <- 12

#normalize for clustering
df_norm = normalize(df, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
# autoplot(prcomp(df_norm), variance_percentage = T, scale = 0)
# 
# tsne.df <- Rtsne(X = df_norm)
# plot(tsne.df$Y)

#CUSTOM DISTANCES ###############################################################################
my_gower = function(x){ daisy(x,metric = "gower", type = list(logratio = 3))}
my_manhattan = function(x) {dist(x, method = "manhattan")}

#################################################################################################

#CUSTOM ALGORITHMS ##############################################################################
#all it need is to spit out cluster assignment. 
dianaHook = function(this_dist, k) {
  tmp = diana(this_dist, diss = TRUE)
  assignment= cutree(tmp, k)
  return(assignment)
}

spectral = function(this_dist, k) {
  tmp <- specc(x = this_dist, centers = k, kernel = "rbfdot", diss = TRUE)
  assignment = tmp@.Data
  return(assignment)
}


#################################################################################################

df_norm_m <- as.matrix(df_norm)
df_norm_m_t <- t(df_norm_m)

# test_alg_list <- c("kmdist", "pam", "hc", "dianaHook", "spectral")
# test_alg_list <- c("kmdist", "pam", "spectral")
# test_alg_list <- c("kmdist", "pam")
test_alg_list <- c("kmdist")
# test_distance_list <- c("pearson", "euclidean", "my_gower", "my_manhattan", "binary", "canberra", "minkowski", "maximum")
# test_distance_list <- c("pearson", "euclidean", "my_gower", "minkowski", "my_manhattan", "canberra")
# test_distance_list <- c("pearson", "euclidean", "my_gower", "minkowski", "my_manhattan", "canberra")
test_distance_list <- c("my_gower")

for (test_alg in test_alg_list) {
  for (test_distance in test_distance_list) {
    
    print(paste0(test_alg, "_", test_distance))

    test_name <- paste0(test_alg, "_", test_distance, "_binary")

    title= paste0(save_dir, test_name)

    # if (test_alg != "spectral") {
    #   results = ConsensusClusterPlus(d = df_norm_m_t, maxK = number_clusters, reps = 100, pItem = 0.5, pFeature = 0.8,
    #                                  title=title, clusterAlg = test_alg, distance = test_distance, seed=120692, plot="png",
    #                                  verbose = TRUE, )
    # 
    #   saveRDS(results, file = paste0(title, "/", test_alg, "_", test_distance, "results.Rds"))
    # 
    #   error_skipper <- function(x){
    #     tryCatch(
    #       expr = {
    #         message(calcICL(results, title = title, plot = "png"))
    #         message("Successfully executed the log(x) call.")
    #       },
    #       error = function(e){
    #         message('Caught an error! but skipping')
    #         print(e)
    #       },
    #       warning = function(w){
    #         message('Caught an warning!')
    #         print(w)
    #       },
    #       finally = {
    #         message('peace out. to the next line')
    #       }
    #     )
    #   }
    # 
    #   ICL_results <- error_skipper()
    # 
    # } else if (test_distance == "pearson") {
    #   dir.create(title)
    #   results <- list()
    #   for (i in 2:number_clusters) {
    #     result <- spectral(df_norm_m, i)
    #     results[[i]] <- result
    #   }
    #   saveRDS(results, paste0(title, "/results.Rds"))
    # 
    # 
    # }
    
    results <- readRDS("kmdist_my_gower_binary/kmdist_my_gowerresults.Rds")
    
      
    for (i in 4:4) {
      if (test_alg != "spectral") {
        clusters <- rownames_to_column(as.data.frame(results[[i]][['consensusClass']]))
        colnames(clusters) <- c("ids", "cluster")
      } else {
        clusters <- as.data.frame(results[[i]])
        colnames(clusters) <- c("cluster")
      }
      
      df_clusters <- cbind(clusters, df)
      
      ######## plotting label vs cluster ##################
      
      #MOTOR GCS LABEL COMPARISON
      df_cluster_plot <- df_clusters %>% dplyr::select(cluster)
      df_cluster_plot <- cbind(pids, df_cluster_plot)
      
      df_gcs <- merge(df_cluster_plot, GCS_label, by = "patientunitstayid")
      gcs.table <- df_gcs %>% dplyr::select(cluster, label)
      gcs.table <- as.data.frame(table(gcs.table))
      
      for (i in 1:length(unique(gcs.table$cluster))) {
        temp <- gcs.table[which(gcs.table$cluster == i), ]
        total <- sum(temp$Freq)
        gcs.table$Freq[which(gcs.table$cluster == i)] <- gcs.table$Freq[which(gcs.table$cluster == i)] / total
      }
      
      A = ggplot(gcs.table, aes(x = cluster, y = Freq, fill = label)) + geom_col(position = position_fill(reverse = T)) + 
        ggtitle(paste0("Neurological proportions n = ", nrow(GCS_label))) + ylab("proportions") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
      
      #DISCHARGE LOCATION -GOS COMPARISON
      df_cluster_plot <- df_clusters %>% dplyr::select(cluster)
      df_cluster_plot <- cbind(pids, df_cluster_plot)
      
      df_gcs <- merge(df_cluster_plot, loc_label, by = "patientunitstayid")
      loc.table <- df_gcs %>% dplyr::select(cluster, label)
      loc.table <- as.data.frame(table(loc.table))
      
      for (i in 1:length(unique(loc.table$cluster))) {
        temp <- loc.table[which(loc.table$cluster == i), ]
        total <- sum(temp$Freq)
        loc.table$Freq[which(loc.table$cluster == i)] <- loc.table$Freq[which(loc.table$cluster == i)] / total
      }
      
      B = ggplot(loc.table, aes(x = cluster, y = Freq, fill = label)) + geom_col(position = position_fill(reverse = T)) + 
        ggtitle(paste0("GOS surrogate proportions n = ", nrow(loc_label))) + ylab("proportions") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 

      
      #SURVIVAL label
      df_cluster_plot <- df_clusters %>% dplyr::select(cluster)
      df_cluster_plot <- cbind(pids, df_cluster_plot)
      
      df_gcs <- merge(df_cluster_plot, surv_label, by = "patientunitstayid")
      surv.table <- df_gcs %>% dplyr::select(cluster, label)
      surv.table <- as.data.frame(table(surv.table))
      
      for (i in 1:length(unique(surv.table$cluster))) {
        temp <- surv.table[which(surv.table$cluster == i), ]
        total <- sum(temp$Freq)
        surv.table$Freq[which(surv.table$cluster == i)] <- surv.table$Freq[which(surv.table$cluster == i)] / total
      }
      
      C = ggplot(surv.table, aes(x = cluster, y = Freq, fill = label)) + geom_col(position = position_fill(reverse = T)) + 
        ggtitle(paste0("Survival proportions n = ", nrow(surv_label))) + ylab("proportions") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
      
      label_plots <- grid.arrange(A, B, C, ncol = 3, nrow = 1, top = textGrob(paste0("label proportions: ", test_name, " k = ", i, "\n"), gp = gpar(fontsize= 30, font = 8))) 
      ggsave(plot = label_plots, paste0(save_dir, "figures/", i, "_label_vs_cluster_", test_name, ".png"), width = 20, height = 8)
      ###############################################################################################################
      
      #og data
      df.pca <- prcomp(df)
      ind <- get_pca_ind(df.pca)
      pca <- as.data.frame(ind$coord)
      pca <- cbind(clusters, pca)
      
      tsne.df <- Rtsne(X = df)
      tsne.df <- tsne.df$Y
      tsne.df <- cbind(tsne.df, clusters)
      
      ## GGPLOT
      
      P <- ggplot(pca) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(cluster)), size = 2) +
        ggtitle("PCA") +
        xlab("pca dimension 1") +
        ylab("pca dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              legend.title=element_blank()) 
      Q <- ggplot(tsne.df) + geom_point(aes(x = tsne.df[,1], y = tsne.df[,2], color = as.factor(cluster)), size = 2) +
        ggtitle("Tsne") +
        xlab("tsne dimension 1") +
        ylab("tsne dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 

      
      #normalized PCA
      df_norm = normalize(df, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
      
      df_clusters <- cbind(clusters, df_norm)

      df.pca_norm <- prcomp(df_norm)
      ind <- get_pca_ind(df.pca_norm)
      df.pca_norm <- as.data.frame(ind$coord)
      df.pca_norm <- cbind(clusters, df.pca_norm)
      
      tsne.df_norm <- Rtsne(X = df_norm)
      tsne.df_norm <- tsne.df_norm$Y
      tsne.df_norm <- cbind(tsne.df_norm, clusters)
      
      ## GGPLOT
      
      R <- ggplot(df.pca_norm) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(cluster)), size = 2) +
        ggtitle("Normalized PCA") +
        xlab("pca dimension 1") +
        ylab("pca dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
      S <- ggplot(tsne.df_norm) + 
        geom_point(aes(x = tsne.df_norm[,1], y = tsne.df_norm[,2], color = as.factor(cluster)), size = 2) +
        ggtitle("Normalized Tsne") +
        xlab("tsne dimension 1") +
        ylab("tsne dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
      
      #standardized PCA
      df_std = normalize(df, method = "standardize", range = c(-1, 1), margin = 1L, on.constant = "quiet")
      
      df_clusters <- cbind(clusters, df_std)
      
      df.pca_std <- prcomp(df_std)
      ind <- get_pca_ind(df.pca_std)
      df.pca_std <- as.data.frame(ind$coord)
      df.pca_std <- cbind(clusters, df.pca_std)
      
      tsne.df_std <- Rtsne(X = df_std)
      tsne.df_std <- tsne.df_std$Y
      tsne.df_std <- cbind(tsne.df_std, clusters)
      
      ## GGPLOT
      
      Y <- ggplot(df.pca_std) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(cluster)), size = 2) +
        ggtitle("Standardized PCA") +
        xlab("pca dimension 1") +
        ylab("pca dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
      Z <- ggplot(tsne.df_std) + 
        geom_point(aes(x = tsne.df_std[,1], y = tsne.df_std[,2], color = as.factor(cluster)), size = 2) +
        ggtitle("Normalized Tsne") +
        xlab("tsne dimension 1") +
        ylab("tsne dimension 2") +
        theme(plot.title = element_text(size = 25),
              axis.text.y = element_text(size = 10), 
              axis.text.x = element_text(size = 10), 
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15), 
              legend.title=element_blank()) 
        
      
      plots <- grid.arrange(P, R, Y, Q, S, Z, ncol = 3, nrow = 2, top = textGrob(paste0(test_name, " k = ", i, "\n"), gp = gpar(fontsize= 30, font = 8))) 
      ggsave(plot = plots, paste0(save_dir, "figures/", i, "_PCA_TSNE_", test_name, ".png"), width = 22.5, height = 15)
      ###############################################################################################################
      
      ### cluster feature figures. 
      df_clusters <- cbind(clusters, df)
      
      cluster_count <- length(unique(df_clusters$cluster))
      df_clusters <- df_clusters %>% dplyr::select(-ids)
      
      temp <- t(combn(seq(from = 1, to = cluster_count, by = 1), 2))
      combination_list <- split(temp, seq(nrow(temp)))

      varnames <- colnames(df_clusters)
      
      #lab
      # SMD_plot(df_clusters, varnames[2:94], "Lab_demo_gcs", combination_list, save_dir)
      # 
      # # varnames <- varnames[-which(varnames %in% lab_names)]
      # 
      # # pts_names <- readRDS("~/Clustering/ptsnames.Rds")
      # SMD_plot(df_clusters, varnames[grep("hr", varnames)], "hr", combination_list, save_dir)
      # SMD_plot(df_clusters, varnames[grep("sao2", varnames)], "sao2", combination_list, save_dir)
      # SMD_plot(df_clusters, varnames[grep("resp", varnames)], "resp", combination_list, save_dir)
      # SMD_plot(df_clusters, varnames[grep("idias", varnames)], "dias", combination_list, save_dir)
      # SMD_plot(df_clusters, varnames[grep("isys", varnames)], "sys", combination_list, save_dir)
      # SMD_plot(df_clusters, varnames[grep("imean", varnames)], "mean", combination_list, save_dir)
      
      # varnames <- varnames[-which(varnames %in% pts_names)]
      
      # med_names <- readRDS("~/Clustering/mednames.Rds")
      # SMD_plot(df_clusters, varnames[350:390], "Med", combination_list, save_dir)
      # 
      # # varnames <- varnames[-which(varnames %in% med_names)]
      # 
      # # demo_names <- readRDS("~/Clustering/demonames_yes_binary.Rds")
      # # demo_names <- readRDS("~/Clustering/demonames_no_binary.Rds")
      # SMD_plot(df_clusters, varnames[c(2:12, 310:349)], "Demo", combination_list, save_dir)
      # 
      # SMD_plot(df_clusters, varnames[c(391:435)], "Diagnosis", combination_list, save_dir)
      # 
      # SMD_plot(df_clusters, varnames[c(92:99)], "GCS", combination_list, save_dir)
      # 
      SMD_plot(df_clusters, varnames, "ALL", combination_list, save_dir)
      
      # med_names <- varnames[334:374]
      # SMD_plot(df_clusters, med_names, "Med", combination_list, save_dir)
      # varnames <- varnames[-which(varnames %in% med_names)]
      
      
      ###### normalized data only: plotting outcome proportions vs overlaid PCA/TSNE plot for EICU. #####
      #
      label_list <- list("GCS" = GCS_label, "GOS" = loc_label, "MORT" = surv_label)
      label_names <- names(label_list)
      
      
      eicu_df <- cbind(pids, df)
      
      plot_dim_redux_prop <- function(df, label_list, cluster_list, label_names, type) {
        temp_count = 1
        for (selected_label in label_list) {
          temp_label_name <- label_names[temp_count]
          
          # #gives cluster pids. 
          # pid_clusters <- cbind(pids, clusters[,-1, drop = F])
          
          #normalize the df to make normalized PCA and Tsne plots. 
          df_norm = normalize(df[,-1], method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
          # #give dataframe pids. 
          df_norm <- cbind(df[,1, drop = F], df_norm)
          
          #merge with cluster
          df_norm <- cbind("cluster" = cluster_list, df_norm)
          
          #merge with label
          df_norm <- merge(selected_label, df_norm, by = "patientunitstayid")
          
          #set aside subset label, cluster, and pids
          templabel <- df_norm$label
          tempcluster <- df_norm$cluster
          temppids <- df_norm$patientunitstayid
          
          #grab just the features
          df_norm <- df_norm %>% dplyr::select(-patientunitstayid, -label, -cluster)
          
          #pca
          df.pca_norm <- prcomp(df_norm)
          ind <- get_pca_ind(df.pca_norm)
          df.pca_norm <- as.data.frame(ind$coord)
          #add back cluster and label
          df.pca_norm <- cbind("clusters" =tempcluster, "label" = templabel, df.pca_norm)
          
          #Tsne
          tsne.df_norm <- Rtsne(X = df_norm)
          tsne.df_norm <- tsne.df_norm$Y
          colnames(tsne.df_norm) <- c("Dim.1", "Dim.2")
          tsne.df_norm <- as.data.frame(tsne.df_norm)
          #add back cluster and label
          tsne.df_norm <- cbind(tsne.df_norm, "clusters" = tempcluster, "label" = templabel)
          
          #3D TSNE
          tsne_3 <- Rtsne::Rtsne(df_norm, check_duplicates = FALSE, dim = 3)$Y
          colnames(tsne_3) <- c("Dim.1", "Dim.2", "Dim.3")
          tsne_3 <- as.data.frame(tsne_3)
          tsne_3 <- cbind(tsne_3, "clusters" = tempcluster, "label" = templabel)
          
          #PCA and TSNE plots with LABEL overlay
          
          PCA <- ggplot(df.pca_norm, aes(x = as.numeric(as.character(Dim.1)), y = as.numeric(as.character(Dim.2)), color = as.factor(label), fill = as.factor(label))) + geom_point(shape = 1, size = 4) + geom_point(size = 4, alpha = 0.6) +
            ggtitle(paste0("PCA: Label = ", temp_label_name, ", ", test_name,"\nn = ", nrow(df_norm), " f = ", ncol(df_norm))) +
            xlab("pca dimension 1") +
            ylab("pca dimension 2") +
            theme(plot.title = element_text(size = 25),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank()) 
          
          TSNE <- ggplot(tsne.df_norm, aes(x = as.numeric(as.character(Dim.1)), y = as.numeric(as.character(Dim.2)), color = as.factor(label), fill = as.factor(label))) + geom_point(shape = 1, size = 4) + geom_point(size = 4, alpha = 0.6) +
            ggtitle(paste0("TSNE: Label = ", temp_label_name, ", ", test_name,"\nn = ", nrow(df_norm), " f = ", ncol(df_norm))) +
            xlab("Tsne dimension 1") +
            ylab("Tsne dimension 2") +
            theme(plot.title = element_text(size = 25),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank()) 
          
          save_ <- paste0(save_dir,"/PCA_TSNE_plots_" ,temp_label_name,"_label_", type, ".png")
          png(save_, width = 16, height = 10, units = "in", res = 300)
          grid.arrange(PCA, TSNE, ncol = 2)
          dev.off()
          
  
          #PCA overlay. 
          temppca <- df.pca_norm
          # levels(temppca$label)
          temppca$clusters_label <- paste0("cluster-", as.character(temppca$clusters), "-", temppca$label)
          
          #TSNE overlay
          temptsne <- tsne.df_norm
          temptsne$clusters_label <- paste0("cluster-", as.character(temptsne$clusters), "-", temptsne$label)
          
          cluster_names <- paste0("cluster-", seq(1, i, 1))
          goodvarnames <- paste0(cluster_names, "-Good")
          badvarnames <- paste0(cluster_names, "-Bad")
          
          colors <- c()
          for (j in 1:i) {
            # colors <- c(colors, hue_pal()(10)[j])   
            colors <- c(colors, brewer.pal(name = "Set1", n = 8)[j]) 
          }
          
          alpha <- c()
          possible_alphas <- seq(from = 1, to = .5, length.out = length(unique(selected_label$label)))
          shape <- c()
          possible_shapes <- c(4, 1) #x, circle
          for (j in 1:length(unique(selected_label$label))) {
            alpha <- c(alpha, rep(possible_alphas[j], i))
            shape <- c(shape, rep(possible_shapes[j], i))
          }
          
          
          varnames <- tibble(
            label = c(badvarnames,goodvarnames),
            Color = rep(colors, length(unique(selected_label$label))),
            alpha = alpha, 
            shape = shape
          )
          
          varnames <- varnames[order(varnames$label), ]
          
          temppca$clusters_label <- as.factor(temppca$clusters_label)
          temppca$clusters_label <- factor(temppca$clusters_label, levels = varnames$label)
          
          title <- paste0("PCA: ", type, " Cluster k = ", i, " for label = ", temp_label_name, " ---", test_name)
          p <- ggplot(temppca) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(clusters_label), alpha = as.factor(clusters_label), shape = as.factor(clusters_label), fill = as.factor(clusters_label)), size = 1.8) +
            ggtitle(title) +
            xlab("pca dimension 1") +
            ylab("pca dimension 2") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank(),
                  legend.position = "none"
            ) + scale_color_manual(values = varnames$Color) +
            scale_alpha_manual(values = varnames$alpha) + 
            scale_shape_manual(values = varnames$shape ) +
            scale_fill_manual(values = varnames$Color)
          
          p_l <- ggplot(temppca) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(clusters_label), alpha = as.factor(clusters_label), shape = as.factor(clusters_label), fill = as.factor(clusters_label)), size = 1.8) +
            ggtitle(title) +
            xlab("pca dimension 1") +
            ylab("pca dimension 2") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank()
                  # legend.position = "none") 
            ) + scale_color_manual(values = varnames$Color) +
            scale_alpha_manual(values = varnames$alpha) + 
            scale_shape_manual(values = varnames$shape ) +
            scale_fill_manual(values = varnames$Color)
          
          title <- paste0("TSNE: ", type, " Cluster k = ", i, " for label = ", temp_label_name, " ---", test_name)
          ts <- ggplot(temptsne) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(clusters_label), alpha = as.factor(clusters_label), shape = as.factor(clusters_label), fill = as.factor(clusters_label)), size = 1.8) +
            ggtitle(title) +
            xlab("pca dimension 1") +
            ylab("pca dimension 2") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank(),
                  legend.position = "none"
            ) + scale_color_manual(values = varnames$Color) +
            scale_alpha_manual(values = varnames$alpha) + 
            scale_shape_manual(values = varnames$shape ) +
            scale_fill_manual(values = varnames$Color)
          
          ts_l <- ggplot(temptsne) + geom_point(aes(x = Dim.1, y = Dim.2, color = as.factor(clusters_label), alpha = as.factor(clusters_label), shape = as.factor(clusters_label), fill = as.factor(clusters_label)), size = 1.8) +
            ggtitle(title) +
            xlab("pca dimension 1") +
            ylab("pca dimension 2") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank()
                  # legend.position = "none") 
            ) + scale_color_manual(values = varnames$Color) +
            scale_alpha_manual(values = varnames$alpha) + 
            scale_shape_manual(values = varnames$shape ) +
            scale_fill_manual(values = varnames$Color)
          
          # plotly_pca <- plotly::plot_ly(temppca, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, color = ~as.factor(clusters), size = 2) %>% layout(title = paste0("3D PCA: ",  test_name, " ", temp_label_name,"\nn = ", nrow(df_norm), " f = ", ncol(df_norm), " k = ", i))
          # htmlwidgets::saveWidget(plotly_pca, paste0(save_dir,"/3D_PCA_k_",i,"_", test_name,"_",temp_label_name,".html")) 
        
          #### OUTCOME PROPORTION PLOTS
          cluster.table <- df.pca_norm %>% dplyr::select(clusters, label)
          
          if (temp_label_name == "GCS") {
            levels(cluster.table$label) <- c("Unfavorable", "Favorable")
          } else if (temp_label_name == "MORT") {
            levels(cluster.table$label) <- c("Expired", "Alive")
          }
          
          cluster.table$clusters <- as.factor(cluster.table$clusters)
          levels(cluster.table$clusters) <- c("b", "c", "a", "d")
          cluster.table$clusters <- relevel(cluster.table$clusters, ref = "d")
          cluster.table$clusters <- relevel(cluster.table$clusters, ref = "c")
          cluster.table$clusters <- relevel(cluster.table$clusters, ref = "b")
          cluster.table$clusters <- relevel(cluster.table$clusters, ref = "a")
          
          
          cluster.table <- as.data.frame(table(cluster.table))
          
          cluster.table$Prop <- NA
          for (cluster_character in (unique(cluster.table$clusters))) {
            temp <- cluster.table[which(cluster.table$clusters == cluster_character), ]
            total <- sum(temp$Freq)
            cluster.table$Prop[which(cluster.table$clusters == cluster_character)] <- cluster.table$Freq[which(cluster.table$clusters == cluster_character)] / total
          }
          
          outcome_clusters <- unique(cluster.table$clusters)
          cluster.table$new_prop <- NA
          new_cluster_table <- rbind()
          for (j in 1:length(outcome_clusters)) {
            temp <- cluster.table[cluster.table$clusters == outcome_clusters[j], ]
            temp <- temp[order(-temp$Prop),]
            for (k in 1:nrow(temp)) {
              temp$new_prop[k] <- sum(temp$Prop[c(k:nrow(temp))], na.rm = T)
            }
            new_cluster_table <- rbind(new_cluster_table, temp)
          }
          
          new_cluster_table$cluster_label <- paste0("cluster-", new_cluster_table$cluster, "-", new_cluster_table$label)
     
          if(temp_label_name == "GCS") {
            temp_variable <- "mGCS"
          } else if (temp_label_name == "MORT") {
            temp_variable <- "Mortality"
          } else if (temp_label_name == "GOS") {
            temp_variable <- "sCPC"
          }
          
          cluster_names <- paste0("cluster-", unique(new_cluster_table$clusters))
          
          if(temp_label_name == "GCS") {
            goodvarnames <- paste0(cluster_names, "-Favorable")
            badvarnames <- paste0(cluster_names, "-Unfavorable")
          } else if (temp_label_name == "MORT") {
            goodvarnames <- paste0(cluster_names, "-Alive")
            badvarnames <- paste0(cluster_names, "-Expired")
          } else if (temp_label_name == "GOS") {
            goodvarnames <- paste0(cluster_names, "-Favorable")
            badvarnames <- paste0(cluster_names, "-Unfavorable")
          }
          
          colors <- c()
          for (j in 1:i) {
            # colors <- c(colors, hue_pal()(10)[j])   
            colors <- c(colors, brewer.pal(name = "Set1", n = 8)[j]) 
          }
          
          alpha <- c()
          possible_alphas <- seq(from = 1, to = .5, length.out = length(unique(selected_label$label)))
          shape <- c()
          possible_shapes <- c(4, 1) #x, circle
          for (j in 1:length(unique(selected_label$label))) {
            alpha <- c(alpha, rep(possible_alphas[j], i))
            shape <- c(shape, rep(possible_shapes[j], i))
          }
          
          
          varnames_pp <- tibble(
            label = c(badvarnames,goodvarnames),
            Color = rep(colors, length(unique(selected_label$label))),
            alpha = alpha, 
            shape = shape
          )
          
          varnames_pp <- varnames_pp[order(varnames_pp$label), ]
          
               
          pp = ggplot(new_cluster_table, aes(fill = as.factor(cluster_label), x = clusters, y = Freq, alpha = as.factor(cluster_label))) + geom_bar(position = position_fill(reverse = F), stat="identity") + 
            ggtitle(paste0(type, " ",temp_variable, " Proportions")) + ylab("proportions") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank(),
                  legend.position = "none") +
            # geom_text(stat='identity',aes(label = Freq), vjust=20) 
            # ylim(0, 1)
            # geom_text(
            #   aes(y = new_prop, label=Freq),
            #   stat='identity',
            #   size=4, 
            #   vjust = 2
            # ) + 
            scale_fill_manual(values = varnames_pp$Color) +
            scale_alpha_manual(values = varnames_pp$alpha)
          
          pp_l = ggplot(new_cluster_table, aes(fill = as.factor(cluster_label), x = clusters, y = Freq, alpha = as.factor(cluster_label))) + geom_bar(position = position_fill(reverse = F), stat="identity") + 
            ggtitle(paste0(type, " ",temp_variable, " Proportions")) + ylab("proportions") +
            theme(plot.title = element_text(size = 15),
                  axis.text.y = element_text(size = 10), 
                  axis.text.x = element_text(size = 10), 
                  axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15), 
                  legend.title=element_blank()) +
            # geom_text(stat='identity',aes(label = Freq), vjust=20) 
            # ylim(0, 1)
            # geom_text(
            #   aes(y = new_prop, label=Freq),
            #   stat='identity',
            #   size=4, 
            #   vjust = 2
            # ) + 
            scale_fill_manual(values = varnames_pp$Color) +
            scale_alpha_manual(values = varnames_pp$alpha)
          
          save_ <- paste0(save_dir,"/PCA_PP_", temp_label_name, "_k_", i , "_", test_name, "_", type,".png")
          png(save_, width = 10, height = 4, units = "in", res = 400)
          grid.arrange(p, pp_l, ncol = 2)
          dev.off()
            
          save_ <- paste0(save_dir,"/TSNE_PP_", temp_label_name, "_k_", i , "_", test_name, "_", type,".png")
          png(save_, width = 10, height = 4, units = "in", res = 400)
          grid.arrange(ts, pp_l, ncol = 2)
          dev.off()
  
          
          temp_count = temp_count + 1
        }
        
      }
      
      plot_dim_redux_prop(eicu_df, label_list, clusters$cluster, label_names, type = "eICU")
      
      
      ###### multiclass classification + all plots associated with MIMIC vs EICU
      eicu_df <- cbind("cluster" = clusters$cluster, df)
      eicu_df$cluster <- paste0("cluster", eicu_df$cluster) 
      eicu_df$cluster <- as.factor(eicu_df$cluster)
      
      control = trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=3,
                             summaryFunction = multiClassSummary,
                             classProbs = T,
                             savePredictions = T)
      
      # modelxgboost <- caret::train(cluster~., data = eicu_df, method="xgbTree", metric = 'ROC',trControl=control)
      # save(modelxgboost,  file = paste('_modelXG_num_cluster_',i, '.Rdata', sep = ''))
      load(file = paste('_modelXG_num_cluster_',i, '.Rdata', sep = ''))
      model.predicted.prob_xgboost <- predict(modelxgboost, newdata = eicu_df, type = "raw")
      
      sink(paste0('_modelXG_num_cluster_',i, ".txt"))
      confusionMatrix(model.predicted.prob_xgboost, as.factor(eicu_df$cluster))
      
      #load in MIMIC as testing. 
      mimic_df <- read.csv("~/Clustering/4450_features/5_14_2023_PTS_LAB_DEMO_MIMIC_imputed.csv")
      mimic_df <- mimic_df %>% dplyr::select("patientunitstayid", all_of(VIFfeatures[-1])) 
      
      mimic_clusters <- predict(modelxgboost, newdata = mimic_df, type = "raw")
      mimic_clusters <- as.integer(gsub('[a-zA-Z]', '', mimic_clusters))
      eicu_df$cluster <- as.integer(gsub('[a-zA-Z]', '', eicu_df$cluster))
      
      #load all mimic labels. 
      m_surv_label <- read.csv("~/Clustering/4450_features/MIMIC_mort_label.csv")
      levels(m_surv_label$label) <- c("Good", "Bad")
      m_surv_label$label <- as.character(m_surv_label$label)
      
      m_gcs_label <- read.csv("~/Clustering/4450_features/MIMIC_GCS_label.csv")
      m_gcs_label$label[which(m_gcs_label$label < 6)] <- "Bad"
      m_gcs_label$label[-which(m_gcs_label$label == "Bad")] <- "Good"
      m_gcs_label$label <- as.character(m_gcs_label$label)
      
      m_loc_label <- read.csv("~/Clustering/4450_features/MIMIC_location_label.csv")
      levels(m_loc_label$label) <- c("Good", "Bad")
      m_loc_label$label <- as.character(m_loc_label$label)
      
      label_list <- list("GCS" = m_gcs_label, "GOS" = m_loc_label, "MORT" = m_surv_label)
      label_names <- names(label_list)
      
      plot_dim_redux_prop(mimic_df, label_list, mimic_clusters, label_names, type = "MIMIC")
      
      mimic_df <- cbind("cluster" = mimic_clusters, mimic_df)
    
      
      SMD_plot_two_types <- function(df1, df2, type1, type2, varnames, category_name, list_combination, save_directory, percent_from_each_end) {
        require(BBmisc)
        ratio <- percent_from_each_end/100
        num_features <- length(varnames) - 1
        feature_from_each_end <- ceiling(num_features * ratio)
        
        combination_list <- list_combination
        category <- category_name
        var_names <- varnames
        
        df1$type = type1
        df2$type = type2
        types <- c(type1, type2)
        
        df_clusters <- rbind(df1, df2)
        
        df_clusters$cluster <- as.factor(df_clusters$cluster)
        levels(df_clusters$cluster) <- c('b', 'c', 'a', 'd')
        df_clusters$cluster <- as.character(df_clusters$cluster)
        
        numclusters <- unique(df_clusters$cluster)
        
        cluster_names <- numclusters[order(numclusters)]

        colors <- c()
        for (j in 1:length(numclusters)) {
          # colors <- c(colors, hue_pal()(10)[j])   
          colors <- c(colors, brewer.pal(name = "Set1", n = 8)[j]) 
        }

        
        color_table <- tibble(
          label = cluster_names,
          Color = colors
        )
        
        for (combinations in combination_list) {
          
          final_var_cmeans <- as.data.frame(var_names[-1])
          colnames(final_var_cmeans) <- "vars"
          for (type in types) {
            tempdf_df_clusters <- df_clusters[which(df_clusters$type == type), ]
            cc <- tempdf_df_clusters[which(tempdf_df_clusters$cluster == cluster_names[combinations[1]] | tempdf_df_clusters$cluster == cluster_names[combinations[2]]), ]
            cc_cluster <- cc %>% dplyr::select(cluster)
            cc <- cc %>% dplyr::select(-cluster)
            cc <- normalize(cc, method = "standardize", range = c(-1, 1), margin = 1L, on.constant = "quiet")
            
            
            c1 <- cc[which(cc_cluster$cluster == cluster_names[combinations[1]]), ]
            c1_type <- c1$type
            c1 <- c1 %>% dplyr::select(-type)
            c2 <- cc[which(cc_cluster$cluster == cluster_names[combinations[2]]), ]
            c2_type <- c2$type
            c2 <- c2 %>% dplyr::select(-type)
            
            cmeans <- as.data.frame(cbind(colMeans(c1), colMeans(c2)))
            colnames(cmeans) <- c("c1", "c2")
            cmeans$c1_minus_c2 <- cmeans$c1 - cmeans$c2
            cmeans$c2_minus_c1 <- cmeans$c2 - cmeans$c1
            
            #plot just lab values
            var_cmeans <- cmeans[which(rownames(cmeans) %in% var_names), ]
            var_cmeans <- rownames_to_column(var_cmeans)
            # colnames(var_cmeans) <- c("vars", "c1", "c2")
            tempnames <- c("c1", "c2", "c1_minus_c2", "c2_minus_c1")
            tempnames <- paste0(type, "_", tempnames)
            colnames(var_cmeans) <- c("vars", tempnames)
            
            final_var_cmeans <- merge(final_var_cmeans, var_cmeans, by = "vars")
            
          }
          
          final_var_cmeans <- final_var_cmeans[order(final_var_cmeans[,2], decreasing = F), ]
          final_var_cmeans$vars <- factor(final_var_cmeans$vars, levels = unique(final_var_cmeans$vars))
          
          final_var_cmeans <- final_var_cmeans[c(1:feature_from_each_end,( nrow(final_var_cmeans) - feature_from_each_end):nrow(final_var_cmeans)),]
          
          final_var_cmeans <- final_var_cmeans[order(final_var_cmeans[,2], decreasing = F), ]
          final_var_cmeans$vars <- factor(final_var_cmeans$vars, levels = unique(final_var_cmeans$vars))
          
          max_lim <- ceiling(max(final_var_cmeans[,c(2:9)]))
          
          if (nrow(final_var_cmeans) > 400) {
            text_size = 3.5
          }else if (nrow(final_var_cmeans) > 200) {
            text_size = 7
          } else if (nrow(final_var_cmeans) > 80) {
            text_size = 15
          } else {
            text_size = 20
          }
          
          cluster_name1a <- paste0(type1, " endotype ", cluster_names[combinations[1]])
          cluster_name2a <- paste0(type1, " endotype ", cluster_names[combinations[2]])
          cluster_name1b <- paste0(type2, " endotype ", cluster_names[combinations[1]])
          cluster_name2b <- paste0(type2, " endotype ", cluster_names[combinations[2]])
          
          cluster_1_color <- color_table$Color[which(color_table$label == cluster_names[combinations[1]])]
          cluster_2_color <- color_table$Color[which(color_table$label == cluster_names[combinations[2]])]
          
          
          g <- ggplot(final_var_cmeans) + 
            geom_point(aes(x = vars, y = eICU_c1, col = cluster_name1a, group = 1), size = 1, alpha = 0.4, color = cluster_1_color) + 
            geom_point(aes(x = vars, y = eICU_c2, col = cluster_name2a, group = 2), size = 1, alpha = 0.4, color = cluster_2_color) +
            geom_line(aes(x = vars, y = eICU_c1_minus_c2, col = cluster_name1a, group = 1), size = 1, color = cluster_1_color) +
            geom_line(aes(x = vars, y = eICU_c2_minus_c1, col = cluster_name2a, group = 2), size = 1, color = cluster_2_color) +
            geom_line(aes(x = vars, y = eICU_c1, col = cluster_name1a, group = 1), size = 1, alpha = 0.5, color = cluster_1_color) +
            geom_line(aes(x = vars, y = eICU_c2, col = cluster_name2a, group = 2), size = 1, alpha = 0.5, color = cluster_2_color) +
            
            geom_point(aes(x = vars, y = MIMIC_c1, col = cluster_name1b, group = 1), size = 1, alpha = 0.4, color = cluster_1_color) + 
            geom_point(aes(x = vars, y = MIMIC_c2, col = cluster_name2b, group = 2), size = 1, alpha = 0.4, color = cluster_2_color) +
            geom_line(aes(x = vars, y = MIMIC_c1_minus_c2, col = cluster_name1b, group = 1), size = 1, linetype = "twodash", color = cluster_1_color) +
            geom_line(aes(x = vars, y = MIMIC_c2_minus_c1, col = cluster_name2b, group = 2), size = 1, linetype = "twodash", color = cluster_2_color) +
            geom_line(aes(x = vars, y = MIMIC_c1, col = cluster_name1b, group = 1), size = 1, alpha = 0.5, linetype = "twodash", color = cluster_1_color) +
            geom_line(aes(x = vars, y = MIMIC_c2, col = cluster_name2b, group = 2), size = 1, alpha = 0.5, linetype = "twodash", color = cluster_2_color) +
            
            geom_hline(yintercept = 0, alpha = 0.3, size = 1) + 
            ylim((-1 * max_lim) ,max_lim) +
            ggtitle(paste0("Standardized Mean Difference\neICU vs MIMIC external validation: \nEndotype ", cluster_names[combinations[1]], " vs ", cluster_names[combinations[2]])) +
            ylab("Standardized variable values") +
            xlab(paste0("variables")) +
            theme(plot.background = element_rect(fill = "white"), 
                  panel.background = element_rect(fill = "white"), 
                  panel.grid.major = element_line(size = 0.3, linetype = "dashed", color = "grey"), 
                  panel.grid.minor = element_line(size = 0.2, linetype = "solid", color = "grey"), 
                  legend.title=element_blank()) +
            theme(plot.title = element_text(size = 30),
                  axis.text.y = element_text(size = text_size), 
                  axis.text.x = element_text(size = 20), 
                  legend.text = element_text(size = 25), 
                  axis.title.x = element_text(size = 25),
                  axis.title.y = element_text(size = 25)) +
            coord_flip() 
          ggsave(plot = g, filename = paste0(save_dir, "/", i, "_", test_name,"_SMD_", category, "_", cluster_names[combinations[1]], "_vs_", cluster_names[combinations[2]],".png"), width = 12, height = 16)
          
        }
      }
      #setting up some metrics for SMD plot
      cluster_count <- length(unique(eicu_df$cluster))
      temp <- t(combn(seq(from = 1, to = cluster_count, by = 1), 2))
      combination_list <- split(temp, seq(nrow(temp)))
      varnames <- colnames(eicu_df)
      test_alg <- "EICU_vs_MIMIC"
      test_distance <- ""
      test_name <- "Kmdist_gower"
      mimic_df <- mimic_df[,-2, drop = F]
      
      SMD_plot_two_types(df1 = eicu_df, df2 = mimic_df, type1 = "eICU", type2 = "MIMIC", varnames, category_name = "EICU_vs_MIMIC", list_combination = combination_list, save_directory = save_dir, percent_from_each_end = 20)
      
      
    }

  }
  
}
