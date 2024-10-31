## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(survML)
library(survival)
library(dplyr)
library(ggplot2)
set.seed(72924)

## ----km, fig.width=8, fig.height=6--------------------------------------------
data(cancer)
km_fit <- survfit(Surv(rfstime, status) ~ 1, data = gbsg)
plot(km_fit, xlab = "Time (days)", ylab = "Recurrence-free survival probability")

## ----data_setup---------------------------------------------------------------
### variables of interest
# rfstime - recurrence-free survival
# status - censoring indicator
# hormon - hormonal therapy treatment indicator
# age - in years
# meno - 1 = premenopause, 2 = post
# size - tumor size in mm
# grade - factor 1,2,3
# nodes - number of positive nodes
# pgr - progesterone receptor in fmol
# er - estrogen receptor in fmol

# create dummy variables and clean data
gbsg$tumgrad2 <- ifelse(gbsg$grade == 2, 1, 0)
gbsg$tumgrad3 <- ifelse(gbsg$grade == 3, 1, 0)
gbsg <- gbsg %>% na.omit() %>% select(-c(pid, grade))

time <- gbsg$rfstime
event <- gbsg$status
X <- gbsg %>% select(-c(rfstime, status)) # remove outcome 

# find column indices of features/feature groups
X_names <- names(X)
age_index <- paste0(which(X_names == "age"))
meno_index <- paste0(which(X_names == "meno"))
size_index <- paste0(which(X_names == "size"))
nodes_index <- paste0(which(X_names == "nodes"))
pgr_index <- paste0(which(X_names == "pgr"))
er_index <- paste0(which(X_names == "er"))
hormon_index <- paste0(which(X_names == "hormon"))
grade_index <- paste0(which(X_names %in% c("tumgrad2", "tumgrad3")), collapse = ",")
tum_index <- paste0(which(X_names %in% c("size", "nodes", "pgr", "er", "tumgrad2", "tumgrad3")),
                    collapse = ",")
person_index <- paste0(which(X_names %in% c("age", "meno", "hormon")), collapse = ",")

feature_group_names <- c("age", "meno.", "size", "nodes",
                         "prog.", "estro.", "hormone", 
                         "grade")
feature_groups <- c(age_index, meno_index, size_index, nodes_index,
                    pgr_index, er_index, hormon_index, grade_index)


## ----fig.width=8, fig.height=8------------------------------------------------
# landmark times for AUC
landmark_times <- c(1000, 2000)

output <- vim(type = "AUC",
              time = time,
              event = event,
              X = X,
              landmark_times = landmark_times,
              large_feature_vector = 1:ncol(X),
              small_feature_vector = (1:ncol(X))[-as.numeric(age_index)],
              conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                        V = 2,
                                                        bin_size = 0.5),
              large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                    V = 2),
              cf_fold_num = 2,
              sample_split = FALSE,
              scale_est = TRUE)
output$result$indx <- rep(age_index, nrow(output$result))
output$result$indx_name <- rep("age", nrow(output$result))
output$result

## -----------------------------------------------------------------------------
# save the objects that we will reuse
saved_conditional_surv_preds <- output$conditional_surv_preds
saved_large_oracle_preds <- output$large_oracle_preds
saved_folds <- output$folds
saved_approx_times <- output$approx_times

pooled_output <- output$result # save the results for age

# iterate over other feature groups
for (i in 2:length(feature_group_names)){
  indx_char <- feature_groups[i]
  indx_name <- feature_group_names[i]
  indx <- as.numeric(strsplit(indx_char, split = ",")[[1]])
  
  output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = 1:ncol(X),
                small_feature_vector = (1:ncol(X))[-indx],
                conditional_surv_preds = saved_conditional_surv_preds,
                large_oracle_preds = saved_large_oracle_preds,
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                sample_split = TRUE,
                scale_est = TRUE)
  
  output$result$indx <- rep(indx_char, nrow(output$result))
  output$result$indx_name <- rep(indx_name, nrow(output$result))
  pooled_output <- rbind(pooled_output, output$result)
}

plot_results <- function(results, plot_title){
  # plot results
  p_auc <- results %>%
    mutate(landmark_time = factor(landmark_time,
                                  levels = c(1000, 2000),
                                  labels = c("1000 days", "2000 days"))) %>%
    arrange(landmark_time, est) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = est, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
        geom_point() +
        theme_bw() +
        xlab("Estimated importance") +
        ylab("Feature group") +
        xlim(c(0,0.5)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_name,
        ) + 
        facet_wrap(~landmark_time, dir = "v", strip.position = "right", scales = "free_y", ncol = 1) + 
        ggtitle(plot_title)+
        theme(strip.background = element_blank(),
              strip.placement = "outside")
    }
  return(p_auc)
}

p_auc <- plot_results(pooled_output, "AUC variable importance relative to full feature vector")
p_auc

## ----fig.width=8, fig.height=5------------------------------------------------
# consider joint importance of all tumor-level and person-level features
feature_group_names2 <- c("tumor", "person")
feature_groups2 <- c(tum_index, person_index)
# repeat the analysis for feature groups
for (i in 1:length(feature_group_names2)){
  indx_char <- feature_groups2[i]
  indx_name <- feature_group_names2[i]
  indx <- as.numeric(strsplit(indx_char, split = ",")[[1]])
  
  output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = 1:ncol(X),
                small_feature_vector = (1:ncol(X))[-indx],
                conditional_surv_preds = saved_conditional_surv_preds,
                large_oracle_preds = saved_large_oracle_preds,
                small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                sample_split = TRUE,
                scale_est = TRUE)
  
  output$result$indx <- rep(indx_char, nrow(output$result))
  output$result$indx_name <- rep(indx_name, nrow(output$result))
  if (!(i == 1)){
    pooled_output <- rbind(pooled_output, output$result)
  } else{
    pooled_output <- output$result
  }
}

p_auc <- plot_results(pooled_output, "AUC variable importance relative to full feature vector (groups)")
p_auc

## ----fig.width=8, fig.height=7------------------------------------------------
# For importance relative to baseline features, the 'small' model uses only person-level (baseline) features
# The 'large' model uses baseline + feature of interest
size_index <- paste0(c(size_index, person_index), collapse = ",")
nodes_index <- paste0(c(nodes_index, person_index), collapse = ",")
pgr_index <- paste0(c(pgr_index, person_index), collapse = ",")
er_index <- paste0(c(er_index, person_index), collapse = ",")
grade_index <- paste0(c(grade_index, person_index), collapse = ",")

feature_group_names <- c("size", "nodes", "prog.", "estro.", "grade")
feature_groups <- c(size_index, nodes_index,
                    pgr_index, er_index, grade_index)

person_index_numeric <- as.numeric(strsplit(person_index, split = ",")[[1]])

for (i in 1:length(feature_group_names)){
  indx_char <- feature_groups[i]
  indx_name <- feature_group_names[i]
  indx <- as.numeric(strsplit(indx_char, split = ",")[[1]])
  
  if (i == 1){
    output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = indx,
                small_feature_vector =  person_index_numeric,
                conditional_surv_preds = saved_conditional_surv_preds,
                large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                sample_split = TRUE,
                scale_est = TRUE)
    saved_small_oracle_preds <- output$small_oracle_preds
  } else{
    output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = indx,
                small_feature_vector =  person_index_numeric,
                conditional_surv_preds = saved_conditional_surv_preds,
                small_oracle_preds = saved_small_oracle_preds,
                large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                sample_split = TRUE,
                scale_est = TRUE)
  }
  
  output$result$indx <- rep(indx_char, nrow(output$result))
  output$result$indx_name <- rep(indx_name, nrow(output$result))
  if (!(i == 1)){
    pooled_output <- rbind(pooled_output, output$result)
  } else{
    pooled_output <- output$result
  }
}

p_auc <- plot_results(pooled_output, "AUC variable importance relative to person-level features")
p_auc


## ----fig.width=8, fig.height=7------------------------------------------------
size_index <- paste0(c(size_index, person_index), collapse = ",")
nodes_index <- paste0(c(nodes_index, person_index), collapse = ",")
pgr_index <- paste0(c(pgr_index, person_index), collapse = ",")
er_index <- paste0(c(er_index, person_index), collapse = ",")
grade_index <- paste0(c(grade_index, person_index), collapse = ",")

feature_group_names <- c("size", "nodes", "prog.", "estro.", "grade")
feature_groups <- c(size_index, nodes_index,
                    pgr_index, er_index, grade_index)

for (i in 1:length(feature_group_names)){
  indx_char <- feature_groups[i]
  indx_name <- feature_group_names[i]
  indx <- as.numeric(strsplit(indx_char, split = ",")[[1]])
  
  if (i == 1){
     output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = (1:ncol(X))[-person_index_numeric],
                small_feature_vector =  (1:ncol(X))[-indx],
                conditional_surv_preds = saved_conditional_surv_preds,
                large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                sample_split = TRUE,
                scale_est = TRUE)
     saved_large_oracle_preds <- output$large_oracle_preds
  } else{
     output <- vim(type = "AUC",
                time = time,
                event = event,
                X = X,
                landmark_times = landmark_times,
                approx_times = saved_approx_times,
                large_feature_vector = (1:ncol(X))[-person_index_numeric],
                small_feature_vector =  (1:ncol(X))[-indx],
                conditional_surv_preds = saved_conditional_surv_preds,
                large_oracle_preds = saved_large_oracle_preds,
                small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm"),
                                                      V = 2),
                cf_folds = saved_folds$cf_folds,
                ss_folds = saved_folds$ss_folds,
                sample_split = TRUE,
                scale_est = TRUE)
  }
  
  output$result$indx <- rep(indx_char, nrow(output$result))
  output$result$indx_name <- rep(indx_name, nrow(output$result))
  if (!(i == 1)){
    pooled_output <- rbind(pooled_output, output$result)
  } else{
    pooled_output <- output$result
  }
}

p_auc <- plot_results(pooled_output, "Adjusted AUC variable importance relative to all tumor-level features")
p_auc

