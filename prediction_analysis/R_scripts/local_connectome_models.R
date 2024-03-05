# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: r-env
#     language: R
#     name: ir
# ---

library(rhdf5) 
library(tidyverse)
library(rsample)

getwd() 

source('functions.R') 

measures_filename <- "Data/unrestricted_mphagen_1_27_2022_20_50_7.csv"
lc_filename <- "Data/final_data_R1.hdf5"
restricted_filename <- "Data/RESTRICTED_arokem_1_31_2022_23_26_45.csv"

res_df <- read.csv(restricted_filename)
measure_df <- read.csv(measures_filename)
pheno_df <- merge(res_df, measure_df, by='Subject') 

lc_df <- h5read(lc_filename, name='loc_conn_features') 

subjects <-  h5read(lc_filename, name='subjects') 

lc_data <- as.data.frame(t(lc_df), row.names=subjects) 

family_df <- hcp_family_data(restricted_filename, dwi_df=lc_data) 

vfold <- group_vfold_cv(family_df, Family_ID, v=5)
train_subjects <- vfold$splits[[1]] %>% analysis %>% pull(Subject)
test_subjects <- vfold$splits[[1]] %>% assessment %>% pull(Subject)

# +
train_lc <- lc_data %>% 
                    rownames_to_column() %>%
                    filter(rowname %in% train_subjects)  %>%
                     mutate(across(contains("V"), scale)) %>% 
                    column_to_rownames()
                     
test_lc <- lc_data %>% 
                rownames_to_column() %>%
                filter(rowname %in% test_subjects) %>% 
                     mutate(across(contains("V"), scale)) %>% 
                    column_to_rownames()

# -

behavioral_list <- c('PicSeq_Unadj', 'Flanker_Unadj', 
                     'ReadEng_Unadj', 
                     'VSPLOT_TC', 'CogFluidComp_Unadj',
                     'CogTotalComp_Unadj', 'CogCrystalComp_Unadj',
                     'Endurance_Unadj', 'NEOFAC_A',
                     'ASR_Extn_T', 'ASR_Totp_T', 'Age_in_Yrs')

results_path <- file.path('Results', 'lc_models', system("git rev-parse HEAD", intern=TRUE ))
dir.create(results_path, recursive=TRUE)

# +
test_rsq_list <- list() 
train_rsq_list <- list()
for (behavior in behavioral_list) { 
    
    print(behavior)
#add in plain lasso/ridge modeling
    y_train <- pheno_df %>% filter(Subject %in% train_subjects) %>% arrange('Subject')%>% 
                column_to_rownames('Subject')  %>% 
                select(behavior) %>% na.omit() 

    y_test <- pheno_df %>% filter(Subject %in% test_subjects) %>% arrange('Subject')%>% 
                column_to_rownames('Subject')  %>% 
                select(behavior) %>% na.omit() 
    
    x_train <- train_lc  %>% rownames_to_column() %>% 
                filter(rowname %in% rownames(y_train))   %>% 
                arrange(rowname) %>% column_to_rownames()

    x_test <- test_lc  %>% rownames_to_column() %>% 
                filter(rowname %in% rownames(y_test))   %>% 
                arrange(rowname) %>% column_to_rownames()

    pca_results <- prcomp(x_train, scale.=FALSE, center=FALSE) 
    
    pca_x_train <- pca_results$x
    
    pca_x_test <-  as.matrix(x_test) %*% pca_results$rotation
    
    stopifnot(identical(rownames(y_train), rownames(x_train))) 
    
    
    pcr_model <- cv.glmnet(x=as.matrix(pca_x_train),
                            y=as.matrix(y_train), 
                            alpha=1)
    
    glm_model <- cv.glmnet(x=as.matrix(x_train),
                            y=as.matrix(y_train), 
                            alpha=1)
    
    
    # pcL_model <- cv.pcLasso(x=as.matrix(x_train), 
    #                        y=as.matrix(y_train), 
    #                        ratio=1) 
        
    test_rsq_list[['glm']][behavior] <- R2(predict(glm_model, 
                newx=as.matrix(x_test),
                s="lambda.min"), 
                as.matrix(y_test)) #use select instead 
    
    test_rsq_list[['pcr']][behavior] <- R2(predict(pcr_model, 
                newx=as.matrix(pca_x_test),
                s="lambda.min"), 
                as.matrix(y_test))
    
    # test_rsq_list[['pcl']][behavior] <- R2(predict(pcL_model, 
    #             xnew=as.matrix(x_test),
    #             s="lambda.min"), 
    #             as.matrix(y_test))

#     print(R2(predict(pca_model, 
#                 xnew=x_test,
#                 s="lambda.min"), 
#                 y_test)) 
    }
