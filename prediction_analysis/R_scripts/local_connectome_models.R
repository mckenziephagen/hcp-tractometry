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
#     display_name: R [conda env:r-env-3]

#     language: R
#     name: conda-env-r-env-3-r
# ---

library(rhdf5) 
library(tidyverse)
library(rsample)
library(glmnet) 
library(caret) 
library(pcLasso) 
library(dtplyr)
library(data.table)

# +
default_args <- list(seed_val = 1, family_sep=FALSE, 
                    bootstrap=FALSE)
args <- R.utils::commandArgs(defaults=default_args, asValues=TRUE)

seed <- args$seed_val
bootstrap <- args$bootstrap


# -

source('functions.R') 

behavioral_list <- c('PicSeq_Unadj', 'Flanker_Unadj', 'CardSort_Unadj',
                     'ReadEng_Unadj', 'DDisc_AUC_200', 'IWRD_TOT',
                     'SCPT_SEN', 'VSPLOT_TC', 'CogFluidComp_Unadj',
                     'CogTotalComp_Unadj', 'CogCrystalComp_Unadj',
                     'Endurance_Unadj', 'NEOFAC_A', 'NEOFAC_O',
                     'NEOFAC_C', 'NEOFAC_N', 'NEOFAC_E',
                     'ASR_Intn_T', 'ASR_Extn_T', 'ASR_Totp_T')

measures_filename <- "Data/unrestricted_mphagen_1_27_2022_20_50_7.csv"
lc_filename <- "Data/final_data_R1.hdf5"
restricted_filename <- "Data/RESTRICTED_arokem_1_31_2022_23_26_45.csv"

res_df <- read.csv(restricted_filename) %>% column_to_rownames('Subject') 
measure_df <- read.csv(measures_filename)%>% column_to_rownames('Subject')
temp_pheno_df <- merge(res_df, measure_df, by="row.names")  %>% column_to_rownames("Row.names")

lc_df <- h5read(lc_filename, name='loc_conn_features') 
subjects <-  h5read(lc_filename, name='subjects') 
lc_data <- as.data.frame(t(lc_df), row.names=subjects) %>% 
                 merge( temp_pheno_df %>% 
                 select(all_of(c('Family_ID', behavioral_list))),
                 by='row.names') 

set.seed(seed)

lc_data %>% select(all_of(c('Family_ID', behavioral_list)))

#test_boot is a bad name for this
if (bootstrap == TRUE) { 
    lc_boot <- lc_data %>% nest(data = everything(), .by=Family_ID) %>%
              dplyr::slice_sample(n = 434, replace=TRUE) %>% 
              unnest( cols=c(data), names_repair = "minimal") 

    lc_data <- lc_boot %>% select(!which(  as.list(lc_boot) %>% duplicated())) %>% rename(subject=Row.names) } 
else { 
    

vfold <- group_vfold_cv(data.frame(lc_data %>% select(c('Family_ID', 'subject'))), Family_ID, v=5)
train_subjects <- vfold$splits[[1]] %>% analysis %>% pull(subject)
test_subjects <- vfold$splits[[1]] %>% assessment %>% pull(subject)

# +
#this can be updated to data.table too 

train_lc <- lc_data %>% 
            filter(subject %in% train_subjects)  
                     
test_lc <- lc_data %>% 
            filter(subject %in% test_subjects) 
# -

results_path <- file.path('Results', 'lc_models',
                          paste0(Sys.Date(), '_', system("git rev-parse --short HEAD", intern=TRUE)), 
                          Reduce(file.path, 'baseline'))
dir.create(results_path, recursive=TRUE)

run_model <- function(function_name='lasso', x_train, y_train, x_test, y_test, alpha=1, theta=1, ...) { 

    if (function_name == 'pcr') { 
        
        pca_results <- prcomp(x_train, scale.=FALSE, center=FALSE) 
        pca_x_train <- pca_results$x
        pca_x_test <-  as.matrix(x_test) %*% pca_results$rotation 
        #sub in predict and check 
        
        pcr_model <- cv.glmnet(x=as.matrix(pca_x_train),
                            y=as.matrix(y_train), 
                            alpha=alpha)


        rsq <- R2(predict(pcr_model, 
                newx=as.matrix(pca_x_test),
                s="lambda.min"), 
                as.matrix(y_test))

        return(list(model=pcr_model, rsq=rsq) )

        } 

    else if (function_name == 'lasso') {
         print("running lasso") 
         lasso_model <- cv.glmnet(x=as.matrix(x_train),
                            y=as.matrix(y_train), 
                            alpha=1)


        r2 <- R2(predict(lasso_model, 
                newx=as.matrix(x_test),
                s="lambda.min"), 
                as.matrix(y_test))

        return(list(model=lasso_model, rsq=r2) )
        }


    else if (function_name == 'pcLasso') { 

        pcL_model <- cv.pcLasso(x=as.matrix(x_train), 
                           y=as.matrix(y_train), 
                           theta=theta) 

        r2 <- R2(predict(pcL_model, 
                xnew=as.matrix(x_test),
                s="lambda.min"), 
                as.matrix(y_test))
        
        return(list(model=pcL_model, rsq=r2) ) 
        }
    else if (function_name == 'sgl') { 

        groups=c()
        for (ii in 1:24) { 
        groups <- append(groups, rep(ii, 100)) } 

        sgl_model <- cv.sparsegl(y = as.matrix(y_train), 
                              x = as.matrix(x_train),
                              group=groups ) 

        r2 <- R2(predict(sgl_model, 
                newx=as.matrix(x_test), s='lambda.min'), 
                as.matrix(y_test))

        return(list(model=sgl_model, rsq=r2) ) 

    } } 


model='lasso'

# +
test_rsq <- list()
models <- list() 
for (behavior in behavioral_list) { 
    
    print(behavior)

    print("prepping y") 
    y_train <- train_lc %>% select(all_of(behavior)) %>% na.omit()

    y_test <- test_lc %>% select(all_of(behavior)) %>% na.omit()
    
    cat(format(Sys.time(), format = "%H:%M:%S"),
        file=paste0(results_path, "/time_file.txt"),sep="\n")
        print("prepping x") 
    
    x_train <- train_lc %>% 
        select(-all_of(behavioral_list)) %>% 
        select(-contains("Family")) %>% 
        filter(rownames(train_lc) %in% rownames(y_train)) %>% 
        as.data.table(keep.rownames=TRUE)
    
    cols <- names(x_train)[names(x_train) %like% "V"]
    x_train[, (cols) := lapply(.SD, function(x) (x - mean(x))/sd(x)),
             .SDcols = cols]

    x_train <- x_train %>% as.data.frame() %>% column_to_rownames('rn')
    
    cat(format(Sys.time(), format = "%H:%M:%S"),
        file=paste0(results_path, "/time_file.txt"),sep="\n", append=TRUE)

     x_test <- test_lc %>% 
        select(-all_of(behavioral_list)) %>% 
        select(-contains("Family")) %>% 
        filter(rownames(test_lc) %in% rownames(y_test)) %>% 
        as.data.table()

    cols <- names(x_train)[names(x_train) %like% "V"]
    x_test[, (cols) := lapply(.SD, function(x) (x - mean(x))/sd(x)),
             .SDcols = cols]                        
    
    cat(format(Sys.time(), format = "%H:%M:%S"),
        file=paste0(results_path, "/time_file.txt"),sep="\n", append=TRUE)
                      
    stopifnot(identical(rownames(y_train), rownames(x_train))) 
    
    
    out <- run_model(function_name=model, x_train=x_train, y_train=y_train,
                                        x_test=x_test, y_test=y_test, theta=theta) 
    
    test_rsq[behavior] <-  out$rsq 

    models[[behavior]] <- out$model
          
} 

# +
#try this next: https://stackoverflow.com/questions/29583665/data-table-alternative-for-dplyr-mutate
# -

save(models, file=paste0(results_path, '/models_', seed, '.csv'))

# +
write.csv(data.frame(test_rsq),
          paste0(results_path, '/test_r2_vals_', seed, '.csv'))



save.image(paste0(results_path, '/model_data_', seed, '.Rdata'))
# -



