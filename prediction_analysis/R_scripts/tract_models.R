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

library(glmnet)
library(caret)
library(dplyr)
library(tidyverse)
library(pcLasso)
library(rsample)
library(R.utils)
library(sparsegl)

source('functions.R') 

# +
default_args <- list(seed_val = 1, 
                     family_sep=TRUE, 
                     model='lasso', 
                     theta=1, 
                     path_spec="notebook_results", 
                     debugging = TRUE,
                     bootstrap = TRUE, 
                     path_prefix = '/pscratch/sd/m/mphagen',

                    sel_metrics = c('dki_fa', 'dki_md', 'dki_mk', 'dki_awf')) 
args <- R.utils::commandArgs(defaults=default_args, asValues=TRUE)

seed <- args$seed_val
family_sep <- args$family_sep
theta <- args$theta
model <- args$model
path_spec <- args$path_spec 
debugging <- args$debugging
bootstrap <- args$bootstrap
path_prefix <- args$path_prefix

sel_metrics <- args$sel_metrics
# -

results_path <- file.path(path_prefix, 'hcp-results', 'tract_models',
              paste0(Sys.Date(), '_', system("git rev-parse --short HEAD", intern=TRUE)), 
                    Reduce(file.path, path_spec))
dir.create(results_path, recursive=TRUE, showWarnings=FALSE)

save(args, file=file.path(results_path, paste0(seed, '_arguments.Rdata' )))

measures_filename <- "Data/unrestricted_mphagen_1_27_2022_20_50_7.csv"
profiles_filename <- "Data/bundle_profiles_all.csv"
restricted_filename <- "Data/RESTRICTED_arokem_1_31_2022_23_26_45.csv"

# +
res_df <- read.csv(restricted_filename) %>% column_to_rownames('Subject') 
measure_df <- read.csv(measures_filename)%>% column_to_rownames('Subject') 
temp_pheno_df <- merge(res_df, measure_df, by="row.names")  %>% column_to_rownames("Row.names")


profiles_data <- read.csv(profiles_filename) %>% 
    select(any_of(c('tractID', 'nodeID', sel_metrics, 'subject')))

profiles_wide <- profiles_data %>% 
    pivot_wider(names_from=c("tractID", "nodeID"), 
                values_from=all_of(sel_metrics))  %>% 
    column_to_rownames('subject')

profiles_wide <- merge(profiles_wide, temp_pheno_df %>% 
                 select(all_of(c('Family_ID', behavioral_list))),
                 by='row.names') 
# -

set.seed(seed)

if (bootstrap == TRUE) { 
    test_boot <- profiles_wide %>%  
              nest(data = everything(), .by=Family_ID) %>%
              dplyr::slice_sample(n = 434, replace=TRUE) 

    counter <- sapply(unique(temp_pheno_df[['Family_ID']]), function(x) 0)

    for (ii in 1:length(test_boot$data)) { 
        fam_id <- unique(test_boot$data[[ii]][['Family_ID']])
        count <- counter[[fam_id]]
        test_boot$data[[ii]][['subject']] <- paste(test_boot$data[[ii]][['Row.names']], count, sep='_')
        counter[[fam_id]] = counter[[fam_id]] + 1 } 

    profiles_boot <- unnest(test_boot, cols=c(data), names_repair = "minimal") 

    hist(unlist(counter), breaks=6) 

    profiles_wide <- profiles_boot %>% 
                      select(!which(  as.list(profiles_boot) %>% duplicated())) %>% 
                      column_to_rownames('subject') %>% select(-'Row.names')}  else { 
    
        profiles_wide <- profiles_wide  %>% column_to_rownames('Row.names') 
                      } 

vfold <- group_vfold_cv(profiles_wide,
                        Family_ID, v=5)
train_subjects <- vfold$splits[[1]] %>% analysis %>% rownames()
test_subjects <- vfold$splits[[1]] %>% assessment %>% rownames()

write.csv(train_subjects, 
          paste0(results_path, 
                 '/', seed, '_train_subjects.csv'), 
         row.names=FALSE)

write.csv(test_subjects, 
          paste0(results_path, 
                 '/', seed, '_test_subjects.csv'), 
         row.names=FALSE)

write.csv(profiles_wide$Family_ID, 
          paste0(results_path, 
                 '/', seed, '_family_IDs.csv'),
         row.names=FALSE)

train_data <- profiles_wide %>% 
    filter(rownames(profiles_wide) %in% train_subjects)                      
test_data <- profiles_wide %>% 
     filter(rownames(profiles_wide) %in% test_subjects) 

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
        if (length(sel_metrics) > 1) { 
            num_groups = 96 } 
        else { 
            num_groups = 24 } 
        for (ii in 1:num_groups) { 
            groups <- append(groups, rep(ii, 100)) 
        } 

        sgl_model <- cv.sparsegl(y = as.matrix(y_train), 
                              x = as.matrix(x_train),
                              group=groups ) 

        r2 <- R2(predict(sgl_model, 
                newx=as.matrix(x_test), s='lambda.min'), 
                as.matrix(y_test))

        return(list(model=sgl_model, rsq=r2) ) } 

    } 


# +
test_rsq <- list()
models <- list() 
for (behavior in behavioral_list) { 
    
    print(behavior)

    print("prepping y") 
    y_train <- train_data %>% select(all_of(behavior)) %>% na.omit()

    y_test <- test_data %>% select(all_of(behavior)) %>% na.omit()

    write.csv(y_test,  
          paste0(results_path, 
                 '/', seed, 'y_test.csv'),
         row.names=FALSE) 
    
     write.csv(y_train,  
          paste0(results_path, 
                 '/', seed, 'y_train.csv'),
         row.names=FALSE)   
    
    print("prepping x") 
    x_train <- train_data %>% 
        select(-behavioral_list) %>% 
        select(-contains("Family")) %>% 
        filter(rownames(train_data) %in% rownames(y_train)) %>% 
        mutate(across(contains("_"), scale))
    
        
     x_test <- test_data %>% 
        select(-behavioral_list) %>% 
        select(-contains("Family")) %>% 
        filter(rownames(test_data) %in% rownames(y_test)) %>% 
        mutate(across(contains("_"), scale))  

    check_scaling(x_train) 
    check_scaling(x_test) 
    stopifnot(identical(rownames(y_train), rownames(x_train))) 
    
    
    out <- run_model(function_name=model, x_train=x_train, y_train=y_train,
                                        x_test=x_test, y_test=y_test, theta=theta) 
    
    test_rsq[behavior] <-  out$rsq 

    models[[behavior]] <- out$model
    
    #make this an if statement
         
        
} 

# +
write.csv(data.frame(test_rsq),
          paste0(results_path, '/test_r2_vals_', seed, '.csv'))

save(models, file=paste0(results_path, '/model_data_', seed, '.Rdata'))

# -

if (debugging == TRUE) { 
    save.image(paste0(results_path, '/image', seed, '.Rdata')) } 


