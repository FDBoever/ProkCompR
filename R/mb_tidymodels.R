#Load necessary packages
list.of.packages <- c("tidyverse", "tidymodels", "ranger", "kknn","janitor",'pryr','ggExtra')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cloud.r-project.org")
invisible(capture.output(lapply(list.of.packages, library, character.only = TRUE, warn.conflicts = FALSE, quietly = T)))


#===========

metdad <- genome.tbl %>%
  filter(genome%in% sorted.cliques$genome) %>%
  left_join(sorted.cliques,by='genome') %>%
  filter(grepl('Marinobacter|Tamil',genome)) %>%
  left_join(out.colors, by='genome') %>%
  filter(!is.na(!!as.name(selcol))) %>%
  filter(quality_class != 'low') %>%
  data.frame()
metdad$group <- paste0('cl',metdad[,colcol])

rownames(metdad) <- metdad$genome
metdad[,selcol] <- as.factor(metdad[,selcol])
metdad[,colcol] <- as.factor(metdad[,colcol])


colors <- metdad %>% select(paste0('col.',colcol),colcol, group) %>% unique() %>% data.frame()
clique.colors <- as.character(colors[,1])
names(clique.colors) <- as.character(colors[,'group'])

#===========


#set random seed for reproducability
set.seed(123)

#split the dataset into a training and testing set
metadata_split <- rsample::initial_split(metdad, prop = 0.9,strate='group') #create a split object where 9/10 proportion of data is used for training
metadata_test <- rsample::testing(metadata_split) #isolate training set into one object
metadata_train <- rsample::training(metadata_split) #isolate testing set into one object
tibble::glimpse(metadata_test) #peek at plant testing set

metadata_recipe <- head(metadata_train) %>%
  recipes::recipe(GC ~ Coding_density + Genome_size + group)

#if the model formula is extremely long, it may be necessary to define the model formula using update_role functions
metadata_recipe2 <- head(metadata_train) %>%
  recipes::recipe() %>%
  recipes::update_role(GC, new_role = "outcome") %>% #
  recipes::update_role(Coding_density, new_role = "predictor")

#"selector functions can be used to select a group of variables at one time
metadata_recipe3 <- head(metadata_train) %>%
  recipes::recipe() %>%
  recipes::update_role(matches("ep."), new_role = "outcome") %>%
  recipes::update_role(contains("N"), new_role = "predictor")

#you can also yse steps to transform variables with step functions (https://recipes.tidymodels.org/reference/index.html#section-step-functions-individual-transformations)
metadata_recipe4 <- head(metadata_train) %>%
  recipes::recipe(GC ~ Coding_density + Genome_size + group) %>%
  recipes::step_center(all_predictors(), -all_outcomes()) %>%
  recipes::step_scale(all_predictors(), -all_outcomes())


#execute the recipe on test data to validate
metadata_recipe <- recipes::prep(metadata_recipe)
#extract data where the recipe was applied
metadata_head <- recipes::juice(metadata_recipe)

#apply recipe to data of interest
metadata_train <- recipes::bake(metadata_recipe, metadata_train)
metadata_test <- recipes::bake(metadata_recipe, metadata_test)
tibble::glimpse(metadata_test)

#Define model
linear_lm <- linear_reg() %>% #define a linear regression model
  set_engine("lm") %>% #decide what package you want to use for the linear regression
  set_mode("regression")  #set regression or classification depending on model specifications

#alternative model arrangement for comparison
linear_glmnet <- linear_reg() %>%
  set_engine("stan") %>%
  set_mode("regression")

#note engines do not have to be R packages, tidymodels supports interfacing with outside programs such as keras or tensorflow. One of the benefits to using R.

#create a workflow
metadata_workflow <- workflows::workflow() %>% #define workflow
  workflows::add_recipe(metadata_recipe) %>% #add recipe specification
  workflows::add_model(linear_lm) #add model specification

#fit workflow to data
metadata_fit <- parsnip::fit(metadata_workflow, metadata_train)

#predict test data from fitting model
metadata_pred <- predict(metadata_fit, metadata_test)
#bind prediction results to test data for easy comparison
metadata_test <- dplyr::bind_cols(metadata_test, metadata_pred)
tibble::glimpse(metadata_test)

#calculate model performance metrics (default)
metadata_metrics <- yardstick::metrics(metadata_test, truth = GC, estimate = .pred)
metadata_metrics

#calculate a single model performance metric
yardstick::rmse(data = metadata_test, truth = GC, estimate = .pred)

#calculate a custom set of model performance metrics
metadata_multi <- yardstick::metric_set(rmse, rsq)
metadata_multi(data = metadata_test, truth = GC, estimate = .pred)

#Split training data for cross validation scheme into 5 folds with 5 repeats
metadata_folds <- rsample::vfold_cv(training(metadata_split), v = 5, repeats = 5)
#define a control object to change settings in the cross validation process
metadata_control <- tune::control_grid(
  verbose = FALSE,  #print out progress while fitting
  allow_par = TRUE, #allow parallel processing
  extract = function(x){extract_model(x)},  #extract the individual fitting model object for each split
  save_pred = TRUE #save model predictions
)

#fit the workflow to the folds object, with control settings
metadata_cv <- tune::fit_resamples(metadata_workflow, metadata_folds, control = metadata_control)
#collect performance metrics from the folds
metadata_metrics <- tune::collect_metrics(metadata_cv)
#collect predictions from the folds
metadata_predictions <- tune::collect_predictions(metadata_cv)

metadata_predictions %>% ggplot(aes(GC, .pred))+geom_point()+geom_smooth(method='lm')

#view
metadata_cv$.extracts[[1]][[1]]








#========================================================================================#
#Definee a new recipe for random forest model
rf_recipe <- head(training(metadata_split)) %>%
  recipes::recipe(group ~ genome + Completeness + Contamination + SS + Coding_density) %>%
  recipes::update_role(genome, new_role = "id") %>%
  recipes::step_string2factor(group)

#Define a new model
rf_tune <- parsnip::rand_forest(trees = 100, mtry = tune(), min_n = tune()) %>% #mtry and min_n parameters are set to "tune()"
  parsnip::set_engine("ranger") %>% #use ranger package
  parsnip::set_mode("classification") #set to classification

#define new workflow with recipe and model
tune_wf <- workflows::workflow() %>%
  workflows::add_recipe(rf_recipe) %>%
  workflows::add_model(rf_tune)

#tune the models
tune_res <- tune::tune_grid(tune_wf, resamples = metadata_folds, grid = 5)

#collect metrics for tuning
tune_metrics <- tune::collect_metrics(tune_res)

#reshape data from wide to long format for plotting
p1data <- tune_metrics %>%
  tidyr::pivot_longer(min_n:mtry, values_to = "value",names_to = "parameter")

#plot data
p1 <- p1data %>%
  ggplot(aes(x = value, y = mean)) +
  geom_point() +
  facet_wrap(parameter~.metric, scales = "free")
p1

#select the best parameters from tuning routine
best_params <- tune::select_best(tune_res, "roc_auc")
#final model using parameters selected from tune results
final_rf <- tune::finalize_model(rf_tune, best_params)

#add final model into workflow
final_wf <- workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)

#fit the finalized model to the testing data
final_res <- final_wf %>% tune::last_fit(metadata_split)
collect_metrics(final_res)



############################################

#set seed to random integer
set.seed(Sys.time())
#create list to store folds
custom_folds <- list()
#create 25 splits of the training data
for(i in 1:25){
  custom_folds[[i]] <- initial_split(training(metadata_split), prop = 0.8)
}
#convert list folds to tibble
custom_folds <- tibble(custom_folds)

#compare object size
pryr::object_size(metadata_folds)
pryr::object_size(custom_folds)

###########################################

#Create Sample data
big_mat <- as.data.frame(matrix(rexp(20000, rate=2), ncol=200))
big_mat$var <- as.factor(sample(c(TRUE, FALSE), nrow(big_mat), replace = TRUE))

#Define PCA Recipe
recipe <- head(big_mat) %>%
  recipes::recipe(var ~ .) %>%
  recipes::step_center(all_predictors()) %>%
  recipes::step_scale(all_predictors()) %>%
  step_pca(all_predictors(), num_comp = 2)

#Define KNN model
model <- nearest_neighbor(neighbors = 5, weight_func = "rectangular") %>%
  set_engine("kknn") %>%
  set_mode("classification")

#Create workflow
workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(model)

#fit model
fit <- fit(workflow, big_mat)
#model predictions
preds <- predict(fit, big_mat)
#bind model predictions to original data
res <- cbind(big_mat, preds) %>%
  select(var, .pred_class)

#calculate performance metrics
res_metrics <- metrics(res, truth = var, estimate = .pred_class)
