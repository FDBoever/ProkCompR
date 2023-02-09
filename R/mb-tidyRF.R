library(tidyverse)
library(tidymodels)
library(ranger)

# Do the random forest with tidymodels ranger instead?
#   think deeply first


#CLASS IMBALANCE

pan.df <- cazy.pan[rownames(metadata),] %>% tibble::rownames_to_column(var="genome") %>% right_join(metadata[,c("BinId","group")], by=c('genome'='BinId')) %>% select(-genome)
pan.df <- pan.df %>% filter(group %in% c('cl1','cl4','cl7','cl14'))

pan_rf.folds <- pan.df %>%
  mutate(group = factor(group)) %>%
  select(-genome) %>%
  bootstraps(strata = group)


rf_spec <- rand_forest(trees = 1000) %>%
  set_mode("classification") %>%
  set_engine("ranger")


#Build workflow
pan_rf.wf <- workflow() %>%
  add_formula(group ~ .) %>%
  add_model(rf_spec)

pan_rf.wf

#Fit
doParallel::registerDoParallel()
pan_rf.rs <- fit_resamples(
  pan_rf.wf,
  resamples = pan_rf.folds,
  control = control_resamples(save_pred = TRUE)
)

pan_rf.rs



#---- done

collect_metrics(pan_rf.rs)


pan_rf.rs %>%
  collect_predictions() %>%
  group_by(id) %>%
  ppv(group, .pred_class)


pan_rf.rs %>%
  collect_predictions() %>%
  group_by(id) %>%
  roc_curve(group, .pred_cl1:.pred_cl14) %>%
  ggplot(aes(1 - specificity, sensitivity, color = id)) +
  geom_abline(lty = 2, color = "gray80", size = 1.5) +
  geom_path(show.legend = FALSE, alpha = 0.6, size = 1.2) +
  facet_wrap(~.level, ncol = 5) +
  coord_equal()


pan_rf.rs %>% collect_predictions() %>%
  conf_mat(group, .pred_class) %>%
  autoplot(type = "heatmap")



glm_spec <- logistic_reg() %>%
  set_engine("glm")
glm_spec

pan_glm.wf <- workflow() %>%
  add_formula(group ~ .) %>%
  add_model(glm_spec)

members_metrics <- metric_set(accuracy, sensitivity, specificity)

doParallel::registerDoParallel()
pan_glm.rs <- fit_resamples(
  pan_glm.wf,
  metrics = members_metrics,
  resamples = pan_rf.folds#,
  #control = control_resamples(save_pred = TRUE)
)

pan_glm.rs





#-----------
#training and testing data
set.seed(1234)
members_split <- initial_split(pan.df, strata = group)
members_train <- training(members_split)
members_test <- testing(members_split)

#To evaluate model performance, we resample
#cross validation folds, each is a fold of data,
##use resamples, to tune, or compare models
set.seed(1234)
members_folds <- vfold_cv(members_train, strata = group)
members_folds


#data preprossesing
#Class imbalance is handled with smote from the themis package
#smote, uses nearest neigbors (upsample)
members_rec <- recipe(group ~ ., data = members_train)# %>%
  #step_medianimpute(age) %>%
  #step_other(peak_id, citizenship) %>%
  #step_dummy(all_nominal(), -died) %>%
  #themis::step_smote(group)


glm_spec <- logistic_reg() %>%
  set_engine("glm")
glm_spec


rf_spec <- rand_forest(trees = 1000) %>%
  set_mode("classification") %>%
  set_engine("ranger")
rf_spec


members_wf <- workflow() %>%
  add_recipe(members_rec)
members_wf




members_metrics <- metric_set(roc_auc, accuracy, sensitivity, specificity)
doParallel::registerDoParallel()
glm_rs <- members_wf %>%
  add_model(glm_spec) %>%
  fit_resamples(
    resamples = members_folds,
    metrics = members_metrics,
    control = control_resamples(save_pred = TRUE)
  )

glm_rs

rf_rs <- members_wf %>%
  add_model(rf_spec) %>%
  fit_resamples(
    resamples = members_folds,
    metrics = members_metrics,
    control = control_resamples(save_pred = TRUE)
  )

rf_rs


members_final <- members_wf %>%
  add_model(glm_spec) %>%
  last_fit(members_split)




#------
# XGBoost
library(XGBoost)


members_split <- initial_split(pan.df, strata = group)
members_train <- training(members_split)
members_test <- testing(members_split)

vb_split <- initial_split(pan.df, strata = group)
vb_train <- training(vb_split)
vb_test <- testing(vb_split)

xgb_spec <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune(),                         ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_spec

xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), vb_train),
  learn_rate(),
  size = 30
)

xgb_grid


#-- run model

doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)

xgb_res




xgb_wf <- workflow() %>%
  add_formula(group ~ .) %>%
  add_model(xgb_spec)

xgb_wf

set.seed(123)
vb_folds <- vfold_cv(vb_train, strata = group)

vb_folds


doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)

xgb_res



collect_metrics(xgb_res)


xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

show_best(xgb_res, "roc_auc")

best_auc <- select_best(xgb_res, "roc_auc")
best_auc

final_xgb <- finalize_workflow(
  xgb_wf,
  best_auc
)

final_xgb


library(vip)

final_xgb %>%
  fit(data = vb_train) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")



final_res <- last_fit(final_xgb, vb_split)
collect_metrics(final_res)

final_res %>%
  collect_predictions() %>%
  roc_curve(group, .pred_cl1:.pred_cl4) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1.5, color = "midnightblue") +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    size = 1.2
  )
