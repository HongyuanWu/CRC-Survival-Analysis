### MUTATUION DRIVEN SURVIVAL ANALYSIS 
#Use mutation's group obtained precedently to analyze survival of each group
#Through this way we can achieve a good accuracy for each group and 
# highlight (emphasize) the effect of a particular group of mutations and of 
# every mutations in the group.

#1) Create temp_dataset both for MSS and MSI.
#   Add a binary column for each mutation where 1 = mutation occurred
#2) Cox Model for each group
# 2.1) For each group make a different model that use mutations as covariates,
#      then test the model with Cox Snell Residuals and Scaled Schoenfeld 
#      Residuals.
# 2.2) For each group make a different model that use mutations as covariates
#      and stratify for some covariates (every covariate that looks not
#      proportional after residuals tests)
#3) Stratified Cox Model
#   Make a unique stratified Cox Model where each mutation of each group
#   is a covariate.
#   Stratify the dataset in according to the membership of the patients to the
#   groups.
#   Then test the model with residuals' tests.

#1) create temp_dataset
MSSsurvival.data_temp = add_mutations_as_column(MSSsurvival.data, 
                                                MSSlift$annotations,
                                                MSSlift$genotypes)
MSIsurvival.data_temp = add_mutations_as_column(MSIsurvival.data,
                                                MSIlift$annotations,
                                                MSIlift$genotypes)
#Cox Model for each group
# - make the model
# - obtain Baseline Cumulative Hazard Function and Survival Function
# - test the model with Cox Snell Residuals (to check overall fit) and Scaled 
#   Schoned Residuals (to check proportionally assumption for each covariate)
#makeCoxModelsAndTests(dataset, patients_of_each_group, mutations_of_each_group) 
# returns a list of list --> MS*groups_cox$(models, survival, 
# baseline_cumulative_hazard, cox_snell_residuals_surv_fun, 
# scaled_schoenfeld_residuals)
#MSS
MSSgroups_cox = makeCoxModelsAndTests(MSSsurvival.data_temp, #dataset
                                      MSSfullgroups, #groups to subset dataset
                                      MSSgroups_grouping1) #groups with mutation
#MSI
MSIgroups_cox = makeCoxModelsAndTests(MSIsurvival.data_temp, #dataset
                                      MSIfullgroups, #groups to subset dataset
                                      MSIgroups_grouping1) #groups with mutation
#For each valid model: 
# - plot functions and tests' results
#MSS
for(i in 1:length(MSSgroups_cox$models)) {
  if(is.list(MSSgroups_cox$models[[i]])){ #is list if is a cox model
    plotCoxSurvFun(MSSgroups_cox$survival[[i]], i = i, 
                   type = "MSS", col = "darkblue")
    plotCoxBslnCmHazFun(MSSgroups_cox$baseline_cumulative_hazard[[i]], i = i,
                        type = "MSS", col = "darkorange1")
    plotCoxSnellResSurvFun(MSSgroups_cox$cox_snell_residuals_surv_fun[[i]], i = i,
                           type = "MSS", col = "forestgreen")
    plotScaledSchoenfeldResiduals(MSSgroups_cox$scaled_schoenfeld_residuals[[i]],
                                  i = i, type = "MSS", col = "darkred")
  }
}
#MSI
for(j in 1:length(MSIgroups_cox$models)) {
  if(is.list(MSIgroups_cox$models[[j]])){ #is list if is a cox model
    plotCoxSurvFun(MSIgroups_cox$survival[[j]], i = j, 
                   type = "MSI", col = "darkblue")
    plotCoxBslnCmHazFun(MSIgroups_cox$baseline_cumulative_hazard[[j]], i = j,
                        type = "MSI", col = "darkorange1")
    plotCoxSnellResSurvFun(MSIgroups_cox$cox_snell_residuals_surv_fun[[j]], i = j,
                           type = "MSI", col = "forestgreen")
    plotScaledSchoenfeldResiduals(MSIgroups_cox$scaled_schoenfeld_residuals[[j]],
                                  i = j, type = "MSI", col = "darkred")
  }
}
#Results not good.

############################## OTHER MODELS ####################################
#Differents Models only for MSI_status = 'MSS'--> NOT GENERAL TESTS

#1) Unique Stratified Cox Model with all mutations as covariate
#add group's information in the dataset
MSSsurvival.data_temp = addGroupsAsColumn(MSSsurvival.data_temp, MSSfullgroups)
#create the formula and ignore the group which have a null cox model
# in the precedent analysis
#flatten of the mutations group list for MSS and remove duplicates, but do it
# only for group with a valid model in the precedent cox model.
mutations = c()
for(i in 1:length(MSSgroups_cox$models)) {
  if(is.list(MSSgroups_cox$models[[i]])) {
    mutations = append(mutations, MSSgroups_grouping1[[i]])
  }
}
mutations = unique(mutations)
#add all mutations as covariate in the formula
formula = paste0("SurvObj ~ ", mutations[1])
for(j in 2:length(mutations)){ 
  formula = paste0(formula, " + ", mutations[j])
}
#add the stratification (by group) in the formula ignoring group with null
# model (null model in precedent cox model)
stratification = " + strata("
c = 0 #check if is first group (add only "group'i'" and not "+ group'i'")
for(i in 1:length(MSSgroups_cox$models)) {
  if(is.list(MSSgroups_cox$models[[i]])) { #if realtive model isn't null
    c = c + 1
    if(c == 1) {
      stratification = paste0(stratification, "group", i)
    } else {
      stratification = paste0(stratification,", group", i)
    }
  }
}
stratification = paste0(stratification, ')')
formula = as.formula(paste0(formula, stratification))
#subset the dataset (NOT GENERAL, ONLY FOR THIS DATASET)
dataset_notnullgroups = subset(MSSsurvival.data_temp, 
                               MSSsurvival.data_temp$group2 == 1 |
                                 MSSsurvival.data_temp$group3 == 1 | 
                                 MSSsurvival.data_temp$group4 == 1 |
                                 MSSsurvival.data_temp$group5 == 1 |
                                 MSSsurvival.data_temp$group6 == 1 |
                                 MSSsurvival.data_temp$group7 == 1)
#make the cox model
MSScox_strata = coxph(formula = formula, data = dataset_notnullgroups)
#Result not good : too much NA covariates --> all survival function is line at 1
#The problem is that a patient can be in more then one groups, so no single
# group could have a good survival function (it will be good only for strata
# that consider more groups)
#To see all the strata and his survival function : summary(survfit(prova))

#2) Unique Stratified Cox Model without covariates (strata = groups)
formula2 = as.formula(paste0("SurvObj ~ 1 ", stratification))
MSScox_strata_2 = coxph(formula = formula2, data = dataset_notnullgroups)
#summary(survfit(prova2)) to see the not null survival function
# plot only it
plot(survfit(MSScox_strata_2), lwd = 1.5,
     xlab = "Time(Days)", ylab = "Survival(%)",
     main = "Survival Function of Cox Model\nstratified by MSS groups",
     col = c(NA, NA, NA, NA, NA,"darkmagenta", "blue", NA, NA, NA, NA, NA, NA,
             NA, NA, NA, NA, NA, "green", NA, NA, "black", NA, "red", NA, 
             "darkorange", NA, NA, NA, NA, NA, NA, NA))
grid(NULL, NULL, lty = 5)
legend("bottomleft", col = c("darkmagenta", "blue", "green"), 
       pch = 1, cex = 0.65, lwd = c(1.5, 1.5, 1.5),
       legend = c("group4 AND group7", "group4 AND group6", 
                  "group4 AND group5 AND group6"))
legend("bottomright", col = c("black", "red", "darkorange"),
       pch = 1, cex = 0.65, lwd = c(1.5, 1.5, 1.5), 
       legend = c("group2 AND group6", "group2 AND group4", 
                  "group2 AND group4 AND group6 AND group7"))
#find a test to do
#test

#observing precedents models ==> KEEP ONLY GROUP 2 AND GROUP 4 
#3) Stratified Cox Model only for group2 and group4
dataset = MSSsurvival.data_temp
dataset$id = sort(sample(1:length(dataset[,1])))
dataset_g2 = dataset[dataset$id %in% MSSfullgroups[[2]],]
dataset_g4 = dataset[dataset$id %in% MSSfullgroups[[4]],]
dataset_g2 = subset(dataset_g2, !(rownames(dataset_g2) %in% rownames(dataset_g4)))
dataset_g2g4 = rbind(dataset_g2, dataset_g4)
dataset_g2g4$SurvObj = Surv(dataset_g2g4$time, dataset_g2g4$status)
prova3 = coxph(formula = SurvObj ~ 1 + strata(group2, group4), 
               data = dataset_g2g4)
#summary(survfit(prova3)) -> 3 strata: only group4, only group2, both group
plot(survfit(prova3, conf.int = T), 
     col = c("red", "blue", "purple"), lwd = 2,
     main = "Survival Function of Cox Model\nstratified by MSS group2 and group4",
     xlab = "Time(Days)", ylab = "Survival(%)")
grid(NULL, NULL, lty = 5)
legend("bottomleft",
       cex = 0.75, pch = 1, col = c("red", "blue", "purple"), lwd = c(1.5,1.5),
       legend = c("group 4", "group 2", "both groups"))
#adding a formula
formula_g2g4 = paste0("SurvObj ~ ", MSSgroups_grouping1[[2]][1])
for(i in 2:length(MSSgroups_grouping1[[2]])) {
  formula_g2g4 = paste0(formula_g2g4, " + ", MSSgroups_grouping1[[2]][i])
}
for(i in 1:length(MSSgroups_grouping1[[4]])) {
  formula_g2g4 = paste0(formula_g2g4, " + ", MSSgroups_grouping1[[4]][i])
}
formula_g2g4 = as.formula(paste0(formula_g2g4, " + strata(group2, group4)"))
cox_g2g4 = coxph(formula = formula_g2g4, data = dataset_g2g4)
plot(survfit(cox_g2g4), lwd = c(2,3,2) , col = c("red", "blue", "purple"),
     xlab = "Time(Days)", ylab = "Survival(%)",
     main = paste0("Survival Function for Cox Model stratified\n",
                   "by MSS group2 and group4 with mutations as covariates"))
grid(NULL, NULL, lty = 5)
legend("bottomright",
       cex = 0.75, pch = 1, col = c("red", "blue", "purple"), lwd = c(1.5,1.5),
       legend = c("group 4", "group 2", "both groups"))
#test !!
#Not good cause of covariates

#Try considering group2 and group4 separately
#cox models
model_g2 = MSSgroups_cox$models[[2]]
model_g4 = MSSgroups_cox$models[[4]]
surv_g2 = MSSgroups_cox$survival[[2]]
surv_g4 = MSSgroups_cox$survival[[4]]
plot(surv_g2, col = "red", xlab = "Time(Days)", ylab = "Survival(%)", 
     main = paste0("Survival Function of Cox Model for MSS groups 2 and 4\n",
                   "with mutations as covariates"))
grid(NULL, NULL, lty = 5)
lines(surv_g4, col = "blue")
legend("bottomleft", lwd = 1.5, pch = 1, col = c("red", "blue"),
       legend = c ("group 2", "group 4"))
#log rank test
group = c()
for(i in 1:length(dataset_g2g4$group2)) {
  if(dataset_g2g4$group2[i] == 1) {
    if(dataset_g2g4$group4[i] == 1) {
      group = append(group, "both")
    } else {
      group = append(group, 2)
    }
  } else {
    group = append(group, 4)
  }
}
dataset_g2g4only = dataset_g2g4
dataset_g2g4only$group = group
log_rank_test_g2g4 = survdiff(formula = SurvObj ~ group, 
                              data = dataset_g2g4)


############################# LEVEL GROUPING ###################################
#group by level (only for MSS cause MSI is too little)
#obtain dataset for each level
level1_dataset = MSSsurvival.data[MSSlevelGroups[[1]],]
level2_dataset = MSSsurvival.data[MSSlevelGroups[[2]],]
level3_dataset = MSSsurvival.data[MSSlevelGroups[[3]],]
#Kaplan-meier for each group
km_level1 = survfit(formula = SurvObj ~ 1, data = level1_dataset)
summary(km_level1)
km_level2 = survfit(formula = SurvObj ~ 1, data = level2_dataset)
summary(km_level2)
km_level3 = survfit(formula = SurvObj ~ 1, data = level3_dataset)
summary(km_level3)
#life tables for each km curve
#lifetab(tis: time points, ninit = # patients,
#        nlost = # censored data, nevent = # events --> deaths )
lifetab_km_level1 = lifetab(tis = append(km_level1$time, NA), 
                            ninit = km_level1$n, nlost = km_level1$n.censor, 
                            nevent = km_level1$n.event)
lifetab_km_level2 = lifetab(tis = append(km_level2$time, NA), 
                            ninit = km_level2$n, nlost = km_level2$n.censor, 
                            nevent = km_level2$n.event)
lifetab_km_level3 = lifetab(tis = append(km_level3$time, NA), 
                            ninit = km_level3$n, nlost = km_level3$n.censor, 
                            nevent = km_level3$n.event)
#write output in a CSV file after checking the existence of the directory
#check directory existence or create it
ifelse(!dir.exists(file.path("CsvLifeTables")), 
       dir.create(file.path("CsvLifeTables")), FALSE)
#write output
write.csv(lifetab_km_level1, file = "CsvLifeTables/LifeTableLevel1.csv")
write.csv(lifetab_km_level2, file = "CsvLifeTables/LifeTableLevel2.csv")
write.csv(lifetab_km_level3, file = "CsvLifeTables/LifeTableLevel3.csv")
#Add a column for the levels in the full dataset
level_column = c()
for(i in 1:length(MSSsurvival.data$time)) {
  if(i %in% MSSlevelGroups[[1]]) {
    level_column = append(level_column, 1)
  } else if(i %in% MSSlevelGroups[[2]]) {
    level_column = append(level_column, 2)
  } else if(i %in% MSSlevelGroups[[3]]) {
    level_column = append(level_column, 3)
  } else { #patient isn't in any level
    level_column = append(level_column, 0)
  }
}
MSSsurvival_levels_dataset = MSSsurvival.data
MSSsurvival_levels_dataset$level = level_column
#Log-Rank Test to compare each group
log_rank_levels = survdiff(formula = SurvObj ~ level,
                    data = MSSsurvival_levels_dataset)
#find p-value
p_value_log_rank_levels = pchisq(log_rank_levels$chisq, 
                                 length(log_rank_levels$n)-1, 
                                 lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for all levels (no levels, level1, level2, level3")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_levels > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level1 and level2
MSSsurvival_levels_1_2_dataset = subset(MSSsurvival_levels_dataset,
                                        level == 1 | level ==2)
log_rank_levels_1_2 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_levels_1_2_dataset)
#find p-value
p_value_log_rank_levels_1_2 = pchisq(log_rank_levels_1_2$chisq, 
                                 length(log_rank_levels_1_2$n)-1, 
                                 lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 1 and 2")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_levels_1_2 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level2 and level3
MSSsurvival_levels_2_3_dataset = subset(MSSsurvival_levels_dataset,
                                        level == 2 | level ==3)
log_rank_levels_2_3 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_levels_2_3_dataset)
#find p-value
p_value_log_rank_levels_2_3 = pchisq(log_rank_levels_2_3$chisq, 
                                     length(log_rank_levels_2_3$n)-1, 
                                     lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 2 and 3")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_levels_2_3 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level1 and level3
MSSsurvival_levels_1_3_dataset = subset(MSSsurvival_levels_dataset,
                                        level == 1 | level == 3)
log_rank_levels_1_3 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_levels_1_3_dataset)
#find p-value
p_value_log_rank_levels_1_3 = pchisq(log_rank_levels_1_3$chisq, 
                                     length(log_rank_levels_1_3$n)-1, 
                                     lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 1 and 3")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_levels_1_3 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#plot all km curves
plot(km_level1, conf.int = FALSE, lwd = 1.5, col = "red",
     xlab = "Time(Days)", ylab = "Survival(%)", 
     main = paste0("Kaplan-Meier estimation of Survival Function\n",
                   "for different levels on MSS graph"))
grid(NULL, NULL, lty = 5)
lines(km_level2, conf.int = FALSE, lwd = 1.5, col = "darkgreen")
lines(km_level3, conf.int = FALSE, lwd = 1.5, col = "blue")
legend("bottomleft", col = c("red","darkgreen","blue"), lwd = c(1.5, 1.5, 1.5),
       legend = c("level1 (depth on graph = 1)", "level2 (depth on graph = 2)",
                  "level3 (depth on graph = 3)"), pch = 1, cex = 0.75)
#Cox model using column level as covariate
cox_levels = coxph(formula = SurvObj ~ level, data = MSSsurvival_levels_dataset)
summary(cox_levels)
#survival
plot(survfit(cox_levels), col = "forestgreen", lwd = 1.5, 
     xlab = "Time(Days)", ylab = "Survival(%)",
     main = "Survival Function of Cox Model with level as covariate")
grid(NULL, NULL, lty = 5)
#Check overall fit with CoxSnellResiduals
l_martingale_residuals = residuals(cox_levels, type = 'martingale')
l_cox_snell_residuals = MSSsurvival_levels_dataset$status - l_martingale_residuals
l_res_surv = survfit(Surv(l_cox_snell_residuals, 
                          MSSsurvival_levels_dataset$status) ~ 1)
plot(l_res_surv$time, -log(l_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model with level as covariate",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
grid(NULL, NULL, lty = 5)
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
#check covariate fitting
#Check proportional hazard assumption : Scaled Schoenfeld Residuals
ssr_cox_levels = cox.zph(cox_levels)
plot(ssr_cox_levels, df = 2, main = "Scaled Schoenfeld Residual test \n of Cox Model with level as covariate",
     col = "darkblue", lwd = 1.5)

#Cox model using every level as a differente covariate
#create the columns
level1_column = c()
level2_column = c()
level3_column = c()
for(i in 1:length(MSSsurvival_levels_dataset$time)) {
  if(i %in% MSSlevelGroups[[1]]) {
    level1_column = append(level1_column, 1)
    level2_column = append(level2_column, 0)
    level3_column = append(level3_column, 0)
  } else if(i %in% MSSlevelGroups[[2]]) {
    level1_column = append(level1_column, 0)
    level2_column = append(level2_column, 1)
    level3_column = append(level3_column, 0)
  } else if(i %in% MSSlevelGroups[[3]]) {
    level1_column = append(level1_column, 0)
    level2_column = append(level2_column, 0)
    level3_column = append(level3_column, 1)
  } else { #patient isn't in any level
    level1_column = append(level1_column, 0)
    level2_column = append(level2_column, 0)
    level3_column = append(level3_column, 0)
  }
}
MSSsurvival_levels_dataset$level1 = level1_column
MSSsurvival_levels_dataset$level2 = level2_column
MSSsurvival_levels_dataset$level3 = level3_column
cox_levels_2 = coxph(formula = SurvObj ~ level1 + level2 + level3,
                     data = MSSsurvival_levels_dataset)
summary(cox_levels_2)
#Check overall fit with CoxSnellResiduals
l2_martingale_residuals = residuals(cox_levels_2, type = 'martingale')
l2_cox_snell_residuals = 
  MSSsurvival_levels_dataset$status - l2_martingale_residuals
l2_res_surv = survfit(Surv(l2_cox_snell_residuals,
                           MSSsurvival_levels_dataset$status) ~ 1)
#plot pl2_res_surv cumulative hazard with reference line x=y to check it
plot(l2_res_surv$time, -log(l2_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model with levels as covariates",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
grid(NULL, NULL, lty = 5)
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
#check covariate fitting
#Check proportional hazard assumption : Scaled Schoenfeld Residuals
ssr_cox_levels_2 = cox.zph(cox_levels_2)
plot(ssr_cox_levels_2, df = 2, main = "Scaled Schoenfeld Residual test \n of Cox Model with level as covariate",
     col = "darkblue", lwd = 1.5)

#Cox Model stratified by level
cox_levels_strata = coxph(formula = SurvObj ~ strata(level),
                          data = MSSsurvival_levels_dataset)
summary(cox_levels_strata)
#survival
surv_levels_strata = survfit(cox_levels_strata)
plot(surv_levels_strata, lwd = 1.5, xlab = "Time(Days)", ylab = "Survival(%)",
     main = "Survival function of Cox Model Stratified by level",
     col = c(NA,"red", "green", "blue"))
grid(NULL, NULL, lty = 5)
legend("bottomleft", col = c("red", "green", "blue"), cex = 0.75, pch = 1,
       legend = c("strata: level1 (depth 1 on the graph)", 
                  "strata: level2 (depth 2 on the graph)", 
                  "strata: level3 (depth 3 on the graph)"))
#Check overall fit with CoxSnellResiduals
l_strata_martingale_residuals = residuals(cox_levels_strata, type = 'martingale')
l_strata_cox_snell_residuals = 
  MSSsurvival_levels_dataset$status - l_strata_martingale_residuals
l_strata_res_surv = survfit(Surv(l_strata_cox_snell_residuals, 
                                 MSSsurvival_levels_dataset$status) ~ 1)
plot(l_strata_res_surv$time, -log(l_strata_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model Stratified by level",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
grid(NULL, NULL, lty = 5)


########################### PATH-LEVELS GROUPING ###############################
#group by level and paths (only for MSS cause MSI is too little)
#obtain dataset for each level
pathlevel1_dataset = MSSsurvival.data[MSSpathLevelGroups[[1]],]
pathlevel2_dataset = MSSsurvival.data[MSSpathLevelGroups[[2]],]
pathlevel3_dataset = MSSsurvival.data[MSSpathLevelGroups[[3]],]
#Kaplan-meier for each group
km_pathlevel1 = survfit(formula = SurvObj ~ 1, data = pathlevel1_dataset)
summary(km_pathlevel1)
km_pathlevel2 = survfit(formula = SurvObj ~ 1, data = pathlevel2_dataset)
summary(km_pathlevel2)
km_pathlevel3 = survfit(formula = SurvObj ~ 1, data = pathlevel3_dataset)
summary(km_pathlevel3)
#life tables for each km curve
#lifetab(tis: time points, ninit = # patients,
#        nlost = # censored data, nevent = # events --> deaths )
lifetab_km_pathlevel1 = lifetab(tis = append(km_pathlevel1$time, NA), 
                            ninit = km_pathlevel1$n,
                            nlost = km_pathlevel1$n.censor, 
                            nevent = km_pathlevel1$n.event)
lifetab_km_pathlevel2 = lifetab(tis = append(km_pathlevel2$time, NA), 
                            ninit = km_pathlevel2$n, 
                            nlost = km_pathlevel2$n.censor, 
                            nevent = km_pathlevel2$n.event)
lifetab_km_pathlevel3 = lifetab(tis = append(km_pathlevel3$time, NA), 
                            ninit = km_pathlevel3$n, 
                            nlost = km_pathlevel3$n.censor, 
                            nevent = km_pathlevel3$n.event)
#write output in a CSV file after checking the existence of the directory
#check directory existence or create it
ifelse(!dir.exists(file.path("CsvLifeTables")), 
       dir.create(file.path("CsvLifeTables")), FALSE)
#write output
write.csv(lifetab_km_pathlevel1, file = "CsvLifeTables/LifeTablePathLevel1.csv")
write.csv(lifetab_km_pathlevel2, file = "CsvLifeTables/LifeTablePathLevel2.csv")
write.csv(lifetab_km_pathlevel3, file = "CsvLifeTables/LifeTablePathLevel3.csv")
#Add a column for the levels in the full dataset
pathlevel_column = c()
for(i in 1:length(MSSsurvival.data$time)) {
  if(i %in% MSSpathLevelGroups[[1]]) {
    pathlevel_column = append(pathlevel_column, 1)
  } else if(i %in% MSSpathLevelGroups[[2]]) {
    pathlevel_column = append(pathlevel_column, 2)
  } else if(i %in% MSSpathLevelGroups[[3]]) {
    pathlevel_column = append(pathlevel_column, 3)
  } else { #patient isn't in any level
    pathlevel_column = append(pathlevel_column, 0)
  }
}
MSSsurvival_pathlevels_dataset = MSSsurvival.data
MSSsurvival_pathlevels_dataset$level = pathlevel_column
#Log-Rank Test to compare each group
log_rank_pathlevels = survdiff(formula = SurvObj ~ level,
                           data = MSSsurvival_pathlevels_dataset)
#find p-value
p_value_log_rank_pathlevels = pchisq(log_rank_pathlevels$chisq, 
                                 length(log_rank_pathlevels$n)-1, 
                                 lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for all levels considering also path (no levels, level1, level2, level3")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_pathlevels > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level1 and level2
MSSsurvival_pathlevels_1_2_dataset = subset(MSSsurvival_pathlevels_dataset,
                                        level == 1 | level ==2)
log_rank_pathlevels_1_2 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_pathlevels_1_2_dataset)
#find p-value
p_value_log_rank_pathlevels_1_2 = pchisq(log_rank_pathlevels_1_2$chisq, 
                                     length(log_rank_pathlevels_1_2$n)-1, 
                                     lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 1 and 2 considering also path")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_pathlevels_1_2 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level2 and level3
MSSsurvival_pathlevels_2_3_dataset = subset(MSSsurvival_pathlevels_dataset,
                                        level == 2 | level ==3)
log_rank_pathlevels_2_3 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_pathlevels_2_3_dataset)
#find p-value
p_value_log_rank_pathlevels_2_3 = pchisq(log_rank_pathlevels_2_3$chisq, 
                                     length(log_rank_pathlevels_2_3$n)-1, 
                                     lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 2 and 3 considering also path")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_pathlevels_2_3 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#compare only level1 and level3
MSSsurvival_pathlevels_1_3_dataset = subset(MSSsurvival_pathlevels_dataset,
                                        level == 1 | level == 3)
log_rank_pathlevels_1_3 = survdiff(formula = SurvObj ~ level,
                               data = MSSsurvival_pathlevels_1_3_dataset)
#find p-value
p_value_log_rank_pathlevels_1_3 = pchisq(log_rank_pathlevels_1_3$chisq, 
                                     length(log_rank_pathlevels_1_3$n)-1, 
                                     lower.tail = FALSE)
#use alpha = 0.1
print("Log Rank Test for levels 1 and 3 considering also path")
print("H0 = \"The groups follow the same distribution\"")
if(p_value_log_rank_levels_1_3 > 0.1) {
  print("Can't reject H0")
} else {
  print("H0 rejected")
}
#plot all km curves
plot(km_pathlevel1, conf.int = FALSE, lwd = 1.5, col = "red",
     xlab = "Time(Days)", ylab = "Survival(%)", 
     main = paste0("Kaplan-Meier estimation of Survival Function\n",
                   "for different levels following any path on MSS graph "))
grid(NULL, NULL, lty = 5)
lines(km_pathlevel2, conf.int = FALSE, lwd = 1.5, col = "darkgreen")
lines(km_pathlevel3, conf.int = FALSE, lwd = 1.5, col = "blue")
legend("bottomleft", col = c("red","darkgreen","blue"), lwd = c(1.5, 1.5, 1.5),
       legend = c("level1 (depth on graph = 1)", "level2 (depth on graph = 2)",
                  "level3 (depth on graph = 3)"), pch = 1, cex = 0.75)
#Cox model using column level as covariate
cox_pathlevels = coxph(formula = SurvObj ~ level,
                   data = MSSsurvival_pathlevels_dataset)
summary(cox_pathlevels)
#survival
plot(survfit(cox_pathlevels), col = "forestgreen", lwd = 1.5, 
     xlab = "Time(Days)", ylab = "Survival(%)",
     main = "Survival Function of Cox Model with level on a path as covariate")
grid(NULL, NULL, lty = 5)
#Check overall fit with CoxSnellResiduals
pl_martingale_residuals = residuals(cox_pathlevels, type = 'martingale')
pl_cox_snell_residuals = MSSsurvival_pathlevels_dataset$status - pl_martingale_residuals
pl_res_surv = survfit(Surv(pl_cox_snell_residuals, 
                                  MSSsurvival_pathlevels_dataset$status) ~ 1)
plot(pl_res_surv$time, -log(pl_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model with level on any path as covariate",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
grid(NULL, NULL, lty = 5)
#check covariates fitting with scaled schoenfeld residuals
ssr_cox_pathlevels = cox.zph(cox_pathlevels)
plot(ssr_cox_pathlevels, df = 2, main = "Scaled Schoenfeld Residual test of Cox Model\nwith level on any path as covariate",
     col = "darkblue", lwd = 1.5)

#Cox model using every level as a differente covariate
#create the columns
pathlevel1_column = c()
pathlevel2_column = c()
pathlevel3_column = c()
for(i in 1:length(MSSsurvival_pathlevels_dataset$time)) {
  if(i %in% MSSpathLevelGroups[[1]]) {
    pathlevel1_column = append(pathlevel1_column, 1)
    pathlevel2_column = append(pathlevel2_column, 0)
    pathlevel3_column = append(pathlevel3_column, 0)
  } else if(i %in% MSSpathLevelGroups[[2]]) {
    pathlevel1_column = append(pathlevel1_column, 0)
    pathlevel2_column = append(pathlevel2_column, 1)
    pathlevel3_column = append(pathlevel3_column, 0)
  } else if(i %in% MSSpathLevelGroups[[3]]) {
    pathlevel1_column = append(pathlevel1_column, 0)
    pathlevel2_column = append(pathlevel2_column, 0)
    pathlevel3_column = append(pathlevel3_column, 1)
  } else { #patient isn't in any level
    pathlevel1_column = append(pathlevel1_column, 0)
    pathlevel2_column = append(pathlevel2_column, 0)
    pathlevel3_column = append(pathlevel3_column, 0)
  }
}
MSSsurvival_pathlevels_dataset$level1 = pathlevel1_column
MSSsurvival_pathlevels_dataset$level2 = pathlevel2_column
MSSsurvival_pathlevels_dataset$level3 = pathlevel3_column
cox_pathlevels_2 = coxph(formula = SurvObj ~ level1 + level2 + level3,
                     data = MSSsurvival_pathlevels_dataset)
summary(cox_pathlevels_2)
#Check overall fit with CoxSnellResiduals
pl2_martingale_residuals = residuals(cox_pathlevels_2, type = 'martingale')
pl2_cox_snell_residuals = 
  MSSsurvival_pathlevels_dataset$status - pl2_martingale_residuals
pl2_res_surv = survfit(Surv(pl2_cox_snell_residuals,
                            MSSsurvival_pathlevels_dataset$status) ~ 1)
#plot pl2_res_surv cumulative hazard with reference line x=y to check it
plot(pl2_res_surv$time, -log(pl2_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model with levels on a path as covariates",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
grid(NULL, NULL, lty = 5)
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
#check covariate fitting
#Check proportional hazard assumption : Scaled Schoenfeld Residuals
ssr_cox_pathlevels_2 = cox.zph(cox_pathlevels_2)
plot(ssr_cox_pathlevels_2, df = 2, main = "Scaled Schoenfeld Residual test \n of Cox Model with level on a path as covariate",
     col = "darkblue", lwd = 1.5)

#Cox Model stratified by level
cox_pathlevels_strata = coxph(formula = SurvObj ~ strata(level),
                          data = MSSsurvival_pathlevels_dataset)
summary(cox_pathlevels_strata)
#survival
surv_pathlevels_strata = survfit(cox_pathlevels_strata)
plot(surv_pathlevels_strata, lwd = 1.5, xlab = "Time(Days)", ylab = "Survival(%)",
     main = "Survival function of Cox Model Stratified by level\n on any path on MSS graph",
     col = c(NA,"red", "green", "blue"))
grid(NULL, NULL, lty = 5)
legend("bottomleft", col = c("red", "green", "blue"), cex = 0.75, pch = 1,
       legend = c("strata: level1 (depth 1 on the graph)", 
                  "strata: level2 (depth 2 on the graph)", 
                  "strata: level3 (depth 3 on the graph)"))
#Check overall fit with CoxSnellResiduals
pl_strata_martingale_residuals = residuals(cox_pathlevels_strata, type = 'martingale')
pl_strata_cox_snell_residuals = 
  MSSsurvival_pathlevels_dataset$status - pl_strata_martingale_residuals
pl_strata_res_surv = survfit(Surv(pl_strata_cox_snell_residuals, 
                           MSSsurvival_pathlevels_dataset$status) ~ 1)
plot(pl_strata_res_surv$time, -log(pl_strata_res_surv$surv), #cumulative hazard of residuals
     type = 'S', lwd = 1.5, col = 'blue',
     main = "Overall fit by Cox Snell Residuals of\n Cox Model Stratified by level on any path",
     xlab = "residuals", ylab = "Cumulative Hazard of residuals")
abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
grid(NULL, NULL, lty = 5)