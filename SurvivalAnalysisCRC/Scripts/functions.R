### FUNCTIONS ###

#function to add descendents to each roots --> return groups
add_descendents = function (glist, gedges, groots) {
  for(i in 1:length(groots)) {#1 NOT 9
    glist = add_descendents_ric(glist, gedges, groots[i], i)
  }
  return(glist)
}
add_descendents_ric = function (glist, gedges, root, index_of_root) {
  arr = which(gedges[,1] == root) #index of descendents from 'root'
  if(length(arr) != 0) {
    #for each descendent recursively find his descendents
    for(j in 1:length(arr)) { 
      subroot = gedges[arr[j],2] #get id from index
      glist[index_of_root] = mapply(c, glist[index_of_root], 
                                    subroot, SIMPLIFY = F)
      glist = add_descendents_ric(glist, gedges, subroot, index_of_root)
    }
  } 
  return(glist)
}

#function to remove samples from ***clinical.data that aren't also in ***lift
#(*** = MSS or MSI)
merge_samples = function(dataset, lift_patients_list){
  data_patients_list = rownames(dataset)
  arr = c()
  for(i in 1:length(data_patients_list)) {
    #if the element is in both dataset --> add it to arr
    if(length(which(lift_patients_list == data_patients_list[i])) != 0){
      arr = append(arr, i)
    }
  }
  #subsetting considering only arr
  dataset = dataset[arr,]
  return(dataset)
}

#function to complete the path of each group
#complete_path = function(groups, edges, roots) {
#   
#}

#function to remove duplicates, OR and XOR from groups
# and substitute mutation's id with mutation's name
simplify_groups = function(group_list, nodes_dataset) {
  for(i in 1:length(group_list)) {
    for(j in 1:length(group_list[[i]])) {
      index = which(nodes_dataset[,1] == group_list[[i]][j])
      if(length(index) > 0) {
        name = nodes_dataset[index, 2]  #mutation's name 
        if(name != "OR" & name != "XOR" & name != "UPOR" & name != "UPXOR") {
          type = paste0('_',nodes_dataset[index, 3])  #type(del., mutat., ...)
          group_list[[i]][j] = paste0(name, type)
        }
        else{
          group_list[[i]][j] = name 
        }
      }
    }
  }
  #remove all XOR, OR, UPXOR and UPOR
  group_list = lapply(group_list, function(x) x[x != "XOR" & x != "OR" & x != "UPOR" & x != "UPXOR"])
  #remove all duplicates
  group_list = lapply(group_list, function(x) x[!duplicated(x)])
  return(group_list)
}

#function to assing group to each patient of a dataset
#mutations = mutations'table [mutations'ID <-> type <-> mutation's name]
#patient_mutation_table = binary matrix [raw: patients, col: mutations'ID]
#for each group, for each patient
# check if is 1 [mutation's ID index - patient's index] in the matrix for
# each mutation in the group.
# OBS: matrix keep ID, so ID --> NAME --> check with group's mutations
assign_patients_to_group = function(dataset, grouplist, 
                                    mutations, patients_mutation_table){
  grouping1 = list(c(),c())
  for(c in 1:length(grouplist)+1){
    grouping1[c] = c();
  }
  mut_tab_col_names = colnames(patients_mutation_table)
  patients_row_names = rownames(patients_mutation_table)
  for(i in 1:length(grouplist)) { #for each group
    for(j in 1:length(dataset[,1])) { #for each patient
      for(z in 1:length(grouplist[[i]])){ #for each mutation in the group
        #split string to obtain name and type of the mutation
        temp = strsplit(grouplist[[i]][z],"_")
        group_mutation_name = temp[[1]][1]
        group_mutation_type = temp[[1]][2]
        #search index of mutation in mutations'table
        mutation_name_index = which(mutations[,2] == group_mutation_name & mutations[,1] == group_mutation_type)
        #convert name to ID
        mutation_id = rownames(mutations)[mutation_name_index]
        #find column relative to mutation ID
        col_index = which(mut_tab_col_names == mutation_id)
        #find row of matrix relative to patient j
        row_index = which(patients_row_names == rownames(dataset)[j])
        #(probably unusued) avoid multiple patients finded or no one finded
        if(length(row_index) != 0 & length(row_index) < 2) {
          #if the patients has the mutation add it to the group
          if(patients_mutation_table[row_index, col_index] == 1) {
            grouping1[i] = mapply(c, grouping1[i], j, SIMPLIFY = F)
          }
        }
      }
    }
  }
  #remove duplicates
  grouping1 = lapply(grouping1, function(x) x[!duplicated(x)])
  return(grouping1)
}

#function to make more useful dataset for mutations driven survival analysyis
# add each mutations as a column of the dataset and set 1 as value for each
# [patient, mutation] where the patient had the mutation
add_mutations_as_column = function(dataset, mutations, pat_mut_matrix) {
  #obtain row names and column names of the matrix
  matrix_patients = rownames(pat_mut_matrix)
  matrix_mutations = colnames(pat_mut_matrix)
  #remove pattern from mutations list
  mutations = subset(mutations, mutations[,1] != "Pattern")
  #for each mutation
  for(mut in 1:length(mutations[,1])) {
    temp = c()
    #for each patient
    for(pat in 1:length(matrix_patients)) {
      #if patients[pat] has the mutations[mut] --> 1
      if(pat_mut_matrix[matrix_patients[pat], rownames(mutations)[mut]] == 1) {
        temp = append(temp, 1)
      } else {
        temp = append(temp, 0)
      }
    }
    #add the column
    col_name = paste0(mutations[mut,2], "_")
    col_name = paste0(col_name, mutations[mut,1])
    dataset[col_name] = temp
  }
  return(dataset)
}

#function to make Cox Model for each group, obtain Survival and Baseline 
# Cumulative Hazard Function, test the model.
makeCoxModelsAndTests = function(dataset, patientsgroups, mutationsgroups) {
  #add ID column to subset the dataset in according to index in patientsgroups
  id = sort(sample(1:length(dataset$time)))
  dataset$id = id
  #define struct for each field (models, survival, hazards, residuals, ...)
  groupsmodels = list(c(), c())
  groupssurv = list(c(), c())
  groupsbchaz = list(c(), c())
  groupscsressurv = list(c(), c())
  groupsssres = list(c(), c())
  #for each group
  for(i in 1:length(mutationsgroups)) { 
    #subset data
    data = dataset[dataset$id %in% patientsgroups[[i]],]
    #create the formula for the Cox Model
    formula = paste0("SurvObj ~ ", mutationsgroups[[i]][1])
    for(j in 2:length(mutationsgroups[[i]])){ 
      formula = paste0(formula, " + ")
      formula = paste0(formula, mutationsgroups[[i]][j])
    }
    formula = as.formula(formula)
    #check if there is at least 1 event in the group's subset of the dataset
    events = FALSE;
    for(e in 1:length(data$status)) {
      if(data$status[e] == 1) {
        events = TRUE
      }
    }
    if(events) {
      #Make, if possible, the Cox Model
      cox = tryCatch({
        coxph(formula = formula, data = data)
      }, warning = function(w) {
        suppressWarnings(coxph(formula = formula, data = data))
      }, error = function(e) {         
        paste0("ERRORE : ", e, 
               "Una mutazione non è presente in nessun paziente") })
      groupsmodels[[i]] = cox
      #obtain Survival Function
      surv = tryCatch({
        survfit(coxph(formula = formula, data = data), conf.int = F)
      }, warning = function(w) {
        survfit(suppressWarnings(coxph(formula = formula, data = data)), 
                conf.int = F)
      }, error = function(e) { 
        paste0("ERRORE : ", e, 
               "Una mutazione non è presente in nessun paziente") })
      groupssurv[[i]] = surv
      #obtain Baseline Cumulative Hazard Function
      bchaz = tryCatch({
        basehaz(coxph(formula = formula, data = data))
      }, warning = function(w) {
        basehaz(suppressWarnings(coxph(formula = formula, data = data)))
      }, error = function(e) { paste0("ERRORE : ",e) })
      groupsbchaz[[i]] = bchaz
      #obtain Cox Snell Residuals Survival Function 
      csressurv = tryCatch({
        martingale_residuals = residuals(coxph(formula = formula, data = data), 
                                         type = 'martingale')
        cox_snell_residuals = data$status - martingale_residuals
        survfit(Surv(cox_snell_residuals, data$status) ~ 1)
      }, warning = function(w) {
        martingale_residuals = residuals(suppressWarnings(coxph(formula = formula,
                                                                data = data)), 
                                         type = 'martingale')
        cox_snell_residuals = data$status - martingale_residuals
        survfit(Surv(cox_snell_residuals, data$status) ~ 1)
      }, error = function(e) {         
        paste0("ERRORE : ", e, 
               "Una mutazione non è presente in nessun paziente") })
      groupscsressurv[[i]] = csressurv
      #obtain Scaled Schoenfeld Residuals
      ssres = tryCatch({
        suppressWarnings(cox.zph(coxph(formula = formula, data = data), 
                                 global = F, transform = 'rank'))
      }, warning = function(w) {
        suppressWarnings(cox.zph(suppressWarnings(coxph(formula = formula, 
                                                        data = data)), 
                                 global = F, transform = 'rank'))
      }, error = function(e) {         
        paste0("ERRORE : ", e, 
               "Una mutazione non è presente in nessun paziente") })
      groupsssres[[i]] = ssres
    } 
    else { #if there aren't events in the subset of the dataset
      groupsmodels[[i]] = "ERRORE: No events to make a Model"
      groupssurv[[i]] = c()
      groupsbchaz[[i]] = c()
      groupscsressurv[[i]] = c()
      groupsssres[[i]] = c()
    }
  }
  return(list("models" = groupsmodels, "survival" = groupssurv,
              "baseline_cumulative_hazard" = groupsbchaz,
              "cox_snell_residuals_surv_fun" = groupscsressurv,
              "scaled_schoenfeld_residuals" = groupsssres))
}

#plot Survival Function of a Cox Model
plotCoxSurvFun = function(surv, i = 0, type = "", col = "red") {
  opt_main = "\n"
  if(i != 0) {
    opt_main = paste0(" group ", i, opt_main)
  }
  if(type != "") {
    opt_main = paste0(" for ", type, opt_main)
  }
  #plot survival function
  plot(surv, col = col, mark.time = TRUE, lwd = 2,
       main = paste0("Survival Function of Cox Model", opt_main,
                     "with mutations as covariate"),
       xlab = "Time(Days)", ylab = "Survival %")
  grid(NULL, NULL, lty = 5)
}

#plot Baseline Cumulative Hazard function from a Cox Model 
plotCoxBslnCmHazFun = function(bchaz, i = 0, type = "MSS", col = "blue") {
  opt_main = ""
  if(i != 0) {
    opt_main = paste0(" group ", i, opt_main)
  }
  if(type != "") {
    opt_main = paste0(" for ",type , opt_main)
  }
  #plot cumulative hazard
  plot(bchaz$hazard ~ bchaz$time, 
       col = col, type = 's',
       xlab = "Time (Days)", ylab = "Cumulative Hazard",
       main = paste0("Baseline Cumulative Hazard Function", opt_main))
  grid(NULL, NULL, lty = 5)
} 

#plot Cox Snell Residuals Surv Function
plotCoxSnellResSurvFun = function(ressurv, i = 0, type = "", col = "blue") {
  opt_main = ""
  if(i != 0) {
    opt_main = paste0(" group ", i, opt_main)
  }
  if(type != "") {
    opt_main = paste0(" for ",type , opt_main)
  }
  plot(ressurv$time, -log(ressurv$surv), #cumulative hazard of residuals
       type = 'S', lwd = 1.5, col = col,
       main = paste0("Cox Snell Residuals", opt_main,
                     ":\n check overall fit of the model"),
       xlab = "residuals", ylab = "Cumulative Hazard of residuals")
  abline(0,1, col = 'red', lwd = 1.5, lty = 2) #reference line
  #Cumulative Hazard has to follow reference line 
  # (if it deviate too much some covariates doesn't respect PH assumption 
  #  and we need Scaled Schoenfeld Residuals to check covariates)
  legend('bottomright',
         cex = 0.75, lty = c(1,2), lwd = c(1.5, 1.5), 
         col = c('forestgreen', 'red'),
         legend = c('cumulative hazard', 'reference line'))
}

#function to plot Scaled Schoenfeld Residuals
plotScaledSchoenfeldResiduals = function(ssr, i = 0, type = "", col = "black") {
  opt_main = ""
  if(i != 0) { opt_main = paste0(" group ", i, opt_main) }
  if(type != "") { opt_main = paste0(" of ",type , opt_main) }
  ind = 1
  repeat{  #for each covariate
    mut = tryCatch({
      ssr[ind]
    }, error = function(e) {1}) #no more covariates
    if(is.numeric(mut)){
      break
    } 
    else {
      valid = TRUE
      #p_arr = parameter array = (rho, chisq, p-value)
      p_arr = eval(parse(text = (as.character(mut))[1])) 
      if(is.nan(p_arr[3])) { #if p-value is NaN --> covariate to inf (NA value in the model)
        valid = FALSE
      } 
      if(valid) {
        suppressWarnings(plot(mut, df = 2, lwd = 1.5, col = col,
                              main = paste0("Scaled Schoenfeld Residual test ", 
                                            opt_main)))
        #df = 2 --> degree of freedom --> linear test
        # so a covariate respect the proportional hazard assumption if is 
        # scaled Schoenfeld residual plots it is a horizontal line (more or less)
      } else {
        print(paste0("ERRORE: p-value NaN for ", rownames(mut$table), opt_main,
                     ": covariate to Inf. or zero residuals"))
      }
    }
    ind = ind + 1
  }
}

#add groups ad column of the dataset
addGroupsAsColumn = function(dataset, groups) {
  #add ID column to subset the dataset in according to index in patientsgroups
  id = sort(sample(1:length(dataset$time)))
  dataset$id = id
  for(i in 1:length(groups)){
    temp = c()
    for(j in 1:length(dataset$id)) {
      if(dataset$id[j] %in% groups[[i]]){
        temp = append(temp, 1)
      } else {
        temp = append(temp, 0)
      }
    }
    temp2 = paste0("group", i)
    dataset$temp2 = temp
    temprn = colnames(dataset)
    temprn[length(temprn)] = temp2
    colnames(dataset) = temprn
  }
  dataset$id = NULL
  return(dataset)
}

#function to find levels of mutation in according to graph's depth
find_levels_mutations = function(roots, nodes, edges) {
  #level1 = roots
  lvls = list(roots, c(), c()) #((level1), (level2), (level3))
  for(i in 1:length(roots)) {
    #all nodes reachable from a root --> level2
    arr = which(edges$from == roots[i])
    #skip logical nodes (OR/XOR/UPOR/UPXOR)
    arr = substitute_logicalnodes(arr, nodes, edges)
    for(j in 1:length(arr)) { 
      #level2 mutations
      node = edges$to[arr[j]]
      lvls[[2]] = append(lvls[[2]], node)
      #all nodes reachable from child 'arr' of a root --> level3
      arr2 = which(edges$from == node)
      if (length(arr2) > 0) {
        #skip logical nodes
        arr2 = substitute_logicalnodes(arr2, nodes, edges)
        for(k in 1:length(arr2)) {
          #level3 mutations
          lvls[[3]] = append(lvls[[3]], edges$to[arr2[k]])
        }
      }
    }
  }
  return(lvls)
}

substitute_logicalnodes = function(index_list, nodes, edges) {
  index_list2 = index_list
  for(a in 1:length(index_list)) {
    node_ind = edges$to[index_list[a]]
    node_name = nodes$name[which(nodes$id == node_ind)]
    if(node_name == "XOR" | node_name == "OR" | node_name == "UPXOR" | node_name == "UPOR") {
      new_ind = find_first_notlogical(node_ind, nodes, edges, c())
      index_list2 = index_list2[! index_list2 == index_list[a]]
      index_list2 = append(index_list2, new_ind)
    } 
  }
  return(index_list2)
}

#find first not logical nodes on the graph 
#new_in = index/es to replace instead of node_ind
find_first_notlogical = function(node_ind, nodes, edges, new_ind) {
  #find descendents of the logical node
  descendents = which(edges$from == node_ind)
  for(i in 1:length(descendents)) {
    #find node name
    node_name = nodes$name[which(nodes$id == edges$to[descendents[i]])]
    #if the descendent is still logical, recursively call the function
    if(node_name == "XOR" | node_name == "OR" | node_name == "UPXOR" | node_name == "UPOR"){
      new_ind = find_first_notlogical(edges$to[descendents[i]] ,nodes, edges, new_ind)
    } else {
      new_ind = append(new_ind, descendents[i])
    }
  }
  return(new_ind)
}

#function to confirm a patient in a lower level only if exist also in upper 
# level and then remove patients of lower level in upper lever
remove_other_levels_patients = function(dataset, levels, mut_levels,
                                        mutations, patients_mutation_table){
  #a patient is in the second level if he had a root mutations and a mutations
  # child of a root --> only if is also in level1
  levels[[2]] = levels[[2]][levels[[2]] %in% levels[[1]]]
  #a patient is in the third level if he had a root mutations, a child of a root
  # mutations and a child of a child of a a root mutations
  levels[[3]] = levels[[3]][levels[[3]] %in% levels[[2]]]
  levels[[3]] = levels[[3]][levels[[3]] %in% levels[[1]]]
  #Patients in lower level cannot be also in upper level (redundancy)
  # remove from level1 patients also in level2 and level3 
  levels[[1]] = levels[[1]][! levels[[1]] %in% levels[[2]]]
  levels[[1]] = levels[[1]][! levels[[1]] %in% levels[[3]]]
  # remove from level2 patients also in level3
  levels[[2]] = levels[[2]][! levels[[2]] %in% levels[[3]]]
  
  #first level contains also patients of possible level4, level5, ...
  # --> remove from level1 all the patients who have also other level mutations
  lev1_new = levels[[1]]
  pmt_mutations = colnames(patients_mutation_table)
  pmt_patients = rownames(patients_mutation_table)
  #remove all XOR_ OR_ ... column mutations from pat_mut_tab
  last_mut = 0
  for(j in 1:length(pmt_mutations)) {
    if(grepl("OR", pmt_mutations[j]) | grepl("XOR", pmt_mutations[j]) |
       grepl("UPOR", pmt_mutations[j]) | grepl("UPXOR", pmt_mutations[j])) {
      break;
    } else {
      last_mut = last_mut + 1;
    }
  }
  patients_mutation_table = patients_mutation_table[,1:last_mut]
  #remove patients
  for(k in 1:length(levels[[1]]))  {
    #find right patient 'i' row --> mutations raw
    mut_row = which(pmt_patients == rownames(dataset)[levels[[1]][k]])
    #count not root mutations for each patients
    not_root = 0
    for(z in 1:last_mut) {
      if(patients_mutation_table[mut_row, z] == 1) {
        mut_index = which(rownames(mutations) == pmt_mutations[z])
        mut_name = mutations[z,2]
        mut_type = mutations[z,1]
        mut = paste0(mut_name, '_', mut_type)
        if(! mut %in% mut_levels[[1]]) {
          not_root = not_root + 1
        }
      }
    }
    if(not_root > 0) {
      lev1_new = lev1_new[! lev1_new == levels[[1]][k]]
    }
  }
  #return levels full of patients
  levels[[1]] = lev1_new
  return(levels)
}

#remove patients from levelGrouping
#if a patient is more then 1 group check that is along a path
remove_patients_not_on_path = function(dataset, groups, levels, mutations,
                                       pm_matrix, nodes, edges) {
  #every patients in lev2 is still in lev1 (analog for lev3 and lev2)
  # --> remove a patient from lev2 is no one of him level2 mutations isn't
  #     on a path starting with one of is root mutations
  #     (analog for lev3 and lev2)
  groups = remove_lev(dataset, groups, levels, mutations, 
                      pm_matrix, nodes, edges, 2)
  groups = remove_lev(dataset, groups, levels, mutations,
                      pm_matrix, nodes, edges, 3)
  return(groups)
}

#remove from level2
remove_lev = function(dataset, groups, levels, mutations,
                       pm_matrix, nodes, edges, num) {
  #matrix row and column names (row = patients, col = mutations)
  pm_mutations = colnames(pm_matrix)
  pm_patients = rownames(pm_matrix)
  data_patients = rownames(dataset)
  #new result
  levgroup_new = groups[[num]]
  #for each patient
  for(i in 1:length(groups[[num]])) {
    #suppose not in level2, at the first path mutation founded
    levpresence = FALSE
    #find his level1 mutations 
    # find patient i in pm_matrix
    pm_patient_i = which(pm_patients == data_patients[groups[[num]][i]])
    #check which root mutations he has
    for(j in 1:length(levels[[num-1]])) {
      #for each mutation he has check if he has also the(/one of the) successive 
      # mutation(/mutations) on path (considering only 2nd level)
      #find relative mutation in matrix --> find index
      pm_mutations_i = which(paste0(mutations[,2], '_', mutations[,1]) == levels[[num-1]][j])
      if(pm_matrix[pm_patient_i, pm_mutations_i] == 1) { 
        #if patient has the mutation check if he has also one of the successive
        # mutations in the second level of the path
        desc = find_nextlevel_descendents(levels[[num-1]][j], nodes, edges)
        if(length(desc) > 0) {
          for(k in 1:length(desc)) {
            pm_mutations2_i = which(paste0(mutations[,2], '_', mutations[,1]) == desc[k])
            if(pm_matrix[pm_patient_i, pm_mutations2_i] == 1) {
              levpresence = TRUE
              break
            } #NO else needed
          }
        }
      }
    }
    if(levpresence == FALSE) {
      levgroup_new = levgroup_new[! levgroup_new == groups[[num]][i]]
    }
  }
  groups[[num]] = levgroup_new
  return(groups)
}

#finde next level mutations from a node at every depth
find_nextlevel_descendents = function(node, nodes, edges){
  #get node name and node type
  node_temp = strsplit(node, '_')
  node_name = node_temp[[1]][1]
  node_type = node_temp[[1]][2]
  #search node index in the nodes dataset
  node_id = nodes$id[which(nodes$name == node_name & nodes$type == node_type)]
  #find descendents 
  desc = which(edges$from == node_id)
  #result array
  desc2 = c()
  if(length(desc) > 0) {
    for(i in 1:length(desc)) {
      id = edges$to[desc[i]]
      name = nodes$name[which(nodes$id == id)]
      if(name == "XOR" | name == "OR" | name == "UPXOR" | name == "UPOR") {
        #replace with right next second level mutations
        mut = find_right_next_mut(id, nodes, edges)
      } else {
        type = nodes$type[which(nodes$id == id)]
        mut = paste0(name, '_', type)
      }
      desc2 = append(desc2, mut)
    }
  }
  return(desc2)
}

#replace logical node with true node
find_right_next_mut = function(id, nodes, edges) {
  #find descendents of the logical node
  descendents = which(edges$from == id)
  new_desc = c()
  #for each descendent
  for(i in 1:length(descendents)) {
    #find node name
    node_name = nodes$name[which(nodes$id == edges$to[descendents[i]])]
    #if the descendent is still logical
    if(node_name == "XOR" | node_name == "OR" | 
       node_name == "UPXOR" | node_name == "UPOR"){
      #recursively call the function
      new_desc = find_right_next_mut(edges$to[descendents[i]], nodes, edges)
    } else {
      node_type = nodes$type[which(nodes$id == edges$to[descendents[i]])]
      mut = paste0(node_name, '_', node_type)
      new_desc = append(new_desc, mut)
    }
  }
  return(new_desc)
}