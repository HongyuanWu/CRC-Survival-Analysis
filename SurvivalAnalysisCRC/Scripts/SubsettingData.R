## SUBSETTING DATA
#Subset survival.data : 
# 1) ignoring samples with "MSI_status" = "MSI-L" --> remove them
# 2) split into MSSsurvival.data and MSI-Hsurvival.data
# 3) remove sample who doesn't have both somatic and CNAs informations
#     (the ones who doesn't appear in both MAF and GISTIC)
#     OBS: done escludign samples that aren't in MSSlift and MSIlift
# 4) for both MSS and MSI dataset split in according to groups 
#     obtained from XmlGraphParser.

#1)removed patients with MSI-L as MSI_status in according to PiCnIc
survival.data = subset(survival.data, MSI_status != "MSI-L")
survival.data = subset(survival.data, MSI_status != "Not Evaluable")
survival.data$MSI_status = sapply(survival.data$MSI_status, as.character)
survival.data$MSI_status = sapply(survival.data$MSI_status, as.factor)

#2)Splitting into MSI and MSS
MSSsurvival.data = subset(survival.data, MSI_status == "MSS")
MSSsurvival.data$MSI_status = NULL
MSIsurvival.data = subset(survival.data, MSI_status == "MSI-H")
MSIsurvival.data$MSI_status = NULL

#3)Remove from survival.data every sample who doesn't have both somatic 
#and CNAs informations.
# --> remove patients in each dataset that doesn't appear in relative lift
MSSliftpatients = rownames(MSSlift$genotypes)
MSSsurvival.data = merge_samples(MSSsurvival.data, MSSliftpatients)
MSIliftpatients = rownames(MSIlift$genotypes)
MSIsurvival.data = merge_samples(MSIsurvival.data, MSIliftpatients)

#4)subset data in according to groups obtained from PiCnIc progression's graph
#GROUPING1 : 
# start from MSS and MSI groups precedently obtained and assign to each
# group all patients who have at least 1 mutation of the group
#MSS
MSSfullgroups = assign_patients_to_group(MSSsurvival.data, MSSgroups_grouping1, 
                                         MSSlift$annotations, #mutations
                                         MSSlift$genotypes) #matrix patient-mutation
#MSI
MSIfullgroups = assign_patients_to_group(MSIsurvival.data, MSIgroups_grouping1,
                                         MSIlift$annotations,#mutations
                                         MSIlift$genotypes)#matrix patient-mutation

#GROUPING2 : 
# each group is a path, add to a path every patient who have at least 
#UNDONE


#ONLY FOR MSS : 

#LEVEL GROUPING : 
#assign to level3 the patients who have a mutations which is a root, 
# a mutations child of a root and a mutations child of a child of a root.
#remove level3 patients from the dataset
#assign to level2 the patients who have a mutations which is a root and a 
# mutations that is a child of a root. (not include level3 patients cause they 
# were removed from the dataset)
#remove level2 patients from the dataset
#assign to level1 the patients who have a root mutations (not include level2 or
# level3 patients)
MSSlevelGroups_temp = assign_patients_to_group(MSSsurvival.data, MSSlevels,
                                          MSSlift$annotations, 
                                          MSSlift$genotypes)
MSSlevelGroups = remove_other_levels_patients(MSSsurvival.data, 
                                              MSSlevelGroups_temp,
                                              MSSlevels, MSSlift$annotations, 
                                              MSSlift$genotypes)

#PATH LEVEL GROUPING
#assign a patient to level1 if he has only a root mutation, to level2 if i has
# also a child mutation of the precedent root and add it to level3 if i has
# also a child mutation of the precedent child of the root
#patients assigned in precedent assign_patients_to_group invocation
#remove patients from level3 and level2 if they aren't on a path
MSSpathLevelGroups = remove_patients_not_on_path(MSSsurvival.data, 
                                                 MSSlevelGroups_temp,
                                                 MSSlevels, MSSlift$annotations, 
                                                 MSSlift$genotypes,
                                                 MSSnodes, MSSedges)
#remove redundancy as done before
MSSpathLevelGroups = remove_other_levels_patients(MSSsurvival.data, 
                                                  MSSpathLevelGroups,
                                                  MSSlevels, MSSlift$annotations, 
                                                  MSSlift$genotypes)
