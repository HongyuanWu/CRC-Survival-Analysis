##CREATE SURVIVAL OBJECT
#create survival object from clinical.data
#A survival object is obtained from time (DEATH time or CENSORING time) and
# status (0 = LIVING, 1 = event occurred = DECEASED )

#Make days_to_death and days_to_last_followup numeric and vital_status chars
clinical.data$days_to_death = sapply(clinical.data$days_to_death, as.numeric)
clinical.data$days_to_last_followup = sapply(clinical.data$days_to_last_followup,
                                             as.numeric)
clinical.data$vital_status = sapply(clinical.data$vital_status, as.character)

#create 'time' replacing 'days_to_death' (represent death times) and 
# 'days_to_last_followup' (represent censoring times ONLY for the patients
# with days_to_death value = NA) with time column.
time = c() #time column
for (i in 1:length(clinical.data$days_to_death)) { 
  if(is.na(clinical.data$days_to_death[i])) { 
    #if days_to_death = NA --> censored data --> days_to_last_followup
    time = append(time, clinical.data$days_to_last_followup[i])
  } else { 
    #death --> days_to_death
    time = append(time, clinical.data$days_to_death[i])
  }
}

#create 'status' replacing 'vital_status' with a new column that has
# 1 for every time that 'vital_status = DECEASED' and 0 in the other case
status = c() #status column
for(i in 1:length(clinical.data$vital_status)) {
  if(identical(clinical.data$vital_status[i],'DECEASED')){
    status = append(status, 1)
  } else {
    status = append(status, 0)
  }
}

#create survival.data and modify it to obtain survival object
survival.data = clinical.data
#remove days_to_death, days_to_last_followup and vital_status
survival.data$days_to_death = NULL
survival.data$days_to_last_followup = NULL
survival.data$vital_status = NULL
#assign 'time' to survival.data
survival.data$time = time
#assign 'status' to survival.data
survival.data$status = status
#create Survival Object
SurvObj = Surv(survival.data$time, survival.data$status)
survival.data$SurvObj = SurvObj


