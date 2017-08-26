##IMPORT DATA FOR SURVIVAL ANALYSIS
#import data from TCGA-COADREAD project.
#Suppose to have ALREADY run PiCnIc on the same samples and simply load 
# resulting data

#import all preprocessed PiCnIc data from TCGA-data
list_of_files = list.files(pattern="*.Rdata", recursive = TRUE)
for(i in 1:length(list_of_files)) {
  if(list_of_files[i] == "TCGA-data/MSI/Rdata-lifted/lifted.Rdata") {
    load(list_of_files[i])
    MSIlift = lift
  } else {
    if(list_of_files[i] == "TCGA-data/MSS/Rdata-lifted/lifted.Rdata") {
      load(list_of_files[i])
      MSSlift = lift
    } else {
      load(list_of_files[i])
    }
  }
}
remove(lift)

#create clinical.data fronm clinical.file
clinical.file = "TCGA-data/Clinical/crc_clinical_sheet.txt"
clinical.data = TCGA.map.clinical.data(
  file = clinical.file, 
  column.samples = 'patient', 
  column.map = c('days_to_death','days_to_last_followup','vital_status'))
#days_to_death OR days_to_last_followup : last time registered 
# (first case: DEAD, second case: CENSORING)
#vital_status = DECEASED or LIVING (--> censored)

#Import informations to sub-typing the CRC tumor (MSI-H/MSI-L/MSS)
cluster.data_temp = read.csv("TCGA-data/Clusters/TCGA-clusters.csv", 
                             header = TRUE, sep = ';')
cluster.data_temp = cluster.data_temp[1:length(clinical.data$days_to_death),]
names = data.frame(lapply(cluster.data_temp['patient'], 
                          as.character), stringsAsFactors=FALSE)
rownames(cluster.data_temp) = names[,1]
cluster.data = cluster.data_temp['MSI_status']
remove(cluster.data_temp)

#add MSI_status informations to clinical.data
clinical.data$MSI_status = cluster.data[,1]
remove(cluster.data)



