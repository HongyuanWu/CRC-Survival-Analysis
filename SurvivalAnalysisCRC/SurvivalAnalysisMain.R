### SURVIVAL ANALYSIS ###
#Use TCGA-COADREAD samples to study Survival Analysis.


### LIBRARIES ###
#check existence of libraries
if(!require('devtools')) install.packages('devtools', dependencies = T, 
                                          repos='http://cran.us.r-project.org')
if(!require('TRONCO')) install_github("BIMIB-DISCo/TRONCO", ref = 'development')
if(!require('survival')) install.packages('survival', dependencies = T)
if(!require('muhaz')) install.packages('muhaz', dependencies = T)
if(!require('KMsurv')) install.packages('KMsurv', dependencies = T)
#import libraries needed
library(devtools) #for tronco
library(TRONCO) #for getting data
library(survival) #for survival analysis
library(muhaz) #for hazard function
library(KMsurv) #for life table


### IMPORT FUNCTIONS ###
#import function used in the program
source('Scripts/functions.R', echo = T)


###`IMPORT CLINICAL DATA ###
#get data from TCGA-COADREAD project (the same used for PiCnIc_COADREAD example)
source('Scripts/GettingData.R', echo = T) 
#it returns clinical.data 


### CREATE SURVIVAL OBJECT ###
#create survival object from the clinical.data
source('Scripts/CreateSurvivalObject.R', echo = T)
#it returns survival.data


### OBTAIN GROUPS TO SUBSET SURVIVAL DATA ###
#Get groups from PiCnIc graph of estimated progressions both for MSI and MSS
source('Scripts/XmlGraphParser.R', echo = T)
#it returns groups based on MSS/MSI and mutation progression's from PiCnIc 

### SUBSET CLINICAL DATA ###
#subset clinical data first by MSS and MSI-H (ignoring MSI-L) and then
# in according to precedent groups
source('Scripts/SubsettingData.R', echo = T)
#it returns subsets in according to precedent groups


### SURVIVAL ANALYSIS FOR MSI AND MSS ###
#Apply Survival Analysis models to MSI and MSS groups to estimate survival
# function, hazard function, ... and compare the groups
source('Scripts/MSI_StatusDrivenSurvivalAnalysis.R', echo = T)


### SURVIVAL ANALYSIS FOR EACH GROUP ###
#Apply Cox Model and Stratified Cox Model to estimate survival for each group
source('Scripts/MutationsDrivenSurvivalAnalysis.R', echo = T)

