## XML GRAPH PARSER
#input : MSI model and MSS model
#output : list of nodes and edges for each model
#algorithm : 
# - export the xml version for each graph
# - call a ruby script that parse the xml file and return 
#   csv lists for each model
# - import the csv files as dataset
# - create groups in according to lists of nodes and edges
#   - grouping1 : each groups contain a root and all his descendents
#   - grouping2 : each groups is a possible path (including alternative part 
#     of the path due to OR/AND/XOR) from a root to a leaf
#     (UNDONE : too few samples to do it)
#   - levelGrouping : group by level using 3 levels
#     - level 3 : child of child of a roots
#     - level 2 : child of a roots
#     - level 1 : roots
#     (level 1 patients CAN'T be also in level 2 or 3, so start assigining
#      patients from level 3, remove them from the dataset and go on bottom up)

#libraries
library(igraph)

#export the XML version for each graph
export.graphml(MSI.models, 'PiCnIcGraph/MSIgraph.xml')
export.graphml(MSS.models, 'PiCnIcGraph/MSSgraph.xml')

#call to a Ruby Script that create CSV file for nodes and edges for both models
system('Scripts/XmlGraphParser.rb', intern=TRUE)

#import the datasets
MSInodes = read.csv(file = 'CsvGraph/MSI/nodes.csv', header = FALSE, sep=",", 
                    col.names = c("id","name","type"))
MSIedges = read.csv(file = 'CsvGraph/MSI/edges.csv', header = FALSE, sep=",", 
                    col.names = c("from","to"))
MSSnodes = read.csv(file = 'CsvGraph/MSS/nodes.csv', header = FALSE, sep=",", 
                    col.names = c("id","name","type"))
MSSedges = read.csv(file = 'CsvGraph/MSS/edges.csv', header = FALSE, sep=",",
                    col.names = c("from","to"))
MSSnodes = data.frame(lapply(MSSnodes, as.character), stringsAsFactors=FALSE)
MSSedges = data.frame(lapply(MSSedges, as.character), stringsAsFactors=FALSE)
MSInodes = data.frame(lapply(MSInodes, as.character), stringsAsFactors=FALSE)
MSIedges = data.frame(lapply(MSIedges, as.character), stringsAsFactors=FALSE)

#find all the roots in the graph
# a node is a root if doesn't have incoming edge
# --> the node's id doesn't appear in the column "to" of edges' dataset
#MSS
MSSroots = c()
for(i in 1:length(MSSnodes$id)) {
  if(length(which(MSSedges$to == MSSnodes$id[i])) == 0)
     MSSroots = append(MSSroots, MSSnodes$id[i])
}
#MSI
MSIroots = c()
for(i in 1:length(MSInodes$id)) {
  if(length(which(MSIedges$to == MSInodes$id[i])) == 0)
    MSIroots = append(MSIroots, MSInodes$id[i])
}

#find all the discendent for each roots
# 1. from array of roots to array of groups 
#    (roots1, ..., rootsN) --> [[root1], ..., [rootN]]
# 2. for each root search descendents and add it to the group

#1. MSS
MSSgroups = list(c(),c())
for(i in 1:length(MSSroots)) {
  MSSgroups[i] = MSSroots[i]
}
#2. MSS
# - grouping1 --> only root and descendents for each group
MSSgroups = add_descendents(MSSgroups, MSSedges, MSSroots)
#Use grouping2 --> every group is a possible path (including alternatives cause
# of presence of OR,XOR,UPOR and UPXOR) from a root to a leaf
#MSSgroups2 = complete_path(MSSgroups, MSSedges, MSSroots)

#1. MSI
MSIgroups = list(c(),c()) 
for(i in 1:length(MSIroots)) {
  MSIgroups[i] = MSIroots[i]
} 
#2. MSI
#Use grouping1 --> only root and descendents for each group
MSIgroups = add_descendents(MSIgroups, MSIedges, MSIroots)
#Use grouping2 --> every group is a possible path (including alternatives cause
# of presence of OR,XOR,UPOR and UPXOR) from a root to a leaf
#MSIgroups2 = complete_path(MSIgroups, MSIedges, MSIroots)

#substitute mutation's id with mutation's name and remove duplicates, OR, XOR, 
# UPOR and UPXOR nodes 
MSSgroups_grouping1 = simplify_groups(MSSgroups, MSSnodes)
MSIgroups_grouping1 = simplify_groups(MSIgroups, MSInodes)
#MSSgroups_grouping2 = simplify_groups(MSSgroups2, MSSnodes)
#MSIgroups_grouping2 = simplify_groups(MSIgroups2, MSInodes)

#grouping2
#paths from each root --> UNDONE

#levelGrouping
#create the 3 levels and assign to each level the mutations
# only for MSS cause MSI is too much little
MSSlevels = find_levels_mutations(MSSroots, MSSnodes, MSSedges)
MSSlevels = simplify_groups(MSSlevels, MSSnodes)

#pathLevelGrouping
#create 3 levels, but populate it differently
#use levelGrouping and subset dataset in different way --> SubsettingData.R