#!/usr/bin/env ruby

### XML GRAPH PARSER ###
#This script is call by the R script nomeancoradascegliere.R and it is used
# to parse a xml visualization of a graph generated with graphml
#It extracts all nodes and edges save it as nodes.csv and edges.csv,
# so the R script can used them to sub-grouping the clinical dataset.

#requires
require 'csv'

#function for parse the xml graph
def xml_parse(workdir, modelfile, csvdir)
  #import graph's xml file --> open a file handler/stream
  file = File.new(workdir+modelfile,"r")

  #read the imported file --> string, remove all \n and space between tags
  graph_string = file.read.gsub(/[\n]/, '').gsub(/    /, '').gsub(/  /, '')

  #modify graph's string to achieve desired format
  #find index of first node and first edge
  first_node_index = graph_string=~/<node/
  first_edge_index = graph_string=~/<edge/
  #last index for edges (last index for nodes = firt index for edges)
  last_edge_index = graph_string=~/<\/graph/
  #split the graph_string into nodes_string and edges_string
  nodes_string = graph_string[first_node_index..first_edge_index-1]
  edges_string = graph_string[first_edge_index..last_edge_index-1]
  #split at every node or edge --> 1 array of nodes and 1 array of edges
  nodes_array = nodes_string.split("</node>")
  edges_array = edges_string.split("</edge>")

  #extract only id and name from every node and put them into an
  # Hash Table --> [ID => NAME]
  i = 0
  nodes = Hash.new
  nodes_array.each do |node|
    node_useful_data = node.split("<data key=\"v_label\">")[0]
    temp = node_useful_data.gsub(/<node id=\"/,'')
    temp = temp.gsub(/\"><data key=\"v_name\">/,' ')
    temp = temp.gsub(/<\/data>/,'')
    temp = temp.split(' ')
    temp_value = Array.new
    temp_value << temp[1] #add node name
    temp_type = node.split("</data>")[2]
    type = temp_type.gsub(/<data key=\"v_type\">/, '')
    temp_value << type #add node type
    nodes[temp[0]]=temp_value #ID => NAME, TYPE
    #nodes[temp[0]]=temp[1]
    i+=1
  end

  #extract only "from ID" and "to ID" from every edge and put them into an
  # Hash Table --> [FROM => TO]
  #But hash table doesn't accept multiple edges from the same nodes
  # use Array of Hashes --> edges = {[FROM => TO], ..., [FROM => TO]}
  #edges = Hash.new
  edges = Array.new
  edges_array.each do |edge|
    edge_useful_data = edge.split(">")[0]
    temp = edge_useful_data.gsub(/<edge source=\"/, '')
    temp = temp.gsub(/\" target=\"/,' ')
    temp = temp.gsub(/\"/,'')
    temp = temp.split(' ')
    #edges[temp[0]]=temp[1]
    edges_temp = Hash.new
    edges_temp[temp[0]]=temp[1]
    edges << edges_temp
  end

  #close the file handler/stream
  file.close

  #export nodes and edges as csv
  CSV.open(workdir+csvdir+"nodes.csv", "wb") {
    |csv| nodes.to_a.each {|elem| csv << elem.flatten}
  }
  CSV.open(workdir+csvdir+"edges.csv", "wb") {
    |csv| edges.each {|elem| elem.to_a.each  {|elem| csv << elem} }
  }
end
#directories
surv_dir = "/home/mattia/Documenti/STAGE/R/SurvivalAnalysisCRC/"

#MAIN
xml_parse(surv_dir,"PiCnIcGraph/MSSgraph.xml", "CsvGraph/MSS/")
xml_parse(surv_dir,"PiCnIcGraph/MSIgraph.xml", "CsvGraph/MSI/")
