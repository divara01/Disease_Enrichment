####Generate omental network####
omental_network_genes <- read.csv("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/omental_network_genes.csv", na.strings="", stringsAsFactors=FALSE)
omental_network_genes$Direction <- NULL

omental_fat_mapped_nodes <- read.delim("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/omental_fat_mapped_nodes.txt", stringsAsFactors=FALSE)

omental_network_genes$From_EntrezID <- omental_fat_mapped_nodes$EntrezID[match(omental_network_genes$From, omental_fat_mapped_nodes$GeneSymbol)]

omental_network_genes$To_EntrezID <- omental_fat_mapped_nodes$EntrezID[match(omental_network_genes$To, omental_fat_mapped_nodes$GeneSymbol)]

omental_network_genes <- na.omit(omental_network_genes) #removes all entries that have NA for either from or to
omental_Gene_IDs_Index <- omental_network_genes

omental_network_genes$From <- NULL
omental_network_genes$To <- NULL

relations_whole_omental <- data.frame(from = omental_network_genes$From_EntrezID, to = omental_network_genes$To_EntrezID)
g_whole_omental <- graph.data.frame(omental_network_genes, directed = TRUE)
vertices_whole_omental <- V(g_whole_omental)

####Generate MM network####
magic_module_genes <- read.delim("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/magic_module_genes_updated.txt", na.strings="", stringsAsFactors=FALSE)
magic_module_genes$EntrezID <- as.character(magic_module_genes$EntrezID)
MM_network_genes <- V(g_whole_omental)$name[which((V(g_whole_omental)$name %in% magic_module_genes$EntrezID))]

MM_network_pos <- match(MM_network_genes, V(g_whole_omental)$name)

g_MM <- induced.subgraph(g_whole_omental, MM_network_pos, impl = "auto")

length(V(g_MM)$name) #458
length(which(!is.na(match(V(g_whole_omental)$name, magic_module_genes$EntrezID)))) #458 magic genes in omental network (with EntrezIds)

####Generate residual omental (RO) network####
RO_network_genes <- V(g_whole_omental)$name[which(!(V(g_whole_omental)$name %in% magic_module_genes$EntrezID))]

RO_network_pos <- match(RO_network_genes, V(g_whole_omental)$name)

g_RO <- induced.subgraph(g_whole_omental, RO_network_pos, impl = "auto")

length(V(g_RO)$name) #4847

####Situation #1-4: Get X neighbors of each gene B subnetwork in RO####
#(Function [GenerateOmentalSubnetworks_InfoFlowOmentalMM] is saved in Functions.Rmd file)

mode_magic_list <- c("in", "out")
mode_neighbors_list <- c("in", "out", "all")
magic_module <- magic_module_genes$EntrezID

for (num_neighbors in 1:5){
  for (mode_magic_index in 1:2){
    mode_magic <- mode_magic_list[mode_magic_index]
    for (mode_neighbors_index in 1:3){
      mode_neighbors <- mode_neighbors_list[mode_neighbors_index]
      functionOutput <- GenerateOmentalSubnetworks_InfoFlowOmentalMM(num_neighbors, mode_magic, mode_neighbors, magic_module)
      temp_df <- functionOutput[[1]]
      RO_neighbors_list <- functionOutput[[2]]
      #Assign relevant title for each data frame
      direction_magic <- ifelse(mode_magic == "out", "FROM", ifelse(mode_magic == "in", "TO", "ALL"))
      direction_neighbor <- ifelse(mode_neighbors == "out", "FROM", ifelse(mode_neighbors == "in", "TO", "ALL"))
      neighbors <- c("C", "D", "E", "F", "G")
      
      colnames(temp_df) <- paste(direction_neighbor, "_", colnames(temp_df), sep = "")
      
      assign(paste("B_", direction_neighbor, "_", neighbors[num_neighbors], "_", direction_magic, "_MM", sep = ""), temp_df)
      rm(temp_df)
    }
  }
}

#cbind all columns to create a summary table
Summary_B_C_FROM_MM <- cbind(B_FROM_C_FROM_MM, B_TO_C_FROM_MM, B_ALL_C_FROM_MM)
Summary_B_D_FROM_MM <- cbind(B_FROM_D_FROM_MM, B_TO_D_FROM_MM, B_ALL_D_FROM_MM)
Summary_B_E_FROM_MM <- cbind(B_FROM_E_FROM_MM, B_TO_E_FROM_MM, B_ALL_E_FROM_MM)
Summary_B_F_FROM_MM <- cbind(B_FROM_F_FROM_MM, B_TO_F_FROM_MM, B_ALL_F_FROM_MM)
Summary_B_G_FROM_MM <- cbind(B_FROM_G_FROM_MM, B_TO_G_FROM_MM, B_ALL_G_FROM_MM)

Summary_B_C_TO_MM <- cbind(B_FROM_C_TO_MM, B_TO_C_TO_MM, B_ALL_C_TO_MM)
Summary_B_D_TO_MM <- cbind(B_FROM_D_TO_MM, B_TO_D_TO_MM, B_ALL_D_TO_MM)
Summary_B_E_TO_MM <- cbind(B_FROM_E_TO_MM, B_TO_E_TO_MM, B_ALL_E_TO_MM)
Summary_B_F_TO_MM <- cbind(B_FROM_F_TO_MM, B_TO_F_TO_MM, B_ALL_F_TO_MM)
Summary_B_G_TO_MM <- cbind(B_FROM_G_TO_MM, B_TO_G_TO_MM, B_ALL_G_TO_MM)

####Situation 5-8: Get X neighbors of each gene B subnetwork in Magic####
#(Function [GenerateMagicSubnetworks_InfoFlowMMOmental] is saved in Functions.Rmd file)

mode_RO_list <- c("in", "out")
mode_neighbors_list <- c("in", "out", "all")

for (num_neighbors in 1:5){
  for (mode_RO_index in 1:2){
    mode_RO <- mode_RO_list[mode_RO_index]
    for (mode_neighbors_index in 1:3){
      mode_neighbors <- mode_neighbors_list[mode_neighbors_index]
      functionOutput <- GenerateMagicSubnetworks_InfoFlowMMOmental(num_neighbors, mode_RO, mode_neighbors)
      temp_df <- functionOutput[[1]]
      MM_neighbors_list <- functionOutput[[2]]
      direction_RO <- ifelse(mode_RO == "out", "TO", ifelse(mode_RO == "in", "FROM", "ALL"))
      direction_neighbor <- ifelse(mode_neighbors == "out", "FROM", ifelse(mode_neighbors == "in", "TO", "ALL"))
      neighbors <- c("P", "Q", "R", "S", "T")
      
      colnames(temp_df) <- paste(direction_neighbor, "_", colnames(temp_df), sep = "")
      
      assign(paste("A_", direction_neighbor, "_", neighbors[num_neighbors], "_", direction_RO, "_RO", sep = ""), temp_df)
      rm(temp_df)
    }
  }
}

#cbind all columns to create a summary table
Summary_A_P_FROM_RO <- cbind(A_FROM_P_FROM_RO, A_TO_P_FROM_RO, A_ALL_P_FROM_RO)
Summary_A_Q_FROM_RO <- cbind(A_FROM_Q_FROM_RO, A_TO_Q_FROM_RO, A_ALL_Q_FROM_RO)
Summary_A_R_FROM_RO <- cbind(A_FROM_R_FROM_RO, A_TO_R_FROM_RO, A_ALL_R_FROM_RO)
Summary_A_S_FROM_RO <- cbind(A_FROM_S_FROM_RO, A_TO_S_FROM_RO, A_ALL_S_FROM_RO)
Summary_A_T_FROM_RO <- cbind(A_FROM_T_FROM_RO, A_TO_T_FROM_RO, A_ALL_T_FROM_RO)

Summary_A_P_TO_RO <- cbind(A_FROM_P_TO_RO, A_TO_P_TO_RO, A_ALL_P_TO_RO)
Summary_A_Q_TO_RO <- cbind(A_FROM_Q_TO_RO, A_TO_Q_TO_RO, A_ALL_Q_TO_RO)
Summary_A_R_TO_RO <- cbind(A_FROM_R_TO_RO, A_TO_R_TO_RO, A_ALL_R_TO_RO)
Summary_A_S_TO_RO <- cbind(A_FROM_S_TO_RO, A_TO_S_TO_RO, A_ALL_S_TO_RO)
Summary_A_T_TO_RO <- cbind(A_FROM_T_TO_RO, A_TO_T_TO_RO, A_ALL_T_TO_RO)
