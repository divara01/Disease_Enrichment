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

####Identify MM -> RO SN connections (situation 1+6)####
#situation 6: Use FROM file for list of A's.  Use TO numbers for SNs
#situation 1: Use FROM file for list of B's.  Use FROM numbers for SNs
listofAs <- vector()
for (i in 1:nrow(Summary_A_P_FROM_RO)){
  unlistRowName <- unlist(strsplit(row.names(Summary_A_P_FROM_RO)[i], split='_', fixed=TRUE))
  listofAs <- append(listofAs, unlistRowName[1])
}

listofBs <- vector()
for (i in 1:nrow(Summary_B_C_FROM_MM)){
  unlistRowName <- unlist(strsplit(row.names(Summary_B_C_FROM_MM)[i], split='_', fixed=TRUE))
  listofBs <- append(listofBs, unlistRowName[1])
}

listofAs_pos <- match(listofAs, V(g_whole_omental)$name)
RO_from_MM <- neighborhood(g_whole_omental, order = 1, nodes = listofAs_pos, mode = "out")

for (list_index in 1: length(RO_from_MM)){
  RO_from_MM[[list_index]] <- V(g_whole_omental)$name[RO_from_MM[[list_index]]]
}
rm(list_index)

RO_FROM_MM_gene_list <- data.frame(matrix(NA, nrow = length(RO_from_MM), ncol = 13))
colnames(RO_FROM_MM_gene_list) <- c("A", "Bs", "num_Bs", "num_Cs", "num_Ds", "num_Es", "num_Fs", "num_Gs", "num_Ps", "num_Qs", "num_Rs", "num_Ss", "num_Ts")

for(list_index in 1:length(RO_from_MM)){
  RO_FROM_MM_gene_list[list_index, 1] <- RO_from_MM[[list_index]][1]
  matched_As <- as.character(na.omit(listofAs[match(RO_from_MM[[list_index]][1], listofAs)]))
  matched_Bs <- as.character(na.omit(listofBs[match(RO_from_MM[[list_index]], listofBs)]))
  A_row_pos <- match(matched_As, listofAs)
  B_row_pos <- match(matched_Bs, listofBs)
  RO_FROM_MM_gene_list$Bs[list_index] <- paste(matched_Bs, sep = ",", collapse = ", ")
  RO_FROM_MM_gene_list$num_Bs[list_index] <- length(matched_Bs)
  RO_FROM_MM_gene_list$num_Cs[list_index] <- paste(Summary_B_C_FROM_MM$FROM_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_FROM_MM_gene_list$num_Ds[list_index] <- paste(Summary_B_D_FROM_MM$FROM_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_FROM_MM_gene_list$num_Es[list_index] <- paste(Summary_B_E_FROM_MM$FROM_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_FROM_MM_gene_list$num_Fs[list_index] <- paste(Summary_B_F_FROM_MM$FROM_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_FROM_MM_gene_list$num_Gs[list_index] <- paste(Summary_B_G_FROM_MM$FROM_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  
  RO_FROM_MM_gene_list$num_Ps[list_index] <- Summary_A_P_FROM_RO$TO_num_genes_in_SN[A_row_pos]
  RO_FROM_MM_gene_list$num_Qs[list_index] <- Summary_A_Q_FROM_RO$TO_num_genes_in_SN[A_row_pos]
  RO_FROM_MM_gene_list$num_Rs[list_index] <- Summary_A_R_FROM_RO$TO_num_genes_in_SN[A_row_pos]
  RO_FROM_MM_gene_list$num_Ss[list_index] <- Summary_A_S_FROM_RO$TO_num_genes_in_SN[A_row_pos]
  RO_FROM_MM_gene_list$num_Ts[list_index] <- Summary_A_T_FROM_RO$TO_num_genes_in_SN[A_row_pos]
}

####Identify Each B and P for each A (situation 1+6)####
listofAs <- vector()
for (i in 1:nrow(Summary_A_P_FROM_RO)){
  unlistRowName <- unlist(strsplit(row.names(Summary_A_P_FROM_RO)[i], split='_', fixed=TRUE))
  listofAs <- append(listofAs, unlistRowName[1])
}

pos_As <- match(listofAs, V(g_whole_omental)$name)

neighbors_of_As_Bs <- neighborhood(g_whole_omental, order = 1, nodes = pos_As, mode = "out")

lengths <- NULL
for (list_index in 1: length(neighbors_of_As_Bs)){
  neighbors_of_As_Bs[[list_index]] <- V(g_whole_omental)$name[neighbors_of_As_Bs[[list_index]]]
  neighbors_of_As_Bs[[list_index]] <- as.vector(na.omit(V(g_RO)$name[match(neighbors_of_As_Bs[[list_index]], V(g_RO)$name)]))
  lengths <- append(lengths, length(neighbors_of_As_Bs[[list_index]]))
}

neighbors_of_As_Ps <- neighborhood(g_whole_omental, order = 1, nodes = pos_As, mode = "in")

for (list_index in 1: length(neighbors_of_As_Ps)){
  neighbors_of_As_Ps[[list_index]] <- V(g_whole_omental)$name[neighbors_of_As_Ps[[list_index]]]
  neighbors_of_As_Ps[[list_index]] <- as.vector(na.omit(V(g_MM)$name[match(neighbors_of_As_Ps[[list_index]], V(g_MM)$name)]))
  lengths <- append(lengths, length(neighbors_of_As_Ps[[list_index]]))
}
rm(list_index)
SNs_Situation1 <- data.frame(matrix(NA, nrow = max(lengths), ncol = length(listofAs)*2))

cols <- NULL
for(i in 1:ncol(SNs_Situation1)){
  var <- ceiling(i/2)
  sub <- ifelse(i %% 2 == 0, "B", "P")
  cols <- append(cols, paste(listofAs[var], sub, sep = "_"))
}
colnames(SNs_Situation1) <- cols
rm(i, var, sub, cols, lengths)

for (count in 1:ncol(SNs_Situation1)){
  val <- count/2
  index <- ceiling(val)
  if (count %% 2 == 0){
    Blength <- length(neighbors_of_As_Bs[[index]])
    for(length in 1:Blength){
      SNs_Situation1[length, count] <- neighbors_of_As_Bs[[index]][length]
    }
  } else{
    Plength <- length(neighbors_of_As_Ps[[index]])
    for(length in 1:Plength){
      SNs_Situation1[length, count] <- neighbors_of_As_Ps[[index]][length]
    }
  } 
}
rm(val, index, count, Blength, Plength)

####Generate Subnetworks for each B and P to each degree (situation 1+6)####
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation1)){
    if(i %% 2 == 0){
      temp_listofBs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "out")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
        length_each_B_neighbors <- length(temp_neighbors_of_Bs[[list_index]])
        for (row in 1:length_each_B_neighbors){
          temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_neighbors_of_Bs[[list_index]][row]
        }
        colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Bs[[list_index]][1], sep = "_")
        temp_df_B_SNs$placeholder <- NA
      }
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "in")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  assign(paste("Situation1_B", neighbors_list[num_neighbors], "SNs", sep = "_"), temp_df_B_SNs)
  assign(paste("Situation1_A", neighbors_list[num_neighbors+5], "SNs", sep = "_"), temp_df_P_SNs)
}

rm(temp_df_P_SNs, temp_df_B_SNs, delete, lengths_P, lengths_B, temp_pos_Bs, temp_pos_Ps)

####Generate SN combining P and each B (situation 1+6)####
#(Function [Generate.MMandROCombined.SN ()] is saved in Functions.Rmd file)
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation1)){
    if(i %% 2 == 0){
      temp_listofBs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "out")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
        length_each_B_neighbors <- length(temp_neighbors_of_Bs[[list_index]])
        for (row in 1:length_each_B_neighbors){
          temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_neighbors_of_Bs[[list_index]][row]
        }
        colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Bs[[list_index]][1], sep = "_")
        temp_df_B_SNs$placeholder <- NA
      }
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "in")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  PAB_Individual_Bs_Subnetworks <- Generate.MMandROCombined.SN(temp_df_B_SNs, temp_df_P_SNs, num_neighbors)
  
  lengths_PAB <- NULL
  for(i in 1:ncol(PAB_Individual_Bs_Subnetworks)){
    lengths_PAB <- append(lengths_PAB, length(which(!is.na(PAB_Individual_Bs_Subnetworks[,i]))))
  }
  delete <- max(lengths_PAB) + 1:nrow(PAB_combined_Subnetworks)
  PAB_Individual_Bs_Subnetworks <- PAB_Individual_Bs_Subnetworks[-delete, ]
  
  filename <- paste("Situation1_PAB_Individual_Bs_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(PAB_combined_Subnetworks, file = filename)
  
  assign(paste("Situation1_B", neighbors_list[num_neighbors], "SNs", sep = "_"), temp_df_B_SNs)
  assign(paste("Situation1_A", neighbors_list[num_neighbors+5], "SNs", sep = "_"), temp_df_P_SNs)
  assign(paste("Situation1_PAB_Individual_Bs_SN_d_", num_neighbors, sep = ""), PAB_Individual_Bs_Subnetworks)
}

####Identify RO -> MM SN connections (situation 2+5)####

#Situation 5: use TO FILE to determine A's, TO number for genes in SN
#Situation 2: use TO FILE to determine B's, FROM number for genes in SN

listofAs <- vector()
for (i in 1:nrow(Summary_A_P_TO_RO)){
  unlistRowName <- unlist(strsplit(row.names(Summary_A_P_TO_RO)[i], split='_', fixed=TRUE))
  listofAs <- append(listofAs, unlistRowName[1])
}

listofBs <- vector()
for (i in 1:nrow(Summary_B_C_TO_MM)){
  unlistRowName <- unlist(strsplit(row.names(Summary_B_C_TO_MM)[i], split='_', fixed=TRUE))
  listofBs <- append(listofBs, unlistRowName[1])
}

listofAs_pos <- match(listofAs, V(g_whole_omental)$name)
RO_to_MM <- neighborhood(g_whole_omental, order = 1, nodes = listofAs_pos, mode = "in")

for (list_index in 1: length(RO_to_MM)){
  RO_to_MM[[list_index]] <- V(g_whole_omental)$name[RO_to_MM[[list_index]]]
}
rm(list_index)

RO_TO_MM_gene_list <- data.frame(matrix(NA, nrow = length(RO_to_MM), ncol = 13))
colnames(RO_TO_MM_gene_list) <- c("A", "Bs", "num_Bs", "num_Cs", "num_Ds", "num_Es", "num_Fs", "num_Gs", "num_Ps", "num_Qs", "num_Rs", "num_Ss", "num_Ts")

for(list_index in 1:length(RO_to_MM)){
  RO_TO_MM_gene_list[list_index, 1] <- RO_to_MM[[list_index]][1]
  matched_As <- as.character(na.omit(listofAs[match(RO_to_MM[[list_index]][1], listofAs)]))
  matched_Bs <- as.character(na.omit(listofBs[match(RO_to_MM[[list_index]], listofBs)]))
  A_row_pos <- match(matched_As, listofAs)
  B_row_pos <- match(matched_Bs, listofBs)
  RO_TO_MM_gene_list$Bs[list_index] <- paste(matched_Bs, sep = ",", collapse = ", ")
  RO_TO_MM_gene_list$num_Bs[list_index] <- length(matched_Bs)
  RO_TO_MM_gene_list$num_Cs[list_index] <- paste(Summary_B_C_TO_MM$TO_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_TO_MM_gene_list$num_Ds[list_index] <- paste(Summary_B_D_TO_MM$TO_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_TO_MM_gene_list$num_Es[list_index] <- paste(Summary_B_E_TO_MM$TO_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_TO_MM_gene_list$num_Fs[list_index] <- paste(Summary_B_F_TO_MM$TO_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  RO_TO_MM_gene_list$num_Gs[list_index] <- paste(Summary_B_G_TO_MM$TO_num_genes_in_SN[B_row_pos], sep = ",", collapse = ", ")
  
  RO_TO_MM_gene_list$num_Ps[list_index] <- Summary_A_P_TO_RO$FROM_num_genes_in_SN[A_row_pos]
  RO_TO_MM_gene_list$num_Qs[list_index] <- Summary_A_Q_TO_RO$FROM_num_genes_in_SN[A_row_pos]
  RO_TO_MM_gene_list$num_Rs[list_index] <- Summary_A_R_TO_RO$FROM_num_genes_in_SN[A_row_pos]
  RO_TO_MM_gene_list$num_Ss[list_index] <- Summary_A_S_TO_RO$FROM_num_genes_in_SN[A_row_pos]
  RO_TO_MM_gene_list$num_Ts[list_index] <- Summary_A_T_TO_RO$FROM_num_genes_in_SN[A_row_pos]
}

####Identify Each B and P for each A (situation 2+5)####
listofAs <- vector()
for (i in 1:nrow(Summary_A_P_TO_RO)){
  unlistRowName <- unlist(strsplit(row.names(Summary_A_P_TO_RO)[i], split='_', fixed=TRUE))
  listofAs <- append(listofAs, unlistRowName[1])
}

pos_As <- match(listofAs, V(g_whole_omental)$name)

neighbors_of_As_Bs <- neighborhood(g_whole_omental, order = 1, nodes = pos_As, mode = "in")

lengths <- NULL
for (list_index in 1: length(neighbors_of_As_Bs)){
  neighbors_of_As_Bs[[list_index]] <- V(g_whole_omental)$name[neighbors_of_As_Bs[[list_index]]]
  neighbors_of_As_Bs[[list_index]] <- as.vector(na.omit(V(g_RO)$name[match(neighbors_of_As_Bs[[list_index]], V(g_RO)$name)]))
  lengths <- append(lengths, length(neighbors_of_As_Bs[[list_index]]))
}

neighbors_of_As_Ps <- neighborhood(g_whole_omental, order = 1, nodes = pos_As, mode = "out")

for (list_index in 1: length(neighbors_of_As_Ps)){
  neighbors_of_As_Ps[[list_index]] <- V(g_whole_omental)$name[neighbors_of_As_Ps[[list_index]]]
  neighbors_of_As_Ps[[list_index]] <- as.vector(na.omit(V(g_MM)$name[match(neighbors_of_As_Ps[[list_index]], V(g_MM)$name)]))
  lengths <- append(lengths, length(neighbors_of_As_Ps[[list_index]]))
}
rm(list_index)
SNs_Situation2 <- data.frame(matrix(NA, nrow = max(lengths), ncol = length(listofAs)*2))

cols <- NULL
for(i in 1:ncol(SNs_Situation2)){
  var <- ceiling(i/2)
  sub <- ifelse(i %% 2 == 0, "B", "P")
  cols <- append(cols, paste(listofAs[var], sub, sep = "_"))
}
colnames(SNs_Situation2) <- cols
rm(i, var, sub, cols, lengths)

for (count in 1:ncol(SNs_Situation2)){
  val <- count/2
  index <- ceiling(val)
  if (count %% 2 == 0){
    Blength <- length(neighbors_of_As_Bs[[index]])
    for(length in 1:Blength){
      SNs_Situation2[length, count] <- neighbors_of_As_Bs[[index]][length]
    }
  } else{
    Plength <- length(neighbors_of_As_Ps[[index]])
    for(length in 1:Plength){
      SNs_Situation2[length, count] <- neighbors_of_As_Ps[[index]][length]
    }
  } 
}
rm(val, index, count, Blength, Plength)

####Generate Subnetworks for each B and P to each degree (situation 2+5)####
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation2)){
    if(i %% 2 == 0){
      i = 4
      temp_listofBs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "in")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
        length_each_B_neighbors <- length(temp_neighbors_of_Bs[[list_index]])
        for (row in 1:length_each_B_neighbors){
          temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_neighbors_of_Bs[[list_index]][row]
        }
        colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation2)[i], temp_neighbors_of_Bs[[list_index]][1], sep = "_")
        temp_df_B_SNs$placeholder <- NA
      }
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "out")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation2)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  assign(paste("Situation2_B", neighbors_list[num_neighbors], "SNs", sep = "_"), temp_df_B_SNs)
  assign(paste("Situation2_A", neighbors_list[num_neighbors+5], "SNs", sep = "_"), temp_df_P_SNs)
}

rm(temp_df_P_SNs, temp_df_B_SNs, delete, lengths_P, lengths_B, temp_pos_Bs, temp_pos_Ps)

####Generate SN combining P and each B (situation2+5)####
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation2)){
    if(i %% 2 == 0){
      temp_listofBs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "in")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
        length_each_B_neighbors <- length(temp_neighbors_of_Bs[[list_index]])
        for (row in 1:length_each_B_neighbors){
          temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_neighbors_of_Bs[[list_index]][row]
        }
        colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation2)[i], temp_neighbors_of_Bs[[list_index]][1], sep = "_")
        temp_df_B_SNs$placeholder <- NA
      }
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "out")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation2)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  BAP_Individual_Ps_Subnetworks <- Generate.MMandROCombined.SN(temp_df_B_SNs, temp_df_P_SNs, num_neighbors)
  
  lengths_BAP <- NULL
  for(i in 1:ncol(BAP_Individual_Ps_Subnetworks)){
    lengths_BAP <- append(lengths_BAP, length(which(!is.na(BAP_Individual_Ps_Subnetworks[,i]))))
  }
  delete <- max(lengths_BAP) + 1:nrow(BAP_Individual_Ps_Subnetworks)
  BAP_Individual_Ps_Subnetworks <- BAP_Individual_Ps_Subnetworks[-delete, ]
  
  filename <- paste("Situation2_BAP_Individual_Ps_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(BAP_Individual_Ps_Subnetworks, file = filename)
  
  assign(paste("Situation2_B", neighbors_list[num_neighbors], "SNs", sep = "_"), temp_df_B_SNs)
  assign(paste("Situation2_A", neighbors_list[num_neighbors+5], "SNs", sep = "_"), temp_df_P_SNs)
  assign(paste("Situation2_BAP_Individual_Ps_SN_d_", num_neighbors, sep = ""), BAP_Individual_Ps_Subnetworks)
}

####Disease Enrichment on Individual Bs and Ps Subnetworks####
filename = "DiseaseEnrichment_Situation2_BAP_Individual_Ps_SN_d_5.csv"
subnetworks <- Situation2_BAP_Individual_Ps_SN_d_5 #NEED TO CHANGE FILE DEGREE HERE

output <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2))

numSN <- length(colnames(subnetworks))
numSamples <- 100 #(enter number of sample per group)
numDiseases <- length(colnames(GeneSetTable))

TableauDiseaseEnrichment <- data.frame(matrix(NA, nrow = 1, ncol = 6))
colnames(TableauDiseaseEnrichment) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")

Disease_sample_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) * numSamples, ncol = ncol(subnetworks)), row.names = NULL)
colnames(Disease_sample_p_values) <- colnames(subnetworks)

output2 <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2, ncol = ncol(subnetworks) * 2))

rows <- NULL
for (i in 1:nrow(output2)){
  var <- ceiling(i/2)
  sub <- ifelse(i %% 2 == 0, "FDR", "pval")
  rows <- append(rows, paste(colnames(GeneSetTable[var]), sub, sep = "_"))
}
rownames(output2) <- rows
rm(rows)

cols <- NULL
for (i in 1:ncol(output2)){
  var <- ceiling(i/2)
  sub <- ifelse(i %% 2 == 0, "b", "a")
  cols <- append(cols, paste(colnames(subnetworks[var]), sub, sep = "_"))
}
colnames(output2) <- cols
rm(cols)

for (k in 1:numSN) {
  for (j in 1:numDiseases) {
    for (i in 1:numSamples) {
      TargetSize <- length(which(!is.na(subnetworks[, k])))
      if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(RO_network_genes[match(GeneSetTable[, j], sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL))]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(GeneSetTable[match(sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
      }
      TestSetSize <- length(which(!is.na(match(RO_network_genes, GeneSetTable[,j]))))
      TotalSize <- length(RO_network_genes)
      
      val <- (j * numSamples) - (numSamples - i)
      
      temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
      
      Disease_sample_p_values[val, k] <- temp_list[[1]]
    }
  }
}

for (i in 1:numSN) {
  Disease_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) , ncol = 6), row.names = NULL)
  colnames(Disease_p_values) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")
  
  for (j in 1:numDiseases) {
    
    if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,i]))) {
      Overlap <- length(which(!is.na(subnetworks[match(GeneSetTable[, j], as.character(subnetworks[, i])), i]))) 
    }
    if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,i]))){
      Overlap <- length(which(!is.na(GeneSetTable[match(as.character(subnetworks[, i]), GeneSetTable[, j]), j]))) 
    }
    
    TargetSize <- length(which(!is.na(subnetworks[, i])))
    TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_RO)$name))))
    
    TotalSize <- length(V(g_RO)$name)
    index <- (i * numSN) - (numSN - j) 
    temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
    Disease_p_values[j, 1] <- colnames(GeneSetTable)[j]
    Disease_p_values[j, 2] <- TestSetSize
    Disease_p_values[j, 3] <- paste(unlist(strsplit(colnames(subnetworks)[i], split='_', fixed=TRUE))[1:4], collapse = "_")
    Disease_p_values[j, 4] <- length(na.omit(subnetworks[,i]))
    Disease_p_values[j, 5] <- temp_list[[1]]
    row <- (2 * j) - 1
    col <- (2 * i) - 1
    output2[row, col] <- temp_list[[1]]
    row2 <- (2 * j) - 1
    col2 <- 2 * i
    output2[row2, col2] <- Overlap
    row3 <- 2 * j
    col3 <- 2 * i
    output2[row3, col3] <- TestSetSize
    
    #FDR = a/(b*c)
    #where: a = # of actual pvalues for that KD that are less than the actual p-values in that KD; b = # of sample p-values for that KD that are less than the actual p-value in that KD * number of samples generated per pathway per KD
    test_pval <- Disease_p_values[j, 5]
    
    if (test_pval <= 0.05){
      a <- length(which(Disease_p_values[, 5] <= test_pval))
      b <- length(which(Disease_sample_p_values[, i] <= test_pval))
      c <- numSamples
      if (b != 0){
        row = 2 * j
        col = (2 * i) - 1
        output2[row, col] <- b / (a * c)
        Disease_p_values[j, 6] <- b / (a * c)
      }
      if (b == 0){
        row = 2 * j
        col = (2 * i) - 1
        output2[row, col] <- "pval highly sig" 
        Disease_p_values[j, 6] <- "pval highly sig"
      }
      rm(test_pval, a, b, c)
    }
  }
  TableauDiseaseEnrichment <- rbind.data.frame(TableauDiseaseEnrichment, Disease_p_values)
}
output <- cbind(output, output2)

#rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)

length(which(Disease_sample_p_values < 0.05))
length(which(Disease_p_values <= 0.05))

TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]
output <- output[,-1]
filename_Tableau <- paste("Tableau_", filename, sep = "")
write.csv(TableauDiseaseEnrichment, file = filename_Tableau)
write.csv(output, file = filename)

####Combine all Bs for each A w/Disease Enrichment (situation1+6)####
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation1)){
    if(i %% 2 == 0){
      #i = 2
      temp_listofBs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "out")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
      }
      temp_Bs <- unique(unlist(temp_neighbors_of_Bs))
      for(row in 1:length(temp_Bs)){
        temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_Bs[row]
      }
      colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation1)[i], sep = "_")
      temp_df_B_SNs$placeholder <- NA
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation1[which(!is.na(SNs_Situation1[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "in")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  AB_combined_Subnetworks <- Generate.MMandROCombined.SN(temp_df_B_SNs, temp_df_P_SNs, num_neighbors)
  
  lengths_AB <- NULL
  for(i in 1:ncol(AB_combined_Subnetworks)){
    lengths_AB <- append(lengths_AB, length(which(!is.na(AB_combined_Subnetworks[,i]))))
  }
  delete <- max(lengths_AB) + 1:nrow(AB_combined_Subnetworks)
  AB_combined_Subnetworks <- AB_combined_Subnetworks[-delete, ]
  
  subnetworks <- AB_combined_Subnetworks
  
  output <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2))
  
  numSN <- length(colnames(subnetworks))
  numSamples <- 100 #(enter number of sample per group)
  numDiseases <- length(colnames(GeneSetTable))
  
  TableauDiseaseEnrichment <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  colnames(TableauDiseaseEnrichment) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")
  
  Disease_sample_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) * numSamples, ncol = ncol(subnetworks)), row.names = NULL)
  colnames(Disease_sample_p_values) <- colnames(subnetworks)
  
  output2 <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2, ncol = ncol(subnetworks) * 2))
  
  rows <- NULL
  for (i in 1:nrow(output2)){
    var <- ceiling(i/2)
    sub <- ifelse(i %% 2 == 0, "FDR", "pval")
    rows <- append(rows, paste(colnames(GeneSetTable[var]), sub, sep = "_"))
  }
  rownames(output2) <- rows
  rm(rows)
  
  cols <- NULL
  for (i in 1:ncol(output2)){
    var <- ceiling(i/2)
    sub <- ifelse(i %% 2 == 0, "b", "a")
    cols <- append(cols, paste(colnames(subnetworks[var]), sub, sep = "_"))
  }
  colnames(output2) <- cols
  rm(cols)
  
  for (k in 1:numSN) {
    for (j in 1:numDiseases) {
      for (i in 1:numSamples) {
        TargetSize <- length(which(!is.na(subnetworks[, k])))
        if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,k]))) {
          Overlap <- length(which(!is.na(RO_network_genes[match(GeneSetTable[, j], sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL))]))) 
        }
        if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
          Overlap <- length(which(!is.na(GeneSetTable[match(sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
        }
        TestSetSize <- length(which(!is.na(match(RO_network_genes, GeneSetTable[,j]))))
        TotalSize <- length(RO_network_genes)
        
        val <- (j * numSamples) - (numSamples - i)
        
        temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
        
        Disease_sample_p_values[val, k] <- temp_list[[1]]
      }
    }
  }
  
  for (i in 1:numSN) {
    Disease_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) , ncol = 6), row.names = NULL)
    colnames(Disease_p_values) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")
    
    for (j in 1:numDiseases) {
      
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,i]))) {
        Overlap <- length(which(!is.na(subnetworks[match(GeneSetTable[, j], as.character(subnetworks[, i])), i]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,i]))){
        Overlap <- length(which(!is.na(GeneSetTable[match(as.character(subnetworks[, i]), GeneSetTable[, j]), j]))) 
      }
      
      TargetSize <- length(which(!is.na(subnetworks[, i])))
      TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_RO)$name))))
      
      TotalSize <- length(V(g_RO)$name)
      index <- (i * numSN) - (numSN - j) 
      temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
      Disease_p_values[j, 1] <- colnames(GeneSetTable)[j]
      Disease_p_values[j, 2] <- TestSetSize
      Disease_p_values[j, 3] <- paste(unlist(strsplit(colnames(subnetworks)[i], split='_', fixed=TRUE))[1:4], collapse = "_")
      Disease_p_values[j, 4] <- length(na.omit(subnetworks[,i]))
      Disease_p_values[j, 5] <- temp_list[[1]]
      row <- (2 * j) - 1
      col <- (2 * i) - 1
      output2[row, col] <- temp_list[[1]]
      row2 <- (2 * j) - 1
      col2 <- 2 * i
      output2[row2, col2] <- Overlap
      row3 <- 2 * j
      col3 <- 2 * i
      output2[row3, col3] <- TestSetSize
      
      #FDR = a/(b*c)
      #where: a = # of actual pvalues for that KD that are less than the actual p-values in that KD; b = # of sample p-values for that KD that are less than the actual p-value in that KD * number of samples generated per pathway per KD
      test_pval <- Disease_p_values[j, 5]
      
      if (test_pval <= 0.05){
        a <- length(which(Disease_p_values[, 5] <= test_pval))
        b <- length(which(Disease_sample_p_values[, i] <= test_pval))
        c <- numSamples
        if (b != 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- b / (a * c)
          Disease_p_values[j, 6] <- b / (a * c)
        }
        if (b == 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- "pval highly sig" 
          Disease_p_values[j, 6] <- "pval highly sig"
        }
        rm(test_pval, a, b, c)
      }
    }
    TableauDiseaseEnrichment <- rbind.data.frame(TableauDiseaseEnrichment, Disease_p_values)
  }
  output <- cbind(output, output2)
  
  #rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)
  
  length(which(Disease_sample_p_values < 0.05))
  length(which(Disease_p_values <= 0.05))
  
  TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]
  output <- output[,-1]
  
  filename <- paste("Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  filename <- paste( "DiseaseEnrichment_Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(output, file = filename)
  
  filename_Tableau <- paste("Tableau_", filename, sep = "")
  write.csv(TableauDiseaseEnrichment, file = filename_Tableau)
  
  assign(paste("Situation1_PAB_combined_SN_d", num_neighbors, sep = "_"), AB_combined_Subnetworks)
}

####Combine all Bs for each A w/Disease Enrichment (situation2+5)####
neighbors_list <- c("C", "D", "E", "F", "G", "P", "Q", "R", "S", "T")
for (num_neighbors in 1:5){
  temp_df_B_SNs <- data.frame(matrix(NA, nrow = length(V(g_RO)$name), ncol = 1))
  temp_df_P_SNs <- data.frame(matrix(NA, nrow = length(V(g_MM)$name), ncol = 1))
  
  for (i in 1:ncol(SNs_Situation2)){
    if(i %% 2 == 0){
      #i = 2
      temp_listofBs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Bs <- match(temp_listofBs, V(g_RO)$name)
      
      temp_neighbors_of_Bs <- neighborhood(g_RO, nodes = temp_pos_Bs, order = num_neighbors, mode = "in")
      for(list_index in 1:length(temp_neighbors_of_Bs)){
        temp_neighbors_of_Bs[[list_index]] <- V(g_RO)$name[temp_neighbors_of_Bs[[list_index]]]
      }
      temp_Bs <- unique(unlist(temp_neighbors_of_Bs))
      for(row in 1:length(temp_Bs)){
        temp_df_B_SNs[row, ncol(temp_df_B_SNs)] <- temp_Bs[row]
      }
      colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation2)[i], sep = "_")
      temp_df_B_SNs$placeholder <- NA
    } else{
      #    i = 3
      temp_listofPs <- SNs_Situation2[which(!is.na(SNs_Situation2[,i])), i]
      temp_pos_Ps <- match(temp_listofPs, V(g_MM)$name)
      
      temp_neighbors_of_Ps <- neighborhood(g_MM, nodes = temp_pos_Ps, order = num_neighbors-1, mode = "out")
      
      for (list_index in 1:length(temp_neighbors_of_Ps)){
        temp_neighbors_of_Ps[[list_index]] <- V(g_MM)$name[temp_neighbors_of_Ps[[list_index]]]
      }
      temp_neighbors_of_Ps <- unique(unlist(temp_neighbors_of_Ps))
      
      length_P_neighbors <- length(temp_neighbors_of_Ps)
      
      for(row in 1:length(temp_neighbors_of_Ps)){
        temp_df_P_SNs[row, ncol(temp_df_P_SNs)] <- temp_neighbors_of_Ps[row]
      }
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation2)[i], temp_neighbors_of_Ps[2], sep = "_")
      temp_df_P_SNs$placeholder <- NA
    }
  }
  temp_df_B_SNs$placeholder <- NULL
  temp_df_P_SNs$placeholder <- NULL
  lengths_B <- NULL
  
  for(i in 1:ncol(temp_df_B_SNs)){
    lengths_B <- append(lengths_B, length(which(!is.na(temp_df_B_SNs[,i]))))
  }
  delete <- (max(lengths_B) + 1):nrow(temp_df_B_SNs)
  temp_df_B_SNs <- temp_df_B_SNs[-delete, ]
  
  lengths_P <- NULL
  for(i in 1:ncol(temp_df_P_SNs)){
    lengths_P <- append(lengths_P, length(which(!is.na(temp_df_P_SNs[,i]))))
  }
  delete <- (max(lengths_P) + 1):nrow(temp_df_P_SNs)
  temp_df_P_SNs <- temp_df_P_SNs[-delete, ]
  
  AB_combined_Subnetworks <- Generate.MMandROCombined.SN(temp_df_B_SNs, temp_df_P_SNs, num_neighbors)
  
  lengths_AB <- NULL
  for(i in 1:ncol(AB_combined_Subnetworks)){
    lengths_AB <- append(lengths_AB, length(which(!is.na(AB_combined_Subnetworks[,i]))))
  }
  delete <- max(lengths_AB) + 1:nrow(AB_combined_Subnetworks)
  AB_combined_Subnetworks <- AB_combined_Subnetworks[-delete, ]
  
  subnetworks <- AB_combined_Subnetworks
  
  output <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2))
  
  numSN <- length(colnames(subnetworks))
  numSamples <- 100 #(enter number of sample per group)
  numDiseases <- length(colnames(GeneSetTable))
  
  TableauDiseaseEnrichment <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  colnames(TableauDiseaseEnrichment) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")
  
  Disease_sample_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) * numSamples, ncol = ncol(subnetworks)), row.names = NULL)
  colnames(Disease_sample_p_values) <- colnames(subnetworks)
  
  output2 <- data.frame(matrix(NA, nrow = ncol(GeneSetTable)*2, ncol = ncol(subnetworks) * 2))
  
  rows <- NULL
  for (i in 1:nrow(output2)){
    var <- ceiling(i/2)
    sub <- ifelse(i %% 2 == 0, "FDR", "pval")
    rows <- append(rows, paste(colnames(GeneSetTable[var]), sub, sep = "_"))
  }
  rownames(output2) <- rows
  rm(rows)
  
  cols <- NULL
  for (i in 1:ncol(output2)){
    var <- ceiling(i/2)
    sub <- ifelse(i %% 2 == 0, "b", "a")
    cols <- append(cols, paste(colnames(subnetworks[var]), sub, sep = "_"))
  }
  colnames(output2) <- cols
  rm(cols)
  
  for (k in 1:numSN) {
    for (j in 1:numDiseases) {
      for (i in 1:numSamples) {
        TargetSize <- length(which(!is.na(subnetworks[, k])))
        if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,k]))) {
          Overlap <- length(which(!is.na(RO_network_genes[match(GeneSetTable[, j], sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL))]))) 
        }
        if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
          Overlap <- length(which(!is.na(GeneSetTable[match(sample(RO_network_genes, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
        }
        TestSetSize <- length(which(!is.na(match(RO_network_genes, GeneSetTable[,j]))))
        TotalSize <- length(RO_network_genes)
        
        val <- (j * numSamples) - (numSamples - i)
        
        temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
        
        Disease_sample_p_values[val, k] <- temp_list[[1]]
      }
    }
  }
  
  for (i in 1:numSN) {
    Disease_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) , ncol = 6), row.names = NULL)
    colnames(Disease_p_values) <- c("Disease", "Disease_num_genes", "SN", "SN_num_genes", "p_value", "FDR")
    
    for (j in 1:numDiseases) {
      
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,i]))) {
        Overlap <- length(which(!is.na(subnetworks[match(GeneSetTable[, j], as.character(subnetworks[, i])), i]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,i]))){
        Overlap <- length(which(!is.na(GeneSetTable[match(as.character(subnetworks[, i]), GeneSetTable[, j]), j]))) 
      }
      
      TargetSize <- length(which(!is.na(subnetworks[, i])))
      TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_RO)$name))))
      
      TotalSize <- length(V(g_RO)$name)
      index <- (i * numSN) - (numSN - j) 
      temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
      Disease_p_values[j, 1] <- colnames(GeneSetTable)[j]
      Disease_p_values[j, 2] <- TestSetSize
      Disease_p_values[j, 3] <- paste(unlist(strsplit(colnames(subnetworks)[i], split='_', fixed=TRUE))[1:4], collapse = "_")
      Disease_p_values[j, 4] <- length(na.omit(subnetworks[,i]))
      Disease_p_values[j, 5] <- temp_list[[1]]
      row <- (2 * j) - 1
      col <- (2 * i) - 1
      output2[row, col] <- temp_list[[1]]
      row2 <- (2 * j) - 1
      col2 <- 2 * i
      output2[row2, col2] <- Overlap
      row3 <- 2 * j
      col3 <- 2 * i
      output2[row3, col3] <- TestSetSize
      
      #FDR = a/(b*c)
      #where: a = # of actual pvalues for that KD that are less than the actual p-values in that KD; b = # of sample p-values for that KD that are less than the actual p-value in that KD * number of samples generated per pathway per KD
      test_pval <- Disease_p_values[j, 5]
      
      if (test_pval <= 0.05){
        a <- length(which(Disease_p_values[, 5] <= test_pval))
        b <- length(which(Disease_sample_p_values[, i] <= test_pval))
        c <- numSamples
        if (b != 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- b / (a * c)
          Disease_p_values[j, 6] <- b / (a * c)
        }
        if (b == 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- "pval highly sig" 
          Disease_p_values[j, 6] <- "pval highly sig"
        }
        rm(test_pval, a, b, c)
      }
    }
    TableauDiseaseEnrichment <- rbind.data.frame(TableauDiseaseEnrichment, Disease_p_values)
  }
  output <- cbind(output, output2)
  
  #rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)
  
  length(which(Disease_sample_p_values < 0.05))
  length(which(Disease_p_values <= 0.05))
  
  TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]
  output <- output[,-1]
  
  
  filename <- paste("Situation2_BAP_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  filename <- paste( "DiseaseEnrichment_Situation2_BAP_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(output, file = filename)
  
  filename_Tableau <- paste("Tableau_", filename, sep = "")
  write.csv(TableauDiseaseEnrichment, file = filename_Tableau)
  
  assign(paste("Situation2_BAP_combined_SN_d", num_neighbors, sep = "_"), AB_combined_Subnetworks)
}


