####Identify MM -> RO SN connections (situation 1 + 6)####

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

####Identify Each B and P for each A (situation 1)####
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

####Generate Subnetworks for each B and P to each degree (situation 1)####
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

####Try to generate SN combining P and each B-WORKS!####
AB_combined_Subnetworks <- data.frame(matrix(NA, nrow = length(V(g_whole_omental)$name), ncol = ncol(Situation1_B_C_SNs)))
colnames(AB_combined_Subnetworks) <- colnames(Situation1_B_C_SNs)

A_list_from_P <- NULL
for(A_index in 1:ncol(Situation1_A_P_SNs)){
  A_list_from_P <- append(A_list_from_P, unlist(strsplit(colnames(Situation1_A_P_SNs)[A_index], split='_', fixed=TRUE))[2])
}

for(B_index in 1:ncol(AB_combined_Subnetworks)){
  Individual_B_from_A <- unlist(strsplit(colnames(Situation1_B_C_SNs)[B_index], split='_', fixed=TRUE))[2]
  P_SN <- as.vector(na.omit(Situation1_A_P_SNs[, match(Individual_B_from_A, A_list_from_P)]))
  B_SN <- as.vector(na.omit(Situation1_B_C_SNs[,B_index]))
  PAB_SN <- union(P_SN, B_SN)
  
  for(input_index in 1:length(PAB_SN)){
    AB_combined_Subnetworks[input_index, B_index] <- PAB_SN[input_index]
  }
}

lengths_AB <- NULL
for(i in 1:ncol(AB_combined_Subnetworks)){
  lengths_AB <- append(lengths_AB, length(which(!is.na(AB_combined_Subnetworks[,i]))))
}
delete <- (max(lengths_AB) + 1):nrow(AB_combined_Subnetworks)
AB_combined_Subnetworks <- AB_combined_Subnetworks[-delete, ]

####Test Function####
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
  
  AB_combined_Subnetworks <- Generate.PAB.SN.Sit1(temp_df_B_SNs, temp_df_P_SNs, num_neighbors)

  lengths_AB <- NULL
  for(i in 1:ncol(AB_combined_Subnetworks)){
    lengths_AB <- append(lengths_AB, length(which(!is.na(AB_combined_Subnetworks[,i]))))
  }
  delete <- max(lengths_AB) + 1:nrow(AB_combined_Subnetworks)
  AB_combined_Subnetworks <- AB_combined_Subnetworks[-delete, ]
  
  filename <- paste("Situation1_AB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  assign(paste("Situation1_B", neighbors_list[num_neighbors], "SNs", sep = "_"), temp_df_B_SNs)
  assign(paste("Situation1_A", neighbors_list[num_neighbors+5], "SNs", sep = "_"), temp_df_P_SNs)
  assign(paste("Situation1_AB_combined_SN_d_", num_neighbors, sep = ""), AB_combined_Subnetworks)
}


####FUNCTIONS####
Generate.PAB.SN.Sit1 <- function(temp_df_B_SNs, temp_df_P_SNs, num_neighbors) {
  temp_combined_Subnetworks <- data.frame(matrix(NA, nrow = length(V(g_whole_omental)$name), ncol = ncol(temp_df_B_SNs)))
  colnames(temp_combined_Subnetworks) <- colnames(temp_df_B_SNs)
  
  A_list_from_P <- NULL
  for(A_index in 1:ncol(temp_df_P_SNs)){
    A_list_from_P <- append(A_list_from_P, unlist(strsplit(colnames(temp_df_P_SNs)[A_index], split='_', fixed=TRUE))[2])
  }
  
  for(B_index in 1:ncol(temp_combined_Subnetworks)){
    #B_index = 250
    Individual_B_from_A <- unlist(strsplit(colnames(temp_df_B_SNs)[B_index], split='_', fixed=TRUE))[2]
    P_SN <- as.vector(na.omit(temp_df_P_SNs[, match(Individual_B_from_A, A_list_from_P)]))
    B_SN <- as.vector(na.omit(temp_df_B_SNs[,B_index]))
    PAB_SN <- union(P_SN, B_SN)
    
    for(input_index in 1:length(PAB_SN)){
      AB_combined_Subnetworks[input_index, B_index] <- PAB_SN[input_index]
    }
    rm(P_SN, B_SN, Individual_B_from_A, PAB_SN)
  }
  
  return(AB_combined_Subnetworks)
}

