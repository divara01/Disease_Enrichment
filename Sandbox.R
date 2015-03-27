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

####Generate SN combining P and each B####
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
        colnames(temp_df_B_SNs)[ncol(temp_df_B_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Bs[[list_index]][1], "g", length(which(!is.na(temp_df_B_SNs[, ncol(temp_df_B_SNs)]))), sep = "_")
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
      colnames(temp_df_P_SNs)[ncol(temp_df_P_SNs)] <- paste("A", colnames(SNs_Situation1)[i], temp_neighbors_of_Ps[2], "g", length(which(!is.na(temp_df_P_SNs[, ncol(temp_df_P_SNs)]))), sep = "_")
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
  assign(paste("Situation1_PAB_combined_SN_d_", num_neighbors, sep = ""), AB_combined_Subnetworks)
}

####Combine all Bs for each A (situation1+6)####
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
  
  Disease_sample_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) * numSamples, ncol = ncol(subnetworks)), row.names = NULL)
  colnames(Disease_sample_p_values) <- colnames(subnetworks)
  
  Disease_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) , ncol = numSN), row.names = colnames(GeneSetTable))
  colnames(Disease_p_values) <- colnames(subnetworks)
  
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
    cols <- append(cols, paste(colnames(Disease_p_values[var]), sub, sep = "_"))
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
    for (j in 1:numDiseases) {
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,i]))) {
        Overlap <- length(which(!is.na(subnetworks[match(GeneSetTable[, j], as.character(subnetworks[, i])), i]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,i]))){
        Overlap <- length(which(!is.na(GeneSetTable[match(as.character(subnetworks[, i]), GeneSetTable[, j]), j]))) 
      }
      
      TargetSize <- length(which(!is.na(subnetworks[, i])))
      TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], RO_network_genes))))
      
      TotalSize <- length(RO_network_genes)
      index <- (i * numSN) - (numSN - j) 
      temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
      Disease_p_values[j, i] <- temp_list[[1]]
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
      
      test_pval <- Disease_p_values[j, i]
      
      if (test_pval <= 0.05){
        a <- length(which(Disease_p_values[, i] <= test_pval))
        b <- length(which(Disease_sample_p_values[, i] <= test_pval))
        c <- numSamples
        if (b != 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- b / (a * c)
        }
        if (b == 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- "pval highly sig"      
        }
        rm(test_pval, a, b, c)
      }
    }
  }
  output <- cbind(output, output2)
  
  rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)
  
  length(which(Disease_sample_p_values < 0.05))
  length(which(Disease_p_values <= 0.05))
  
  output <- output[,-1]  
  
  filename <- paste("Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  filename <- paste( "DiseaseEnrichment_Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(output, file = filename)
  
  assign(paste("Situation1_PAB_combined_SN_d", num_neighbors, sep = "_"), AB_combined_Subnetworks)
}

####Combine all Bs for each A (situation2+5)####
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
  
  Disease_sample_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) * numSamples, ncol = ncol(subnetworks)), row.names = NULL)
  colnames(Disease_sample_p_values) <- colnames(subnetworks)
  
  Disease_p_values <- data.frame(matrix(NA, nrow = ncol(GeneSetTable) , ncol = numSN), row.names = colnames(GeneSetTable))
  colnames(Disease_p_values) <- colnames(subnetworks)
  
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
    cols <- append(cols, paste(colnames(Disease_p_values[var]), sub, sep = "_"))
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
    for (j in 1:numDiseases) {
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,i]))) {
        Overlap <- length(which(!is.na(subnetworks[match(GeneSetTable[, j], as.character(subnetworks[, i])), i]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) > length(na.omit(subnetworks[,i]))){
        Overlap <- length(which(!is.na(GeneSetTable[match(as.character(subnetworks[, i]), GeneSetTable[, j]), j]))) 
      }
      
      TargetSize <- length(which(!is.na(subnetworks[, i])))
      TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], RO_network_genes))))
      
      TotalSize <- length(RO_network_genes)
      index <- (i * numSN) - (numSN - j) 
      temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
      Disease_p_values[j, i] <- temp_list[[1]]
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
      
      test_pval <- Disease_p_values[j, i]
      
      if (test_pval <= 0.05){
        a <- length(which(Disease_p_values[, i] <= test_pval))
        b <- length(which(Disease_sample_p_values[, i] <= test_pval))
        c <- numSamples
        if (b != 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- b / (a * c)
        }
        if (b == 0){
          row = 2 * j
          col = (2 * i) - 1
          output2[row, col] <- "pval highly sig"      
        }
        rm(test_pval, a, b, c)
      }
    }
  }
  output <- cbind(output, output2)
  
  rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)
  
  length(which(Disease_sample_p_values < 0.05))
  length(which(Disease_p_values <= 0.05))
  
  output <- output[,-1]  
  
  filename <- paste("Situation2_BAP_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  filename <- paste( "DiseaseEnrichment_Situation2_BAP_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(output, file = filename)
  
  assign(paste("Situation2_BAP_combined_SN_d", num_neighbors, sep = "_"), AB_combined_Subnetworks)
}
