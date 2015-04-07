####Combine all Bs for each A w/Disease Enrichment (situation1+6)####
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
  
 # rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)
  
  length(which(Disease_sample_p_values < 0.05))
  length(which(Disease_p_values <= 0.05))
  
  TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]
  output <- output[,-1]
  
  assign("test_PW_DIS", output)
  
  temp_df_sig <- Pairwise.Disease.Table.Sig(test_PW_DIS)
  temp_df_nom_sig <- Pairwise.Disease.Table.Nom.Sig(test_PW_DIS)
  
  filename_pairwise_sig <- paste("Dis_Sig_Pairwise_Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(temp_df_sig, file = filename_pairwise_sig)
  
  filename_pairwise_nom_sig <- paste("Dis_Nom_Sig_Pairwise_Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(temp_df_nom_sig, file = filename_pairwise_nom_sig)
  
  filename <- paste("Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(AB_combined_Subnetworks, file = filename)
  
  filename <- paste( "DiseaseEnrichment_Situation1_PAB_combined_SN_d_", num_neighbors, ".csv", sep = "")
  write.csv(output, file = filename)
  
  filename_Tableau <- paste("Tableau_", filename, sep = "")
  write.csv(TableauDiseaseEnrichment, file = filename_Tableau)
  
  assign(paste("Situation1_PAB_combined_SN_d", num_neighbors, sep = "_"), AB_combined_Subnetworks)
}

####Create Pairwise table for Significant Enrichment Scores####
Pairwise.Disease.Table.Sig <- function(test_PW_DIS) {
  for(i in 1:nrow(test_PW_DIS)){
    for (j in 1:ncol(test_PW_DIS)){

      if((i %% 2 == 0) && (j %% 2 != 0)){ #Puts NA for any FDR values that are > 0.1
        test_PW_DIS[i, j] <- ifelse(test_PW_DIS[i, j] > 0.100000001, ifelse(test_PW_DIS[i, j] == "pval highly sig", test_PW_DIS[i, j], "NA"), test_PW_DIS[i, j])
        test_PW_DIS[i-1, j] <- ifelse(test_PW_DIS[i,j] == "NA", NA, test_PW_DIS[i-1, j]) #Comment out for nominally significant
      }
      
      if(i %% 2 != 0){ #Puts NA for any pval that is > 0.05
        test_PW_DIS[i, j] <- ifelse(as.numeric(test_PW_DIS[i, j]) > 0.05, NA, test_PW_DIS[i, j])
      }
      
      if(j %% 2 == 0){ #Puts NA for all counts
        test_PW_DIS[i, j] <- NA
      }
    }
  }
  
  PVAL_df <- data.frame(matrix(NA, nrow = ncol(GeneSetTable), ncol = ncol(test_PW_DIS)/2), row.names = colnames(GeneSetTable))
  
  for(i in 1:ncol(test_PW_DIS)){
    for (j in 1:nrow(test_PW_DIS)){
      if(i %% 2 != 0 && j %% 2 != 0){
        temp_row <- ceiling(j/2)
        temp_col <- ceiling(i/2)
        PVAL_df[temp_row, temp_col] <- test_PW_DIS[j, i]
        colnames(PVAL_df)[temp_col] <- colnames(test_PW_DIS)[i]
      }
    }
  }
  
  PVAL_df <- ifelse(is.na(PVAL_df), 0, 1)
  
  test <- as.data.frame(PVAL_df %*% t(PVAL_df))
  
  return(test)
}

####Create Pairwise table for Nominally Significant Enrichment Scores####
Pairwise.Disease.Table.Nom.Sig <- function(test_PW_DIS) {
  for(i in 1:nrow(test_PW_DIS)){
    for (j in 1:ncol(test_PW_DIS)){
      
      if((i %% 2 == 0) && (j %% 2 != 0)){ #Puts NA for any FDR values that are > 0.1
        test_PW_DIS[i, j] <- ifelse(test_PW_DIS[i, j] > 0.100000001, ifelse(test_PW_DIS[i, j] == "pval highly sig", test_PW_DIS[i, j], "NA"), test_PW_DIS[i, j])
      }
      
      if(i %% 2 != 0){ #Puts NA for any pval that is > 0.05
        test_PW_DIS[i, j] <- ifelse(as.numeric(test_PW_DIS[i, j]) > 0.05, NA, test_PW_DIS[i, j])
      }
      
      if(j %% 2 == 0){ #Puts NA for all counts
        test_PW_DIS[i, j] <- NA
      }
      
    }
  }
  
  PVAL_df <- data.frame(matrix(NA, nrow = ncol(GeneSetTable), ncol = ncol(test_PW_DIS)/2), row.names = colnames(GeneSetTable))
  
  for(i in 1:ncol(test_PW_DIS)){
    for (j in 1:nrow(test_PW_DIS)){
      if(i %% 2 != 0 && j %% 2 != 0){
        temp_row <- ceiling(j/2)
        temp_col <- ceiling(i/2)
        PVAL_df[temp_row, temp_col] <- test_PW_DIS[j, i]
        colnames(PVAL_df)[temp_col] <- colnames(test_PW_DIS)[i]
      }
    }
  }
  
  PVAL_df <- ifelse(is.na(PVAL_df), 0, 1)
  
  test <- as.data.frame(PVAL_df %*% t(PVAL_df))
  
  return(test)
}