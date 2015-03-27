---
  title: "Functions"
author: "A. Divaraniya"
date: "October 8, 2014"
output: html_document
---
#************************************************************************************
####Stouffer's Method####
#Use - To calculate p-values for each gene in AD across various tissues 
Stouffer.test <- function(p, w) {
  if (missing(w)) {
    w <- rep(1, length(p)) / length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- (qnorm(1-p))* direction #Incorporate sign here (multiplied) USE FUNCTION SIGN
  Z <- sum( w * Zi) / sqrt(sum(w ^ 2))
  directionZ <- sign(Z)
  p.val <- 1 - pnorm(Z)
  return(c(Z = Z, p.value = p.val, direction = directionZ))
}
#************************************************************************************

#************************************************************************************

####Calculate Over Enrichment Score Function####
#Use - To calculatate enrichment score and return over enrichment p-value from Fisher Test
OverEnrichment <- function(Overlap, TargetSize, TestSetSize, TotalSize) {
  
  a <- Overlap #number matching genes in disease signature to each subnetwork
  b <- TestSetSize - a #number genes in disease signature - number matching genes in disease signature to each subnetwork
  c <- TargetSize - a #number genes in each subnetwork - number matching genes in disease signature to each subnetwork
  d <- TotalSize - TestSetSize - c #total number genes in network - number genes in disease signature - c
  
  ContinTable <- matrix(data =c(a,b,c,d),nrow = 2, ncol =2, byrow = T, dimnames = NULL)
  
  FisherTestOverEnrich <- fisher.test(ContinTable, alternative = "greater")
  FisherTestUnderEnrich <- fisher.test(ContinTable, alternative = "less")
  FisherResultsOverEnrich <- as.data.frame(t(cbind(FisherTestOverEnrich)))
  FisherResultsUnderEnrich <- as.data.frame(t(cbind(FisherTestUnderEnrich)))
  #rm(a, b, c, d)
  
  return(c(FisherTestOverEnrich$p.value, FisherTestOverEnrich$estimate))
}
#************************************************************************************

#************************************************************************************

#####GenerateOmentalSubnetworks_InfoFlowOmentalMM Function####
#Use - Determine Direction of Information Flow From and To Omental From and To Magic Module (use this for permutations as well)
GenerateOmentalSubnetworks_InfoFlowOmentalMM <- function(num_neighbors, mode_magic, mode_neighbors, magic_module){
  
  RO_network_genes <- V(g_whole_omental)$name[which(!(V(g_whole_omental)$name %in% magic_module))]
  
  RO_network_pos <- match(RO_network_genes, V(g_whole_omental)$name)
  
  g_RO_func <- induced.subgraph(g_whole_omental, RO_network_pos, impl = "auto")
  
  #1. Identify first degree neighbors of magic module genes
  core_magic_pos_on_omental <- match(magic_module, V(g_whole_omental)$name)
  core_magic_pos_on_omental <- core_magic_pos_on_omental[!is.na(core_magic_pos_on_omental)]
  
  magic_neighbors_list <- neighborhood(g_whole_omental, order = 1, nodes = core_magic_pos_on_omental, mode = mode_magic)
  
  #2. Convert positions into Entrez Gene ID
  for (list_index in 1: length(magic_neighbors_list)){
    magic_neighbors_list[[list_index]] <- V(g_whole_omental)$name[magic_neighbors_list[[list_index]]]
  }
  rm(list_index)
  
  #3. Unlist magic_neighbors_list
  unlist_magic_neighbors <- unique(unlist(magic_neighbors_list))
  
  #4. Remove Magic Module genes from RO neighbors list
  RO_neighbors <- unlist_magic_neighbors[which(!(unlist_magic_neighbors %in% magic_module))] #369 This is my list of Bs
  
  #5a. Get neighbors of these RO genes (this will be your subnetwork for each B out to num_neighbors)
  RO_pos_on_omental <- match(RO_neighbors, V(g_RO_func)$name)
  
  RO_neighbors_list <- neighborhood(g_RO_func, order = num_neighbors, nodes = RO_pos_on_omental, mode = mode_neighbors)
  
  #6a. Convert RO_neighbors_list to Entrez ID
  for (list_index in 1:length(RO_neighbors_list)){
    RO_neighbors_list[[list_index]] <- V(g_RO_func)$name[RO_neighbors_list[[list_index]]]
  }
  rm(list_index)
  
  for(i in 1:length(RO_neighbors_list)){
    RO_neighbors_list[[i]] <- RO_neighbors_list[[i]][which(is.na(match(RO_neighbors_list[[i]], magic_module)))]
  }
  rm(i)
  
  #Get first degree neighbor of each gene within each B subnetwork
  temp_df <- data.frame(matrix(NA, nrow = 1, ncol = 9))
  colnames(temp_df) <- c("num_genes_in_SN", "MM_from_SN_OSN",  "RO_from_SN_OSN",  "MM_from_SN/RO_OSN",  "RO_from_SN_ISN", "MM_to_SN_OSN", "RO_to_SN_OSN",  "MM_to_SN/RO_OSN", "RO_to_SN_ISN")
  
  for (i in 1:length(RO_neighbors_list)){
    x <- RO_neighbors_list[[i]]
    x <- x[-1]
    
    pos <- match(x, V(g_whole_omental)$name)
    num_genes_SN <- length(pos)  
    if (num_genes_SN == 0){
      num_magic_to_OSN = 0
      num_RO_to_OSN = 0
      ratio_to_OSN = 0
      num_magic_from_OSN = 0
      num_RO_from_OSN = 0
      ratio_from_OSN = 0
      num_RO_from_ISN = 0
      num_RO_to_ISN = 0
    }
    
    if(num_genes_SN != 0){
      neighbors_x_list_from <- neighborhood(g_whole_omental, order = 1, nodes = pos, mode = "out")
      neighbors_x_list_to <- neighborhood(g_whole_omental, order = 1, nodes = pos, mode = "in")
      ####      
      #FROM#
      ###
      for (list_index in 1: length(neighbors_x_list_from)){
        neighbors_x_list_from[[list_index]] <- V(g_whole_omental)$name[neighbors_x_list_from[[list_index]]]
      }
      
      unlist_neighbors_x_list_from <- unlist(neighbors_x_list_from)
      unlist_neighbors_x_list_from_OSN <- unlist_neighbors_x_list_from[!(unlist_neighbors_x_list_from %in% RO_neighbors_list[[i]])] 
      unlist_neighbors_x_list_from_ISN <- unlist_neighbors_x_list_from[(unlist_neighbors_x_list_from %in% RO_neighbors_list[[i]])] 
      
      num_magic_from_OSN <-  length(which(unlist_neighbors_x_list_from_OSN %in% magic_module))
      num_RO_from_OSN <- length(which(!(unlist_neighbors_x_list_from_OSN %in% magic_module))) 
      ratio_from_OSN = ifelse(num_RO_from_OSN == 0, ifelse(num_magic_from_OSN == 0, NaN, Inf), round(num_magic_from_OSN/num_RO_from_OSN, 4))
      
      num_RO_from_ISN <- length(which(!(unlist_neighbors_x_list_from_ISN %in% magic_module))) 
      
      ###
      #TO#
      ###
      for (list_index in 1: length(neighbors_x_list_to)){
        neighbors_x_list_to[[list_index]] <- V(g_whole_omental)$name[neighbors_x_list_to[[list_index]]]
      }
      
      unlist_neighbors_x_list_to <- unlist(neighbors_x_list_to)
      unlist_neighbors_x_list_to_OSN <- unlist_neighbors_x_list_to[!(unlist_neighbors_x_list_to %in% RO_neighbors_list[[i]])] 
      unlist_neighbors_x_list_to_ISN <- unlist_neighbors_x_list_to[(unlist_neighbors_x_list_to %in% RO_neighbors_list[[i]])]       
      
      num_magic_to_OSN <-  length(which(unlist_neighbors_x_list_to_OSN %in% magic_module))
      num_RO_to_OSN <- length(which(!(unlist_neighbors_x_list_to_OSN %in% magic_module))) 
      ratio_to_OSN = ifelse(num_RO_to_OSN == 0, ifelse(num_magic_to_OSN == 0, NaN, Inf), round(num_magic_to_OSN/num_RO_to_OSN, 4))
      
      num_RO_to_ISN <- length(which(!(unlist_neighbors_x_list_to_ISN %in% magic_module))) 
    }
    
    vect <- c(num_genes_SN, num_magic_from_OSN, num_RO_from_OSN, ratio_from_OSN, num_RO_from_ISN, num_magic_to_OSN, num_RO_to_OSN, ratio_to_OSN, num_RO_to_ISN)
    rm(num_genes_SN, num_magic_from_OSN, num_RO_from_OSN, ratio_from_OSN, num_RO_from_ISN, num_magic_to_OSN, num_RO_to_OSN, ratio_to_OSN, num_RO_to_ISN)
    temp_df <- rbind.data.frame(temp_df, vect)
    
    row.names(temp_df)[nrow(temp_df)] <- paste(RO_neighbors_list[[i]][1], "d", num_neighbors, sep= "_")
  }
  
  temp_df <- temp_df[-1,]

  
  return(list(temp_df, RO_neighbors_list, RO_neighbors))
}
#************************************************************************************

#************************************************************************************

####GenerateMagicSubnetworks_InfoFlowMMOmental Function####
#Use - Situations #5-8: Determine Direction of Information Flow From and To Omental From and To Magic Module
GenerateMagicSubnetworks_InfoFlowMMOmental <- function(num_neighbors, mode_RO, mode_neighbors){
  #1. Identify first degree neighbors of RO genes
  RO_pos_on_omental <- match(RO_network_genes, V(g_whole_omental)$name)
  RO_pos_on_omental <- RO_pos_on_omental[!is.na(RO_pos_on_omental)]
  
  RO_neighbors_list <- neighborhood(g_whole_omental, order = 1, nodes = RO_pos_on_omental, mode = mode_RO)
  
  #2. Convert positions into Entrez Gene ID
  for (list_index in 1: length(RO_neighbors_list)){
    RO_neighbors_list[[list_index]] <- V(g_whole_omental)$name[RO_neighbors_list[[list_index]]]
  }
  rm(list_index)
  
  #3. Unlist RO_neighbors_list
  unlist_RO_neighbors <- unique(unlist(RO_neighbors_list)) #5070 unique genes
  
  #4. Remove RO genes from magic neighbors list
  MM_neighbors <- unlist_RO_neighbors[which((unlist_RO_neighbors %in% magic_module_genes$EntrezID))] #223 This is my list of As
  
  #5a. Get neighbors of these magic genes (this will be your subnetwork for each A out to num_neighbors)
  MM_pos_on_omental <- match(MM_neighbors, V(g_MM)$name)
  
  MM_neighbors_list <- neighborhood(g_MM, order = num_neighbors, nodes = MM_pos_on_omental, mode = mode_neighbors)
  
  #6a. Convert MM_neighbors_list to Entrez ID
  for (list_index in 1:length(MM_neighbors_list)){
    MM_neighbors_list[[list_index]] <- V(g_MM)$name[MM_neighbors_list[[list_index]]]
  }
  rm(list_index)
  
  #Get first degree neighbor of each gene within each Bm subnetwork
  temp_df <- data.frame(matrix(NA, nrow = 1, ncol = 9))
  colnames(temp_df) <- c("num_genes_in_SN", "MM_from_SN_OSN",  "RO_from_SN_OSN",  "MM_from_SN/RO_OSN",  "RO_from_SN_ISN", "MM_to_SN_OSN", "RO_to_SN_OSN",  "MM_to_SN/RO_OSN", "RO_to_SN_ISN")
  
  for (i in 1:length(MM_neighbors_list)){
    x <- MM_neighbors_list[[i]]
    x <- x[-1]
    
    pos <- match(x, V(g_whole_omental)$name)
    num_genes_SN <- length(pos)  
    if (num_genes_SN == 0){
      num_magic_to_OSN = 0
      num_RO_to_OSN = 0
      ratio_to_OSN = 0
      num_magic_from_OSN = 0
      num_RO_from_OSN = 0
      ratio_from_OSN = 0
      num_MM_from_ISN = 0
      num_MM_to_ISN = 0
    }
    
    if(num_genes_SN != 0){
      neighbors_x_list_from <- neighborhood(g_whole_omental, order = 1, nodes = pos, mode = "out")
      neighbors_x_list_to <- neighborhood(g_whole_omental, order = 1, nodes = pos, mode = "in")
      ###    
      #FROM#
      ###
      for (list_index in 1: length(neighbors_x_list_from)){
        neighbors_x_list_from[[list_index]] <- V(g_whole_omental)$name[neighbors_x_list_from[[list_index]]]
      }
      
      unlist_neighbors_x_list_from <- unlist(neighbors_x_list_from)
      unlist_neighbors_x_list_from_OSN <- unlist_neighbors_x_list_from[!(unlist_neighbors_x_list_from %in% MM_neighbors_list[[i]])] 
      unlist_neighbors_x_list_from_ISN <- unlist_neighbors_x_list_from[(unlist_neighbors_x_list_from %in% MM_neighbors_list[[i]])] 
      
      num_magic_from_OSN <-  length(which(unlist_neighbors_x_list_from_OSN %in% magic_module_genes$EntrezID))
      num_RO_from_OSN <- length(which(!(unlist_neighbors_x_list_from_OSN %in% magic_module_genes$EntrezID))) 
      ratio_from_OSN = ifelse(num_RO_from_OSN == 0, ifelse(num_magic_from_OSN == 0, NaN, Inf), round(num_magic_from_OSN/num_RO_from_OSN, 4))
      
      num_MM_from_ISN <- length(unlist_neighbors_x_list_from_ISN %in% magic_module_genes$EntrezID) 
      
      ###
      #TO#
      ###
      for (list_index in 1: length(neighbors_x_list_to)){
        neighbors_x_list_to[[list_index]] <- V(g_whole_omental)$name[neighbors_x_list_to[[list_index]]]
      }
      
      unlist_neighbors_x_list_to <- unlist(neighbors_x_list_to)
      unlist_neighbors_x_list_to_OSN <- unlist_neighbors_x_list_to[!(unlist_neighbors_x_list_to %in% MM_neighbors_list[[i]])] 
      unlist_neighbors_x_list_to_ISN <- unlist_neighbors_x_list_to[(unlist_neighbors_x_list_to %in% MM_neighbors_list[[i]])]       
      
      num_magic_to_OSN <-  length(which(unlist_neighbors_x_list_to_OSN %in% magic_module_genes$EntrezID))
      num_RO_to_OSN <- length(which(!(unlist_neighbors_x_list_to_OSN %in% magic_module_genes$EntrezID))) 
      ratio_to_OSN = ifelse(num_RO_to_OSN == 0, ifelse(num_magic_to_OSN == 0, NaN, Inf), round(num_magic_to_OSN/num_RO_to_OSN, 4))
      
      num_MM_to_ISN <- length(unlist_neighbors_x_list_to_ISN %in% magic_module_genes$EntrezID) 
    }
    
    vect <- c(num_genes_SN, num_magic_from_OSN, num_RO_from_OSN, ratio_from_OSN, num_MM_from_ISN, num_magic_to_OSN, num_RO_to_OSN, ratio_to_OSN, num_MM_to_ISN)
    rm(num_genes_SN, num_magic_from_OSN, num_RO_from_OSN, ratio_from_OSN, num_MM_from_ISN, num_magic_to_OSN, num_RO_to_OSN, ratio_to_OSN, num_MM_to_ISN)
    temp_df <- rbind.data.frame(temp_df, vect)
    
    row.names(temp_df)[nrow(temp_df)] <- paste(MM_neighbors_list[[i]][1], "d", num_neighbors, sep= "_")
  }
  
  temp_df <- temp_df[-1,]
  
  #Sort by highest connections to MM
  #   temp_df <- temp_df[order(temp_df[,4],  decreasing = TRUE),]
  
  return(list(temp_df, MM_neighbors_list, MM_neighbors))
}
#************************************************************************************

####Generate.MMandROCombined.SN Function####
#Generates subnetworks for Situation one MM > P > A > B > RO
Generate.MMandROCombined.SN <- function(temp_df_B_SNs, temp_df_P_SNs, num_neighbors) {
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
      temp_combined_Subnetworks[input_index, B_index] <- PAB_SN[input_index]
      
    }
    colnames(temp_combined_Subnetworks)[B_index] <- paste(colnames(temp_df_B_SNs)[B_index], "g", length(PAB_SN), sep="_")
    rm(P_SN, B_SN, Individual_B_from_A, PAB_SN)
  }
  
  return(temp_combined_Subnetworks)
}
#************************************************************************************
