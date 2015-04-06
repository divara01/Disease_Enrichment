####Load Omental Co-Expression Data####
all_omental_coexpression <- read.delim("~/Google Drive/PhD Research/Research Files/Data Acquisition/human_mgh_omental_all-Q4_tommodinfo.txt", stringsAsFactors=FALSE)
omental_coexpression_gene_modules <- all_omental_coexpression[,c(3, 5, 23)]
length(na.omit(omental_coexpression_gene_modules[,2]))
table(omental_coexpression_gene_modules$module)
#black         blue        brown         cyan        green  greenyellow         grey      magenta 
#  240         1257         1238           61          712           90         1850          178 
#midnightblue         pink       purple          red       salmon          tan    turquoise       yellow 
#          59          190          140          411           71           72         2731          859 

####Identify Omental Genes in Coexpression Data####
omental_modules <- omental_coexpression_gene_modules[na.omit(match(V(g_whole_omental)$name, omental_coexpression_gene_modules$LocusLinkID)),]
table(omental_modules$module)
#black         blue        brown         cyan        green  greenyellow         grey      magenta 
#   14          544          769           47          494           19          389           85 
#midnightblue         pink       purple          red       salmon          tan    turquoise       yellow 
#           5          124           77          275           44            2         1123          351 

####Identify MM Genes in Coexpression Data####
magic_module_genes <- read.delim("~/Google Drive/PhD Research/Joel Dudley/Research Files/Network Analysis/magic_module_genes_updated.txt", na.strings="", stringsAsFactors=FALSE)

magic_module_coexpression <- omental_coexpression_gene_modules[na.omit(match(V(g_MM)$name, omental_coexpression_gene_modules$LocusLinkID)),]
table(magic_module_coexpression$module)
#     blue     brown     green       red turquoise    yellow 
#       1         4       442         1         2         2 

####Identify RO Genes in Coexpression Data####
RO_coexpression_modules <- omental_coexpression_gene_modules[match(V(g_RO)$name, omental_coexpression_gene_modules$LocusLinkID),]

RO_coexpression_genes <- RO_coexpression_modules$LocusLinkID

table(RO_coexpression_modules$module)
#black         blue        brown         cyan        green  greenyellow         grey      magenta 
#   14          543          765           47           52           19          389           85 
#midnightblue         pink       purple          red       salmon          tan    turquoise       yellow 
#           5          124           77          274           44            2         1121          349 

####Store Whole Omental genes according to module in a DF (generate SN) + Disease Enrichment on each module####
modules <- unique(omental_modules$module)

coexpression_modules <- data.frame(matrix(NA, nrow = max(table(omental_modules$module)), ncol = length(modules)))
colnames(coexpression_modules) <- modules

for(col in 1:ncol(coexpression_modules)){
  mod <- colnames(coexpression_modules)[col]
  row = 1
  for(index in 1:nrow(omental_modules)){
    temp_mod <- omental_modules$module[index]
    if (temp_mod == mod){
      coexpression_modules[row, col] <- omental_modules$LocusLinkID[index]
      row = row + 1
    }
  }
  len <- length(which(!is.na(coexpression_modules[,col])))
  colnames(coexpression_modules)[col] <- paste(colnames(coexpression_modules)[col], "g", len, sep = "_")
}

filename <- "Coexpression_Module_Entire_omental_Disease_Enrichments.csv"
subnetworks <- coexpression_modules

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
        Overlap <- length(which(!is.na(V(g_whole_omental)$name[match(GeneSetTable[, j], sample(V(g_whole_omental)$name, TargetSize, replace = FALSE, prob = NULL))]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(GeneSetTable[match(sample(V(g_whole_omental)$name, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
      }
      TestSetSize <- length(which(!is.na(match(V(g_whole_omental)$name, GeneSetTable[,j]))))
      TotalSize <- length(V(g_whole_omental)$name)
      
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
    TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_whole_omental)$name))))
    
    TotalSize <- length(V(g_whole_omental)$name)
    index <- (i * numSN) - (numSN - j) 
    temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
    Disease_p_values[j, 1] <- colnames(GeneSetTable)[j]
    Disease_p_values[j, 2] <- TestSetSize
    Disease_p_values[j, 3] <- colnames(subnetworks)[i]
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
output <- output[, -1]
TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]

rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)

length(which(Disease_sample_p_values < 0.05))
length(which(Disease_p_values <= 0.05))

filename_Tableau <- paste("Tableau_", filename, sep = "")
write.csv(TableauDiseaseEnrichment, file = filename_Tableau)

write.csv(output, file = filename)

####Store MM genes according to module in a DF (generate SN) + Disease Enrichment on each module####
filename = "MM_CoexpressionModuleDiseaseEnrichments.csv"
modules <- unique(magic_module_coexpression$module)

coexpression_modules <- data.frame(matrix(NA, nrow = max(table(magic_module_coexpression$module)), ncol = length(modules)))
colnames(coexpression_modules) <- modules

for(col in 1:ncol(coexpression_modules)){
  mod <- colnames(coexpression_modules)[col]
  row = 1
  for(index in 1:nrow(magic_module_coexpression)){
    temp_mod <- magic_module_coexpression$module[index]
    if (temp_mod == mod){
      coexpression_modules[row, col] <- magic_module_coexpression$LocusLinkID[index]
      row = row + 1
    }
  }
  len <- length(which(!is.na(coexpression_modules[,col])))
  colnames(coexpression_modules)[col] <- paste(colnames(coexpression_modules)[col], "g", len, sep = "_")
}

subnetworks <- coexpression_modules

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
        Overlap <- length(which(!is.na(V(g_MM)$name[match(GeneSetTable[, j], sample(V(g_MM)$name, TargetSize, replace = FALSE, prob = NULL))]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(GeneSetTable[match(sample(V(g_MM)$name, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
      }
      TestSetSize <- length(which(!is.na(match(V(g_MM)$name, GeneSetTable[,j]))))
      TotalSize <- length(V(g_MM)$name)
      
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
    TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_MM)$name))))
    
    TotalSize <- length(V(g_MM)$name)
    index <- (i * numSN) - (numSN - j) 
    temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
    Disease_p_values[j, 1] <- colnames(GeneSetTable)[j]
    Disease_p_values[j, 2] <- TestSetSize
    Disease_p_values[j, 3] <- colnames(subnetworks)[i]
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
    print(paste("j=",j, sep = ""))
  }
  print(paste("i=",i, sep = ""))
  TableauDiseaseEnrichment <- rbind.data.frame(TableauDiseaseEnrichment, Disease_p_values)
}

output <- cbind(output, output2)
output <- output[, -1]
TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]

rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)

length(which(Disease_sample_p_values < 0.05))
length(which(Disease_p_values <= 0.05))

filename_Tableau <- paste("Tableau_", filename, sep = "")
write.csv(TableauDiseaseEnrichment, file = filename_Tableau)

write.csv(output, file = filename)

####Store RO genes according to module in a DF (generate SN) + Disease Enrichment on each module####
filename = "RO_CoexpressionModuleDiseaseEnrichments.csv"
modules <- unique(RO_coexpression_modules$module)

coexpression_modules <- data.frame(matrix(NA, nrow = max(table(RO_coexpression_modules$module)), ncol = length(modules)))
colnames(coexpression_modules) <- modules

for(col in 1:ncol(coexpression_modules)){
  mod <- colnames(coexpression_modules)[col]
  row = 1
  for(index in 1:nrow(RO_coexpression_modules)){
    temp_mod <- RO_coexpression_modules$module[index]
    if (temp_mod == mod){
      coexpression_modules[row, col] <- RO_coexpression_modules$LocusLinkID[index]
      row = row + 1
    }
  }
  len <- length(which(!is.na(coexpression_modules[,col])))
  colnames(coexpression_modules)[col] <- paste(colnames(coexpression_modules)[col], "g", len, sep = "_")
}

subnetworks <- coexpression_modules

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
        Overlap <- length(which(!is.na(V(g_RO)$name[match(GeneSetTable[, j], sample(V(g_RO)$name, TargetSize, replace = FALSE, prob = NULL))]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(GeneSetTable[match(sample(V(g_RO)$name, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
      }
      TestSetSize <- length(which(!is.na(match(V(g_RO)$name, GeneSetTable[,j]))))
      TotalSize <- length(V(g_RO)$name)
      
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
    Disease_p_values[j, 3] <- colnames(subnetworks)[i]
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
    print(paste("j=",j, sep = ""))
  }
  print(paste("i=",i, sep = ""))
  TableauDiseaseEnrichment <- rbind.data.frame(TableauDiseaseEnrichment, Disease_p_values)
}

output <- cbind(output, output2)
output <- output[, -1]
TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]

rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)

length(which(Disease_sample_p_values < 0.05))
length(which(Disease_p_values <= 0.05))

filename_Tableau <- paste("Tableau_", filename, sep = "")
write.csv(TableauDiseaseEnrichment, file = filename_Tableau)

write.csv(output, file = filename)

####Load combined SNs for each disease in PAB/BAP analysis####
#These are SN that combine all SN that are enriched for each disease
#tag <- "PAB"
#CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/PAB_CombinedSN_SN_d1.csv", na.strings="", stringsAsFactors=FALSE)
#degree <- 1
#CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/PAB_CombinedSN_SN_d3.csv", na.strings="", stringsAsFactors=FALSE)
#degree <- 3
#CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/PAB_CombinedSN_SN_d5.csv", na.strings="", stringsAsFactors=FALSE)
#degree <- 5

tag <- "BAP"
#CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/BAP_CombinedSN_SN_d1.csv", na.strings="", stringsAsFactors=FALSE)
#degree <- 1
#CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/BAP_CombinedSN_SN_d3.csv", na.strings="", stringsAsFactors=FALSE)
#degree <- 3
CombinedSN <- read.csv("~/Google Drive/PhD Research/Research Files/Data Acquisition/BAP_CombinedSN_SN_d5.csv", na.strings="", stringsAsFactors=FALSE)
degree <- 5

GeneSetTable <- data.frame(matrix(NA, nrow = nrow(CombinedSN), ncol = ncol(CombinedSN)))
colnames(GeneSetTable) <- colnames(CombinedSN)
for(i in 1:ncol(CombinedSN)){
  unique_genes <- unique(CombinedSN[,i])
  for(j in 1:length(unique_genes)){
    GeneSetTable[j, i] <- unique_genes[j]
  }
}

#Enrich coexpression modules for SN enriched for diseases#
modules <- unique(omental_modules$module)

coexpression_modules <- data.frame(matrix(NA, nrow = max(table(omental_modules$module)), ncol = length(modules)))
colnames(coexpression_modules) <- modules

for(col in 1:ncol(coexpression_modules)){
  mod <- colnames(coexpression_modules)[col]
  row = 1
  for(index in 1:nrow(omental_modules)){
    temp_mod <- omental_modules$module[index]
    if (temp_mod == mod){
      coexpression_modules[row, col] <- omental_modules$LocusLinkID[index]
      row = row + 1
    }
  }
  len <- length(which(!is.na(coexpression_modules[,col])))
  colnames(coexpression_modules)[col] <- paste(colnames(coexpression_modules)[col], "g", len, sep = "_")
}

filename <- paste(tag, "_Disease_Specific_SN_Coexpression_Module_Disease_Enrichments_d_", degree, ".csv", sep = "")
subnetworks <- coexpression_modules

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
        Overlap <- length(which(!is.na(V(g_whole_omental)$name[match(GeneSetTable[, j], sample(V(g_whole_omental)$name, TargetSize, replace = FALSE, prob = NULL))]))) 
      }
      if (length(na.omit(GeneSetTable[,j])) <= length(na.omit(subnetworks[,k]))) {
        Overlap <- length(which(!is.na(GeneSetTable[match(sample(V(g_whole_omental)$name, TargetSize, replace = FALSE, prob = NULL), GeneSetTable[, j]), j])))
      }
      TestSetSize <- length(which(!is.na(match(V(g_whole_omental)$name, GeneSetTable[,j]))))
      TotalSize <- length(V(g_whole_omental)$name)
      
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
    TestSetSize <- length(which(!is.na(match(GeneSetTable[,j], V(g_whole_omental)$name))))
    
    TotalSize <- length(V(g_whole_omental)$name)
    index <- (i * numSN) - (numSN - j) 
    temp_list <- OverEnrichment(Overlap, TargetSize, TestSetSize, TotalSize)
    Disease_p_values[j, 1] <- paste(colnames(GeneSetTable)[j], "g", length(na.omit(GeneSetTable[,j])), sep = "_")
    Disease_p_values[j, 2] <- TestSetSize
    Disease_p_values[j, 3] <- colnames(subnetworks)[i]
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
output <- output[, -1]
TableauDiseaseEnrichment <- TableauDiseaseEnrichment[-1, ]

rm(i, j, k, val, Overlap, TargetSize, TestSetSize, TotalSize)

length(which(Disease_sample_p_values < 0.05))
length(which(Disease_p_values <= 0.05))

filename_Tableau <- paste("Tableau_", filename, sep = "")
write.csv(TableauDiseaseEnrichment, file = filename_Tableau)

write.csv(output, file = filename)