source('ScaleR.R')
library(dplyr)
library(varhandle)

# Import test data and select subset (11 columns and 12 rows)
m <- 'lcms-disease'# write data type here

if(m == 'rna-age'){
  data <- t(read.table(file = '/ScaleR/Data_Pub/TCGA-OV.htseq_fpkm.tsv', sep = '\t', header = FALSE))
  colnames(data) <- data[1,]
  data <- data[-c(1),]
  data <- as.data.frame(data)
  data[,-1] <- lapply(data[,-1], as.numeric)

  Pheno.Data <- read.table(file = '/ScaleR/Data_Pub//OV_clinicalMatrix', sep = '\t', header = TRUE)
  Clin <- NULL
  Clin$ID <- Pheno.Data$sampleID
  Clin$Age <- Pheno.Data$age_at_initial_pathologic_diagnosis
  Clin$stage <- Pheno.Data$clinical_stage
  Clin$grade <- Pheno.Data$neoplasm_histologic_grade
  Clin <- as.data.frame(Clin)
  Clin <- (na.omit(Clin))

  data$Ensembl_ID <- gsub('.{1}$', '', data$Ensembl_ID)

  Data.Full <- merge(Clin, data, by.x='ID', by.y='Ensembl_ID')
  age.data <- na.omit(Data.Full)
  x <- age.data[,5:dim(data)[2]]
  y <- age.data$Age
  gamma = seq(0,1,by=0.05)
  alpha = seq(0,1,by=0.1)
  ncomp = c(1:10)

  print("Running PLS...")
  results <- ScaleR(x, y, model='pls', gamma=gamma, nfolds=10, ncomp=ncomp)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rna-seq_age_pls.csv")

  print("Running Lasso...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=1, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rna-seq_age_lasso.csv")

  print("Running Ridge...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=0, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rna-seq_age_ridge.csv")

  print("Running Elastic Nets...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=alpha, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rna-seq_age_enets.csv")
}

if(m == 'micro-age'){
  Micro_Data <- t(read.table(file = '/ScaleR/Data_Pub/HT_HG-U133A', sep = '\t', header = FALSE))
  colnames(Micro_Data) <- Micro_Data[1,]
  Micro_Data <- Micro_Data[-c(1),]
  Micro_Data <- as.data.frame(Micro_Data)
  Micro_Data[,-1] <- lapply(Micro_Data[,-1], as.numeric)

  Pheno.Data <- read.table(file = '/ScaleR/Data_Pub//OV_clinicalMatrix', sep = '\t', header = TRUE)
  Clin <- NULL
  Clin$ID <- Pheno.Data$sampleID
  Clin$Age <- Pheno.Data$age_at_initial_pathologic_diagnosis
  Clin$stage <- Pheno.Data$clinical_stage
  Clin$grade <- Pheno.Data$neoplasm_histologic_grade
  Clin <- as.data.frame(Clin)
  Clin <- (na.omit(Clin))

  Data.Full <- merge(Clin, Micro_Data, by.x='ID', by.y='sample')
  Data.Full <- na.omit(Data.Full)

  x <- Data.Full[,5:dim(Data.Full)[2]]
  y <- Data.Full$Age
  gamma = seq(0,1,by=0.05)
  alpha = seq(0,1,by=0.1)
  ncomp = c(1:10)

  print("Running PLS...")
  results <- ScaleR(x, y, model='pls', gamma=gamma, nfolds=10, ncomp=ncomp)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/micro_age_pls.csv")

  print("Running Lasso...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=1, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/micro_age_lasso.csv")

  print("Running Ridge...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=0, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/micro_age_ridge.csv")

  print("Running Elastic Nets...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=alpha, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/micro_age_enets.csv")
}

if(m == 'rppa-age'){
  RPPA.Data <- t(read.table(file = '/ScaleR/Data_Pub/RPPA_RBN', sep = '\t', header = FALSE))
  colnames(RPPA.Data) <- RPPA.Data[1,]
  RPPA.Data <- RPPA.Data[-c(1),]
  Pheno.Data <- read.table(file = '/ScaleR/Data_Pub//OV_clinicalMatrix', sep = '\t', header = TRUE)
  Clin <- NULL
  Clin$ID <- Pheno.Data$sampleID
  Clin$Age <- Pheno.Data$age_at_initial_pathologic_diagnosis
  Clin$stage <- Pheno.Data$clinical_stage
  Clin$grade <- Pheno.Data$neoplasm_histologic_grade
  Clin <- as.data.frame(Clin)
  Clin <- (na.omit(Clin))
  Data.Full <- merge(Clin, RPPA.Data, by.x='ID', by.y='Sample_description')
  Data.Full <- na.omit(Data.Full)
  Data.Full <- unfactor(Data.Full)
  x <- Data.Full[,5:dim(Data.Full)[2]]
  y <- Data.Full$Age

  gamma = seq(0,1,by=0.05)
  alpha = seq(0,1,by=0.1)
  ncomp = c(1:10)

  print("Running PLS...")
  results <- ScaleR(x, y, model='pls', gamma=gamma, nfolds=10, ncomp=ncomp)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rppa_age_pls.csv")

  print("Running Lasso...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=1, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rppa_age_lasso.csv")

  print("Running Ridge...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=0, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rppa_age_ridge.csv")

  print("Running Elastic Nets...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=alpha, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/rppa_age_enets.csv")
}

if(m == 'lcms-disease'){

  LCMS.Data <- read.table(file = '/Users/sammoney-kyrle/gdrive/Imperial/scaleR/ScaleR/Data_Pub/m_MTBLS2266_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv', sep = '\t', header = TRUE)
  LCMS <- t(LCMS.Data[21:ncol(LCMS.Data)])
  colnames(LCMS) <-  paste(LCMS.Data$retention_time, LCMS.Data$mass_to_charge, sep='_')
  LCMS <- LCMS[-c(1),]
  LCMS <- as.data.frame(LCMS)
  LCMS$ID <- gsub('.{5}$', '', rownames(LCMS))
  LCMS$ID <- gsub('X', '', LCMS$ID)

  Pheno <- read.delim(file = '/Users/sammoney-kyrle/gdrive/Imperial/scaleR/ScaleR/Data_Pub//s_MTBLS2266.txt')

  Clin <- NULL
  Clin$Sample.Name <- Pheno$Sample.Name
  Clin$Disease <- Pheno$Characteristics.Disease.
  Clin <- as.data.frame(Clin)
  Data.Full <- merge(Clin, LCMS, by.x='Sample.Name', by.y='ID')
  normalize.medFC <- function(X) {
    # Perform median fold change normalisation
    #           X - data set [Variables & Samples]
    medSam <- apply(X,1,median)
    count<-dim(X)[2]
    medSam[which(medSam==0)] <- 0.01
    for (i in 1:count) {
        medFCiSmpl <- X[,i]/medSam
       # sel                               <- which(!is.na(as.character(medFCiSmpl)))
        X[,i]                             <- X[,i]/median(medFCiSmpl)
        }
    return (X)
  }

  Data.Full[3:ncol(Data.Full)] <- normalize.medFC(Data.Full[3:ncol(Data.Full)])
  Data.Full <- as.data.frame(t(na.omit(t(Data.Full))))
  Data.Full <- unfactor(Data.Full)
  Data.Full$Disease <- as.numeric(as.factor(Data.Full$Disease))
  Data.Full <-Data.Full[!Data.Full$Disease == "3", ]
  Data.Full[,3:dim(Data.Full)[2]] <- apply(Data.Full[,3:dim(Data.Full)[2]], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
  x <- Data.Full[,3:dim(Data.Full)[2]]
  y <- Data.Full$Disease

  gamma = seq(0,1,by=0.05)
  alpha = seq(0,1,by=0.1)
  ncomp = c(1:10)

  print("Running PLS...")
  results <- ScaleR(x, y, model='pls', gamma=gamma, nfolds=10, ncomp=ncomp)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/lcms_disease_pls.csv")

  print("Running Lasso...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=1, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/lcms_disease_lasso.csv")

  print("Running Ridge...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=0, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/lcms_disease_ridge.csv")

  print("Running Elastic Nets...")
  results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=alpha, nfolds=10)
  print("Writing results...")
  write.csv(results, "/ScaleR/results/lcms_disease_enets.csv")
}
