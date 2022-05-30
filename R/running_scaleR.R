source('ScaleR.R')

# Import test data and select subset (11 columns and 12 rows)
data <- t(read.table(file = 'Data_Pub/TCGA-OV.htseq_fpkm.tsv', sep = '\t', header = FALSE))
colnames(data) <- data[1,]
data <- data[-c(1),]
data <- as.data.frame(data)
# data <- data[1:12,]
data[,-1] <- lapply(data[,-1], as.numeric)

Pheno.Data <- read.table(file = 'Data_Pub//OV_clinicalMatrix', sep = '\t', header = TRUE)
Clin <- NULL
Clin$ID <- Pheno.Data$sampleID
Clin$Age <- Pheno.Data$age_at_initial_pathologic_diagnosis
Clin$stage <- Pheno.Data$clinical_stage
Clin$grade <- Pheno.Data$neoplasm_histologic_grade
Clin <- as.data.frame(Clin)
Clin <- (na.omit(Clin))

data$Ensembl_ID <- gsub('.{1}$', '', data$Ensembl_ID)

Data.Full <- merge(Clin, data, by.x='ID', by.y='Ensembl_ID')
data <- na.omit(Data.Full)
x <- data[,5:dim(data)[2]]
y <- data$Age
gamma = seq(0,1,by=0.05)
alpha = seq(0,1,by=0.1)
ncomp = 1:10


print("Running PLS...")
results <- ScaleR(x, y, model='pls', gamma=gamma, nfolds=10, ncomp=ncomp)
print("Writing results...")
write.csv(results, "results/pls_rna-seq_age_results.csv")

print("Running Lasso...")
results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=1, nfolds=10)
print("Writing results...")
write.csv(results, "/results/lasso_rna-seq_age_results.csv")

print("Running Ridge...")
results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=0, nfolds=10)
print("Writing results...")
write.csv(results, "/results/ridge_rna-seq_age_results.csv")

print("Running Elastic Nets...")
results <- ScaleR(x, y, model='glm', gamma=gamma, alpha=alpha, nfolds=10)
print("Writing results...")
write.csv(results, "/results/enets_rna-seq_age_results.csv")
