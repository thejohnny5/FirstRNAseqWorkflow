#Import DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
#BiocManager::install(version="3.13")
if (!"DESeq2" %in% rownames(installed.packages()))
  BiocManager::install("DESeq2")
suppressPackageStartupMessages(library('DESeq2'))

infile = snakemake@input[["raw_counts"]]


conds = scan(infile, nlines = 1, what = character())
conds = as.vector(unlist(conds))
#Read in counts
counts = read.table(infile, header=TRUE, sep="\t", quote="", row.names=1, skip = 1)

# Some tools generate the estimated counts as real numbers
# DESeq 2 allows only integers. We need to convert real numbers to rounded integers.
numeric_idx = sapply(counts, mode) == 'numeric'
counts[numeric_idx] = round(counts[numeric_idx], 0)

#Build the dataset
samples = names(counts)
condition = factor(conds)
colData = data.frame(samples=samples, condition=condition)

#Create DESeq2 dataset
dds = DESeqDataSetFromMatrix(countData=counts, colData=colData, design = ~condition)

#Set the reference chanmel
reference_channel=conds[1]
#dds$condition = relevel(dds$condition, reference_channel)
condition_names = levels(dds$condition)
#condition_names = c("EB", "EC", "EN", "M", "BiPS2")
# Run deseq2.
dds = DESeq(dds, betaPrior = FALSE)

#Set up a matrix to get average counts for each feature
meanMatrix = array(numeric(), c(nrow(counts), length(condition_names)))
rownames(meanMatrix) = rownames(counts)
colnames(meanMatrix) = paste(condition_names, '_baseMean', sep="")

#Get normalized counts for each sample
normed = counts(dds, normalized=TRUE)
normed = round(normed, 1)
start = 1
stop = 0


for (i in condition_names){
  mean_columns = colData$samples[colData$condition == i]
  mean_columns = as.character(mean_columns)
  print(i)
  if (length(mean_columns)>1){
    print(i)
    meanMatrix[, paste(i, "_baseMean", sep="")] = rowMeans(normed[,mean_columns])
  }else{
    meanMatrix[, paste(i, "_baseMean", sep="")] = normed[,mean_columns]
  }
}

prefix = c("Geneid\t")


cat(prefix, file=snakemake@output[["norm_counts_mean"]], sep='\t', append=F)
#Save the baseMeans of all features
write.table(data.frame(meanMatrix), snakemake@output[["norm_counts_mean"]], quote=FALSE, sep='\t', append=T)


############################################
# Pairwise analysis between all conditions #
############################################
my_data = list()
out_dir = snakemake@output[["pairwise_dir"]]

dir.create(out_dir, recursive=T)

for (i in 1:length(condition_names)){
  if (i==length(condition_names)){break()}
  else{
    for (j in (i+1):length(condition_names)){
      #Print Conditions being compared
      #if (condition_names[j] == 'BiPS2'){
      print(condition_names[i])
      print(condition_names[j])
      print("")
      #Set name of two conditions being compared
      name = paste(condition_names[i], condition_names[j], sep='.')
      
      #obtain DeSeq results of pairwise comparison in dataframe
      my_data[[name]] = data.frame(results(dds, c("condition", condition_names[i], condition_names[j])))
      
      #Rename values
      names(my_data[[name]])[names(my_data[[name]])=="pvalue"] = "PValue"
      names(my_data[[name]])[names(my_data[[name]])=="padj"] = "FDR"
      
      #Return fold change from log2foldchange
      my_data[[name]]$foldChange = 2 ^ my_data[[name]]$log2FoldChange
      
      #Adjust p value
      my_data[[name]]$Padj = p.adjust(my_data[[name]]$PValue, method="hochberg")
      
      #Sort for output
      my_data[[name]] = my_data[[name]][with(my_data[[name]], order(PValue, -foldChange)), ]
      
      #calculate false positive
      my_data[[name]]$falsePos = 1:nrow(my_data[[name]]) * my_data[[name]]$FDR
      
      #Rearrange column order
      my_data[[name]] = my_data[[name]][c(1, 7, 2, 3, 4, 5, 8, 6, 9)]
      
      #Ensure Pvalue and Padj are in scientific notation
      my_data[[name]]$PValue = format(my_data[[name]]$PValue, scientific = TRUE)
      #my_data[[name]]$Padj = format(my_data[[name]]$Padj, scientific = TRUE)
      
      #Rename columns to contain the conditions being compared
      colnames(my_data[[name]]) = paste(colnames(my_data[[name]]), name, sep='-')
      
      #Save the files for each pairwise comparison
      cat(prefix, file=paste(out_dir, '/', name, '.tsv', sep=''), sep='\t', append=F)
      write.table(my_data[[name]], paste(out_dir, '/', name, '.tsv', sep = ''), quote=FALSE, sep='\t', append=T)}
    #}
  }
}

#Write counts to tsv
#Write double header
conds = c('', conds, '\nGeneid\t')
cat(conds, file = snakemake@output[["norm_counts"]], sep='\t', append=F)
write.table(normed, snakemake@output[["norm_counts"]], quote=FALSE, sep='\t', append=T)

#Write log2 to tsv
logs = log2(normed + 1)
cat(conds, file = snakemake@output[["norm_counts_log2"]], sep='\t', append=F)
write.table(logs, snakemake@output[["norm_counts_log2"]], quote=FALSE, sep='\t', append=T)

log2meanMat = log2(meanMatrix + 1)
cat(prefix, file=snakemake@output[["norm_counts_mean_log2"]], sep='\t', append=F)
write.table(log2meanMat, snakemake@output[["norm_counts_mean_log2"]], quote=FALSE, sep='\t', append=T)

