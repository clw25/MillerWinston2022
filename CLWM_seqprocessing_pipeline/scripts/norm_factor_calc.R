library(tidyr)
library(dplyr)
#read in sample summary file
summary_file = read.table(snakemake@input[[1]], stringsAsFactors = F, header = F, sep=" ")
#remove path and .bam extension
samples <- data.frame(lapply(summary_file,function(x){
  gsub(".bam","",x)
}))
samples <- data.frame(lapply(samples,function(x){
  gsub("fastq/aligned/","",x)
}))

#separate samples into indivdual columns
colnames(samples) = c("Sample","Reads")
samples = samples %>% separate(Sample, c('genotype','replicate','type'), sep = "-")
samples = samples %>% separate(type, c('type','species'), sep = "_")

#add additional columns for normalization factors (alpha)
summary_file = as.data.frame(summary_file[,1])
summary_file[,(ncol(summary_file)+1):(ncol(summary_file)+2)] = NA
colnames(summary_file)= c('file','alpha_SI','alpha_RPM')

#calculate normalization factors
for (i in 1:nrow(samples))
#for RPM calculation
{alpha_rpm = (1/as.numeric(as.character(samples$Reads[i]))) * 1000000
summary_file$alpha_RPM[i] = alpha_rpm
#for spike-in (SI) normalization
#inputs = lib size normaliztion
{if ((samples[i,3]=='input'| samples[i,3]=='RNHinput') & samples[i,4]=='Scer')
{alpha = ((1/as.numeric(as.character(samples$Reads[i]))) * 1000000)} 
#IP or RNH = 1/(reads-IP-spike-in * (reads-input-experimental/reads-input-spikein))
  else if ((samples[i,3]=='IP' | samples[i,3]=='RNH') & samples[i,4]=='Scer')
{IP_Spom = as.numeric(as.character(samples$Reads[samples$genotype == samples$genotype[i] & samples$replicate == samples$replicate[i] & samples$type == samples$type[i] & samples$species =='Spom']))
  input_Spom= as.numeric(as.character(samples$Reads[samples$genotype == samples$genotype[i] & samples$replicate == samples$replicate[i] & samples$type == 'input' & samples$species =='Spom']))
  input_Scer= as.numeric(as.character(samples$Reads[samples$genotype == samples$genotype[i] & samples$replicate == samples$replicate[i] & samples$type == 'input' & samples$species =='Scer']))
  alpha = ((1/(IP_Spom*(input_Scer/input_Spom)))*1000000)}
  else {alpha= (0)}
  summary_file$alpha_SI[i] = alpha
}}

write.table(summary_file[summary_file$alpha_SI != 0,], snakemake@output[[1]], row.names=F, quote=F, sep="\t")

