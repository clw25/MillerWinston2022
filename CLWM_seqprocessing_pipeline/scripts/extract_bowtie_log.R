

samples = read.table("barcodes.tsv", stringsAsFactors = F, header = F, sep='\t')

results = NULL
stats = NULL
for (i in 1:nrow(samples))
{
  #read log file for each sample 
  filename = paste("logs/alignment.",samples[i,1],".log",sep='')
  log = read.table(filename, header = F, sep = '\t', stringsAsFactors = F)
  #get total reads
  temp = strsplit(log[1,],' ')
  temp2 = unlist(temp)
  total = as.numeric(temp2[1])
  value = NULL
  percent = NULL
  #get unmapped, unique mappers, and multipmappers both the number of reads (value) and percent of total (percent)
  for (j in 3:5){
    temp = strsplit(log[j,],'\\(')
    temp2 = unlist(temp)
    value = cbind(value, as.numeric(temp2[1]))
    percent = cbind(percent, round(((as.numeric(temp2[1])/total)*100), digits = 2))
  }
  #combine all samples together into results table
  results = rbind(results,cbind(total,value,percent))
}
#add sample names and column names
results = cbind(samples[1],results) 
colnames(results) = c('Sample','total', 'unmapped', 'unique','multimappers','unmapped %', 'unique %', 'multimapper %')

write.table(results, snakemake@output[[1]], row.names = F, quote = F, sep = '\t')
