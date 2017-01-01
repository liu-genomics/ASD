# parallele version for computing CG using multiple cores. 
# after loading library(data.table)
#The following object is masked from ‘package:reshape’:
  
#  melt
library(parallel)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
print(args)
infile <- args[1] # the bed files that needs to have CG content calculted 
node_n <- min(as.numeric(args[2]), detectCores()-2) # the second argument is the number of nodes that will be used, maximal is the available minus 2. 



#calculate CG content for 50bp, 100bp 200bp and 500bp window, 
calculate_cg <-function(seq){
  letter = as.data.frame(strsplit(c(seq),split=""))
  GC_count = length(letter[letter[,1] == "C"|letter[,1] == "c" | letter[,1] == "G" | letter[,1] == "g",1])
  GC_pct = GC_count/length(letter[,1]) 
  GC_pct
}

temp = read.delim(infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
print("Finished reading input file.")

cl <- makeCluster(node_n)
clusterExport(cl, "calculate_cg")
cg_pct = data.frame(temp[,1],as.vector(unlist(parLapply(cl, temp[,2], calculate_cg))))

stopCluster(cl)

print("Finished calculating cg content")

write.table(cg_pct, paste(infile, ".cg", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
quit("no")