#### args <- c("/Sample_data/sample_purity.txt","/intermediate/snv_input.txt", "/intermediate/middle/multiplicity.txt", "/Sample_data/Res/")

suppressMessages(library(dummies))
args <- commandArgs(trailingOnly=TRUE)
estimated_purity <- unlist(lapply(1:4, function(x){
	suma <- read.table(sprintf("/semifinal/res%d/subclonal_structure.txt", x))
	return(max(suma[,3]))
	}))
purity <- unlist(read.table(args[1]))
diff <- abs(estimated_purity - purity)
ind <- max(which.min(diff))

suma <- read.table(sprintf('/semifinal/res%d/subclonal_structure.txt',ind))
assign <- unlist(read.table(sprintf('/semifinal/res%d/mutation_assignments.txt',ind)))
all_mut <- read.table(args[2])
mut_out <- cbind(all_mut[,1:2],rep(-1,nrow(all_mut)),rep(0,nrow(all_mut)))
allMutComb <- paste(all_mut[,1], all_mut[,2], sep='_')
multi <- read.table(args[3])
assignedMutComb <- paste(multi[,1], multi[,2], sep = '_') 
correspondence <- match(assignedMutComb,allMutComb)
prev <- suma[match(assign,suma[,1]),3]
unique_assign <- unique(assign)
new_index <- 1:length(unique_assign)
new_assign <- new_index[match(assign, unique_assign)]
mut_out[correspondence,3] <- new_assign
mut_out[correspondence,4] <- prev
mut_out[which(mut_out[,3]==-1),3] <- length(unique_assign)+1

res_1B <- length(unique_assign)
res_1C <- matrix(0,nrow=length(unique(mut_out[,3])), ncol=3)
res_1C[,1] <- sort(unique(mut_out[,3]))
for(i in 1:nrow(res_1C)){
	index <- which(mut_out[,3]==res_1C[i,1])
	res_1C[i,2] <- length(index)
	res_1C[i,3] <- mean(mut_out[index,4])
}
res_2A <- mut_out[,3]
res_2B <- dummy(res_2A)
write.table(res_2B, file='/intermediate/2B_matrix',sep=',',row.names=F, col.names=F)
cmd <- sprintf('python3.6 compute_2B.py /intermediate/2B_matrix %s/2B.txt.gz', args[4])
a <- system(cmd,intern=TRUE)
write.table(res_1B, file=sprintf('%s/1B.txt',args[4]),sep='\t',row.names=F, col.names=F)
write.table(res_1C, file=sprintf('%s/1C.txt',args[4]),sep='\t',row.names=F, col.names=F)
write.table(res_2A, file=sprintf('%s/2A.txt',args[4]),sep='\t',row.names=F, col.names=F)