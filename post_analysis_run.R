#------------------------------------------------------------#
# This script takes care of the postprocessing for CliP
# the postprocessing mostly taking care of the hamonization of the results for both downsampled/non-downsampled
# and make some further filtering which may not be necessary and you can turn them off by setting filtering to 0.
# Typically you do not need to do the filtering, but when the average coverage is low, say 30-40X, you may want to do this extra filtering.
# Initialized by Kaixian Yu
# Date: 05/04/2018
# Email: kaixiany@163.com
#------------------------------------------------------------# 
# The script takes commandline arguments: CliP_results_dir preprocessed_output_dir Output_dir lambda filtering_flag
# The results folder should contain results for only one sample and the non-downsampled version is always given higher priority. 
# The postprocessing is done per lam, the results of the same sample for different lambda can be stored in the same folder.
args <- commandArgs(trailingOnly = T)
library(parallel)
input.prefix        <- args[1]
preprocessed.prefix <- args[2]
output.prefix       <- args[3]
lam                 <- args[4]
filtering           <- as.numeric(args[5])
threshold           <- 0.1
clonalFrac          <- 0.3
least.ratio         <- 0.05
SAFETY              <- 0.05
N.REP               <- 30
RATIO.TH            <- 0.90
if(!dir.exists(output.prefix)){
	dir.create(output.prefix)
}
outfile1 <- sprintf("%s/lam%s_phi.txt", input.prefix, lam)
outfile2 <- sprintf("%s/lam%s_label.txt", input.prefix, lam)
outfile3 <- sprintf("lam%s_rep", lam)
ifsubsample <- list.files(input.prefix,pattern = outfile3,full.names=T)
if( file.exists(outfile1) || file.exists(outfile2) ){
	if(!file.exists(outfile1)){
		r <- unlist(read.table(sprintf("%s/r.txt", preprocessed.prefix), header=F))
		n <- unlist(read.table(sprintf("%s/n.txt", preprocessed.prefix), header=F))
		minor <- unlist(read.table(sprintf("%s/minor.txt",  preprocessed.prefix), header=F))
		total <- unlist(read.table(sprintf("%s/total.txt",  preprocessed.prefix), header=F))
		theta.hat <- r / n
		phi.hat <- ( 2) / ( minor / theta.hat - total + 2 )
		label <- unlist(read.table(outfile2,header=F))
		all.label <- unique(label)
		cluster.phi <- unlist(lapply(all.label,function(x){ 
					ind <- which(label == x)
					return( sum(phi.hat[ind]*n[ind])/sum(n[ind]) )
				}))
		phi <- unlist(lapply(label,function(x){ cluster.phi[which(all.label == x)] }))
		
	}
	
	if(!file.exists(outfile2)){
		phi <- unlist(read.table(outfile1,header=F))
		cluster.phi <- unique(phi)
		label <- unlist(lapply(phi,function(x){which(cluster.phi == x)-1} ))
		all.label <- unique(label)
		
	}
	label <- unlist(read.table(outfile2,header=F))
	all.label <- unique(label)
	phi <- unlist(read.table(outfile1,header=F))
	cluster.phi <- unique(phi)

} else if (length(ifsubsample) >0){
	tmp <- lapply(ifsubsample,function(x){read.table(x,header=F, sep="\t")})
	if(ncol(tmp[[1]])<3){
		tmp <- lapply(ifsubsample,function(x){read.table(x,header=F, sep=",")})
	}
	pooled.beta      <- list()
	pooled.pssm      <- list()
	pooled.beta[[1]] <- unlist(lapply(tmp,function(x){ unlist(x[,3]) } ))
	pooled.pssm[[1]] <- unlist(lapply(tmp,function(x){ unlist(x[,2])/sum(x[,2]) } ))
	re.cluster <- T
	while(re.cluster){
		re.cluster <- F
		if(length(unique(pooled.beta[[1]]))>1){
		consensus.ccf <- lapply(pooled.beta,function(x){
		#first determine if 1 cluster is enough
		flag  <- T
		N.clu <- 2
		while(flag){
			tmp.kmeans <-kmeans(x,N.clu,nstart = N.REP)
			#print(N.clu)
			N.clu <- N.clu + 1
			if(N.clu == 3){
				tmp.center  <- tmp.kmeans$centers
				if(abs(tmp.center[1]-tmp.center[2]) <= SAFETY ){
					flag <- F
					return(list(center = mean(x),size=length(x), cluster = tmp.kmeans$cluster))
				} else{
					if(tmp.kmeans$betweenss/tmp.kmeans$totss >= RATIO.TH){
					flag <- F
					return(list(center = as.numeric(tmp.kmeans$centers),size=tmp.kmeans$size, cluster = tmp.kmeans$cluster))
					}
				}
			} else{
				if(tmp.kmeans$betweenss/tmp.kmeans$totss >= RATIO.TH){
				flag <- F
				return(list(center = as.numeric(tmp.kmeans$centers),size=tmp.kmeans$size, cluster = tmp.kmeans$cluster))
				}
			}
		}
	})
	ind <- which(consensus.ccf[[1]]$size == 1)
	if(length(ind)>0)
	{
	  #value.to.exclude <- consensus.ccf[[1]]$center[ind]
	  ccf.ind          <- unlist(lapply(ind,function(x){which(consensus.ccf[[1]]$cluster == x)}))
	  #ccf.ind <- unlist(lapply(value.to.exclude, function(x){ which.min(abs(pooled.beta[[1]]-x)) } ))
	  re.cluster <- T
	  pooled.beta[[1]]<- pooled.beta[[1]][-ccf.ind]
	  pooled.pssm[[1]]<- pooled.pssm[[1]][-ccf.ind]
	}
	} else{
	  consensus.ccf[[1]]$center <- unique(pooled.beta[[1]])
	  consensus.ccf[[1]]$cluster <- rep(1,length(pooled.beta[[1]]))
	}
	}
	consensus.phi         <- unlist(consensus.ccf[[1]]$center)
	r <- unlist(read.table(sprintf("%s/r.txt", preprocessed.prefix), header=F))
	n <- unlist(read.table(sprintf("%s/n.txt", preprocessed.prefix), header=F))
	minor <- unlist(read.table(sprintf("%s/minor.txt",  preprocessed.prefix), header=F))
	total <- unlist(read.table(sprintf("%s/total.txt",  preprocessed.prefix), header=F))
	theta.hat <- r / n
	phi.hat <- ( 2) / ( minor / theta.hat - total + 2 )
	No.mutations <- length(r)
	consensus.pssm        <- unlist(tapply(pooled.pssm[[1]],consensus.ccf[[1]]$cluster, median))
	consensus.pssm        <- consensus.pssm / sum(consensus.pssm) 
	consensus.sd          <- sqrt( consensus.pssm * No.mutations *  (min(consensus.phi,0.99) )*(1 - min(0.99,consensus.phi)) )						
	Distance              <- Reduce("cbind",lapply(1:length(consensus.phi),function(x){ abs(phi.hat- consensus.phi[x]) / consensus.sd[x] }))
	if(length(ncol(Distance))>0){
		label <- apply(Distance, 1, which.min)	
	} else {
		label <- rep(1,length(Distance))							
	}
	cluster.phi <- consensus.phi
	all.label <- 1:length(cluster.phi)
	phi <- cluster.phi[label]
} else {
	#return(sprintf("No results found for %s",x))		
	stop(sprintf('No results found at the given location: %s', input.prefix))
}

#### Some filtering ####
suma <- matrix(0,ncol=3, nrow=length(cluster.phi))
for(i in 1:length(cluster.phi)){
	suma[i,1] = all.label[i]
	suma[i,2] = length(which(label==all.label[i]))
	suma[i,3] = cluster.phi[i]
}
if(filtering == 1){
	if(length(cluster.phi)==1){
		write.table(suma, file = sprintf("%s/subclonal_structure.txt",
                   output.prefix), quote = F, sep = "\t", 
            col.names = F, row.names = F )
		write.table(label, file = sprintf("%s/mutation_assignments.txt",
                                output.prefix), quote = F, sep = "\t", 
           col.names = F, row.names = F )
		   
		cat(sprintf('There is only one clone in this sample, no further filtering is needed.'))
	}else{
		suma <- suma[order(suma[,3],decreasing = F),]
		least.mut <- least.ratio * length(phi)
	   refine <- F
	   No.cls <- dim(suma)[1]
	   # filter 1: super clusters
	   if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac){
	     refine <- T
	   }
	   while(refine){
	     refine <- F
	     ind     <- c(which(label == suma[No.cls, 1]), which(label == suma[No.cls-1,1]))
	     new.phi <- (suma[No.cls, 2]*suma[No.cls, 3] + 
	                 suma[No.cls-1, 2]*suma[No.cls-1, 3]) / 
	                (suma[No.cls, 2] + suma[No.cls-1, 2])
	     new.label  <- suma[No.cls, 1]
	     phi[ind]   <- new.phi
	     label[ind] <- new.label
	     No.cls     <- No.cls - 1
	     suma <- matrix(NA, ncol = 3, nrow = No.cls)
	     suma[, 3] <- sort(unique(phi),decreasing = F)
	     suma[, 1] <- unlist(lapply(suma[, 3],
	                               function(y){
	                                 unique(label[which(phi == y)])
	                               }))
	     suma[, 2] <- unlist(lapply(suma[, 1],
	                                function(y){
	                                  length(which(label == y))
	                                }))
	     if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac){
	       refine <- T
	     }
	   }
   
	   # filter 2: small clone
	   if(suma[No.cls,2] / sum(suma[,2]) < clonalFrac && No.cls > 2){
	     refine <- T
	   }
	   while(refine){
	     refine <- F
	     ind     <- c(which(label == suma[No.cls, 1]), which(label == suma[No.cls-1,1]))
	     new.phi <- (suma[No.cls, 2]*suma[No.cls, 3] + 
	                   suma[No.cls-1, 2]*suma[No.cls-1, 3]) / 
	       (suma[No.cls, 2] + suma[No.cls-1, 2])
	     new.label  <- suma[No.cls, 1]
	     phi[ind]   <- new.phi
	     label[ind] <- new.label
	     No.cls     <- No.cls - 1
	     suma <- matrix(NA, ncol = 3, nrow = No.cls)
	     suma[, 3] <- sort(unique(phi),decreasing = F)
	     suma[, 1] <- unlist(lapply(suma[, 3],
	                                function(y){
	                                  unique(label[which(phi == y)])
	                                }))
	     suma[, 2] <- unlist(lapply(suma[, 1],
	                                function(y){
	                                  length(which(label == y))
	                                }))
	     if(suma[No.cls,2] / sum(suma[,2]) < clonalFrac && No.cls > 2){
	       refine <- T
	     }
	   }
	   # filter 3: too close
	   if( No.cls > 2 ){
	     diff <- suma[2:No.cls,3]-suma[1:(No.cls-1),3]
	     if(min(diff) < threshold ){
	       refine <- T
	     }
	   }
	     while(refine){
	       min.ind <- which.min(diff) + 1
	       refine <- F
	       ind     <- c(which(label == suma[min.ind, 1]), which(label == suma[min.ind-1,1]))
	       new.phi <- (suma[min.ind, 2]*suma[min.ind, 3] + 
	                     suma[min.ind-1, 2]*suma[min.ind-1, 3]) / 
	         (suma[min.ind, 2] + suma[min.ind-1, 2])
	       new.label  <- suma[min.ind, 1]
	       phi[ind]   <- new.phi
	       label[ind] <- new.label
	       No.cls     <- No.cls - 1
	       suma <- matrix(NA, ncol = 3, nrow = No.cls)
	       suma[, 3] <- sort(unique(phi),decreasing = F)
	       suma[, 1] <- unlist(lapply(suma[, 3],
	                                  function(y){
	                                    unique(label[which(phi == y)])
	                                  }))
	       suma[, 2] <- unlist(lapply(suma[, 1],
	                                  function(y){
	                                    length(which(label == y))
	                                  }))
	       if( No.cls > 2 ){
	         diff <- suma[2:No.cls,3]-suma[1:(No.cls-1),3]
	         if(min(diff) < threshold ){
	           refine <- T
	         }
	       }
	       
	       
	     }
	   #### filter 4: small subclone (as defined by at least least.mut mutations)
	   if( length(which(suma[,2]<= least.mut)) > 0  && No.cls > 2 && which(suma[,2]<= least.mut)[1] != nrow(suma)){
	     refine <- T
	   }
	   while(refine){
	     refine <- F
		 small.subc <- which(suma[,2] <= least.mut)[1]
		 if(small.subc == 1){
			target.ind <- 2					 
		 } else {
			distances <- abs(c(suma[small.subc,3]-suma[small.subc-1,3], suma[small.subc,3]-suma[small.subc+1,3]))
			target.ind <- small.subc+which.min(distances) * 2 - 3
		 }
		 
	     ind     <- c(which(label == suma[small.subc, 1]), which(label == suma[target.ind,1]))
	     new.phi <- (suma[small.subc, 2]*suma[small.subc, 3] + 
	                   suma[target.ind, 2]*suma[target.ind, 3]) / 
	       (suma[small.subc, 2] + suma[target.ind, 2])
	     new.label  <- suma[small.subc, 1]
	     phi[ind]   <- new.phi
	     label[ind] <- new.label
	     No.cls     <- No.cls - 1
	     suma <- matrix(NA, ncol = 3, nrow = No.cls)
	     suma[, 3] <- sort(unique(phi),decreasing = F)
	     suma[, 1] <- unlist(lapply(suma[, 3],
	                                function(y){
	                                  unique(label[which(phi == y)])
	                                }))
	     suma[, 2] <- unlist(lapply(suma[, 1],
	                                function(y){
	                                  length(which(label == y))
	                                }))
	     if(length(which(suma[,2]<= least.mut)) > 0 && No.cls > 2 && which(suma[,2]<= least.mut)[1] != nrow(suma) ){
	       refine <- T
	     }
	   }
	   write.table(suma, file = sprintf("%s/subclonal_structure.txt",
                   output.prefix), quote = F, sep = "\t", 
            col.names = F, row.names = F )
		write.table(label, file = sprintf("%s/mutation_assignments.txt",
		                                output.prefix), quote = F, sep = "\t", 
		           col.names = F, row.names = F )
				   
             
	}
}   

