#------------------------------------------------------------#
# This script takes care of the postprocessing for CliP.
# The postprocessing procedure mostly takes care of the harmonization of the results from the penalized model.
# Authors: Kaixian Yu, Yujie Jiang
# Email: yjiang13@mdanderson.org
#------------------------------------------------------------# 

#------------------------------------------------------------#
# Notation:
# CP: cellular prevalence
# CCF: cancer cell fraction, which is the proportion of tumor cells bearing the mutation
# Clonal mutation: mutations belonging to the initiating tumor cell and are expected to occur in every cell in the tumor
# Subclonal mutation: mutations which arose in descendant subpopulations
# Clonal fraction: number of clonal mutation divided by the number of total mutation
# Subclonal fraction: number of subclonal mutation divided by the number of total mutation
# Super cluster: a clone with CCF > 1, may indicate germline contamination or purity estimation errors
#------------------------------------------------------------#

#------------------------------------------------------------#
# Subjective variables:
# threshold: controls the lower limit of CP difference between any two clusters
# clonalFrac: controls the lower limit of clonal cluster in terms of percentage of mutations
# clonalFrac_superCluster: controls the lower limit of clonal super cluster in terms of percentage of mutations
# least.ratio: controls the lower limit of cluster size in terms of percentage of mutations 
#------------------------------------------------------------#



args <- commandArgs(trailingOnly = T)
library(parallel)
input.prefix        <- args[1]
preprocessed.prefix <- args[2]
output.prefix       <- args[3]
filtering           <- as.numeric(args[4])
threshold           <- 0.1
clonalFrac          <- 0.15
clonalFrac_superCluster <- 0.4
least.ratio         <- 0.05
SAFETY              <- 0.05
N.REP               <- 30
RATIO.TH            <- 0.90
# The 11 lambda values in the default lambda list
Lambda_list = c(0.01, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)

preliminary_files = list.files(input.prefix, pattern="*_label")

if (length(preliminary_files) == 0){
   stop("There is an error here.")
}
# else{
#    preliminary_files = stringi::stri_replace_all_fixed(preliminary_files, pattern = "_label.txt", replacement = "")
#    preliminary_files = stringi::stri_replace_all_fixed(preliminary_files, pattern = "lam", replacement = "")
#    Lambda_list = as.numeric(preliminary_files)
# }

if(!dir.exists(output.prefix)){
	dir.create(output.prefix)
}

# If user does not specify the lambda value 
if(length(args) < 5){
  # Loop over all 11 lambda values in the default lambda list
  for (lam in Lambda_list){
    print(lam)
    # Extract the output data
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
      
      # If doing subsampling
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
	  consensus.ccf <- lapply(pooled.beta,function(x){
            return(list(center = 0, size= 0, cluster = 0))
          })
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
    
    #### Start the filtering ####
    # The updated postprocessing pipeline includes 4 steps
    suma <- matrix(0,ncol=3, nrow=length(cluster.phi))
    for(i in 1:length(cluster.phi)){
      suma[i,1] = all.label[i]
      suma[i,2] = length(which(label==all.label[i]))
      suma[i,3] = cluster.phi[i]
    }
    if(filtering == 1){
      if(length(cluster.phi)==1){
        Mut_position <- read.table(sprintf("%s/multiplicity.txt", preprocessed.prefix), header=F)[c(1,2)]
        Mut_position <- cbind(Mut_position,label)
        colnames(Mut_position)[1]<-"chromosome_index"
        colnames(Mut_position)[2]<-"position"
        colnames(Mut_position)[3]<-"cluster_index"
        
        colnames(suma) <- c("cluster_index", "num_SNV", "cellular_prevalence")
        write.table(suma, file = sprintf("%s/subclonal_structure_lam%s.txt",
                                         output.prefix, lam), quote = F, sep = "\t", 
                    col.names = T, row.names = F )
        write.table(Mut_position, file = sprintf("%s/mutation_assignments_lam%s.txt",
                                          output.prefix, lam), quote = F, sep = "\t", 
                    col.names = T, row.names = F )
        
        cat(sprintf('There is only one clone in this sample, no further filtering is needed.'))
      }else{
        suma <- suma[order(suma[,3],decreasing = F),]
        least.mut <- least.ratio * length(phi)
        refine <- F
        No.cls <- dim(suma)[1]
        # Filter 1: deal with superclusters
        # If the max CP > 1 & the sample has > 2 clusters & the current clonal fraction <= 0.4:
        # Merge the two clusters with the largest CP values
        if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac_superCluster){
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
          if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac_superCluster){
            refine <- T
          }
        }
        
        # Filter 2: deal with small clones, i.e., the clonal cluster has a small number of mutations
        # If the sample has > 2 clusters & the current clonal fraction <= 0.15:
        # Merge the two clusters with the largest CP values
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
        # Filter 3: deal with adjacent clusters, i.e., cluster with similar CP values
        # If the sample has > 2 clusters & if CP values between any two clusters < 0.1
        # Merge those two clusters together
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
        # Filter 4: deal with small subclones, i.e., subclonal clusters with a small number of mutations
        # If a subclone has #SNV < 0.05 * total mutaton
        # Merge this subclone with its closest cluster
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
        
        rev_cluster_index <- rev(suma[,1]+1)
        index_dict <- list()
        for (i in c(1:(length(rev_cluster_index)))) {
          index_dict[[rev_cluster_index[i]]] <- i
        }
        
        #for (i in 1:length(label)){
        #  label[i] = index_dict[[label[i]]]
        #}
        
        label = sapply(label+1,
                       function(y){
                         y = index_dict[[y]] - 1
                       })
        
        Mut_position <- read.table(sprintf("%s/multiplicity.txt", preprocessed.prefix), header=F)[c(1,2)]
        Mut_position <- cbind(Mut_position,label)
        colnames(Mut_position)[1]<-"chromosome_index"
        colnames(Mut_position)[2]<-"position"
        colnames(Mut_position)[3]<-"cluster_index"
        
        colnames(suma) <- c("cluster_index", "num_SNV", "cellular_prevalence")
        
        suma[,"cluster_index"] = sort(0:(length(suma[,"cluster_index"])-1), decreasing = T)
        suma = suma[order(suma[,"cluster_index"]),]
        
        write.table(Mut_position, file = sprintf("%s/mutation_assignments_lam%s.txt",
                                          output.prefix, lam), quote = F, sep = "\t", 
                    col.names = T, row.names = F )
        write.table(suma, file = sprintf("%s/subclonal_structure_lam%s.txt",
                                         output.prefix, lam), quote = F, sep = "\t", 
                    col.names = T, row.names = F )
        
        
      }
    }   
  }
  
} else{ # If user specifies the lambda value 
  lam  <- args[5]
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
    
    # If doing subsampling
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
	consensus.ccf <- lapply(pooled.beta,function(x){
          return(list(center = 0, size= 0, cluster = 0))
        })
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
  
  #### Start the filtering ####
  # The updated postprocessing pipeline is in 4 steps
  suma <- matrix(0,ncol=3, nrow=length(cluster.phi))
  for(i in 1:length(cluster.phi)){
    suma[i,1] = all.label[i]
    suma[i,2] = length(which(label==all.label[i]))
    suma[i,3] = cluster.phi[i]
  }
  if(filtering == 1){
    if(length(cluster.phi)==1){
      Mut_position <- read.table(sprintf("%s/multiplicity.txt", preprocessed.prefix), header=F)[c(1,2)]
      Mut_position <- cbind(Mut_position,label)
      colnames(Mut_position)[1]<-"chromosome_index"
      colnames(Mut_position)[2]<-"position"
      colnames(Mut_position)[3]<-"cluster_index"
      
      colnames(suma) <- c("cluster_index", "num_SNV", "cellular_prevalence")
      write.table(suma, file = sprintf("%s/subclonal_structure_lam%s.txt",
                                       output.prefix, lam), quote = F, sep = "\t", 
                  col.names = T, row.names = F )
      write.table(Mut_position, file = sprintf("%s/mutation_assignments_lam%s.txt",
                                        output.prefix, lam), quote = F, sep = "\t", 
                  col.names = T, row.names = F )
      
      cat(sprintf('There is only one clone in this sample, no further filtering is needed.'))
    }else{
      suma <- suma[order(suma[,3],decreasing = F),]
      least.mut <- least.ratio * length(phi)
      refine <- F
      No.cls <- dim(suma)[1]
      # Filter 1: deal with superclusters
      # If the max CP > 1 & the sample has > 2 clusters & the current clonal fraction <= 0.4:
      # Merge the two clusters with the largest CP values
      if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac_superCluster){
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
        if(max(suma[, 3]) > 1 && No.cls > 2 && suma[No.cls,2] / sum(suma[,2]) < clonalFrac_superCluster){
          refine <- T
        }
      }
      
      # Filter 2: deal with small clones, i.e., the clonal cluster has a small number of mutations
      # If the sample has > 2 clusters & the current clonal fraction <= 0.15:
      # Merge the two clusters with the largest CP values
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
      # Filter 3: deal with adjacent clusters, i.e., cluster with similar CP values
      # If the sample has > 2 clusters & if CP values between any two clusters < 0.1
      # Merge those two clusters together
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
      # Filter 4: deal with small subclones, i.e., subclonal clusters with a small number of mutations
      # If a subclone has #SNV < 0.05 * total mutaton
      # Merge this subclone with its closest cluster
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

      rev_cluster_index <- rev(suma[,1]+1)
      index_dict <- list()
      for (i in c(1:(length(rev_cluster_index)))) {
        index_dict[[rev_cluster_index[i]]] <- i
      }


      label = sapply(label+1,
                     function(y){
                       y = index_dict[[y]] - 1
                     })
      
      Mut_position <- read.table(sprintf("%s/multiplicity.txt", preprocessed.prefix), header=F)[c(1,2)]
      Mut_position <- cbind(Mut_position,label)
      colnames(Mut_position)[1]<-"chromosome_index"
      colnames(Mut_position)[2]<-"position"
      colnames(Mut_position)[3]<-"cluster_index"
      
      colnames(suma) <- c("cluster_index", "num_SNV", "cellular_prevalence")
      
      suma[,"cluster_index"] = sort(0:(length(suma[,"cluster_index"])-1), decreasing = T)
      suma = suma[order(suma[,"cluster_index"]),]
      
      write.table(Mut_position, file = sprintf("%s/mutation_assignments_lam%s.txt",
                                        output.prefix, lam), quote = F, sep = "\t", 
                  col.names = T, row.names = F )
      write.table(suma, file = sprintf("%s/subclonal_structure_lam%s.txt",
                                       output.prefix, lam), quote = F, sep = "\t", 
                  col.names = T, row.names = F )   
      
    }
  }   
}

