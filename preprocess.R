#------------------------------------------------------------#
# This script takes care of the preprocessing for CliP
# Generating the input file for CliP and de-identify all the samples
# The input files for each sample should be in one directory
# Different samples should be in different folder
# Initialized by Kaixian Yu
# Date: 05/04/2018
# Email: kaixiany@163.com
#------------------------------------------------------------#
# The script takes commandline arguments: input_SNV input_CNV purity_file Output_dir meta_file
args = commandArgs(trailingOnly=TRUE)
# debug use
# args <- c('/Users/kaixiany/Working/CliP/Sample_data/sample_SNV.txt', '/Users/kaixiany/Working/CliP/Sample_data/sample_cnv_ver1.txt', '/Users/kaixiany/Working/CliP/Sample_data/sample_purity.txt', '/Users/kaixiany/Working/CliP/Sample_data/intermediate/', '/Users/kaixiany/Working/CliP/approximate_meta_data.Rdata')
snv.file      <- args[1]
cnv.file      <- args[2]
purity.file   <- args[3]
output.prefix <- args[4]
meta.file     <- args[5]
# if the input files do not exist then quit
if(!file.exists(snv.file)){
	stop(sprintf('The Input SNV file: %s Does not exist.',snv.file))
}
if(!file.exists(cnv.file)){
	stop(sprintf('The Input CNV file: %s Does not exist.',snv.file))
}
if(!file.exists(purity.file)){
	stop(sprintf('The Input Purity file: %s Does not exist.',purity.file))
}
# if the output directory does not exist, then try to create the full path
if(!dir.exists(output.prefix)){
	dir.create(output.prefix, recursive = TRUE)
}
# some utility functions used to conduct the linear approximation

theta = function(w,bv,cv,cn){
  #return((exp(w)*bv)/(cn+exp(w)*cv))
  return((exp(w)*bv)/( (1+exp(w))*cn*(1-purity) + (1+exp(w))*cv*purity) )
}
LinearEvaluate <- function(x,a,b){
  return(a*x+b)
}
LinearApproximate <- function(bv,cv,cn,diag.plot = F){
  w            <- -40:40/10
  actual.theta <- theta(w,bv,cv,cn) 
  No.w         <- length(w) 
  diff         <- rep(0,(No.w-3)*(No.w-2)/2)
  coef         <- matrix(0,ncol=6,nrow = (No.w-3)*(No.w-2)/2)
  w.cut        <- matrix(0,ncol=2,nrow = (No.w-3)*(No.w-2)/2)
  k            <- 1
  approximated.theta <- list()
  for( i in 2:(length(w)-2)){
    for(j in (i+1):(length(w)-1) ){
      w.cut[k,] <- c(w[i],w[j])
      coef[k,1] <- (actual.theta[i]-actual.theta[1])/(w[i]-w[1])
      coef[k,2] <- actual.theta[1]-w[1]*coef[k,1]
      coef[k,3] <- (actual.theta[j]-actual.theta[i])/(w[j]-w[i])
      coef[k,4] <- actual.theta[i]-w[i]*coef[k,3]
      coef[k,5] <- (actual.theta[No.w]-actual.theta[j])/(w[No.w]-w[j])
      coef[k,6] <- actual.theta[No.w]-w[No.w]*coef[k,5]
      approximated.theta[[k]] <- c(LinearEvaluate(w[1:i],coef[k,1],coef[k,2]),
                                   LinearEvaluate(w[(i+1):j],coef[k,3],coef[k,4]),
                                   LinearEvaluate(w[(j+1):No.w],coef[k,5],coef[k,6]))
      diff[k]  <- max(abs(actual.theta-approximated.theta[[k]]))
      k        <- k+1
    }
  }
  K <- which.min(diff)
  if(diag.plot){
    plot(theta(w,bv,cv,cn)~w, ylim=c(0,1),xlim=c(-5,5),type="l")
    par(new = T)
    plot(approximated.theta[[K]]~w,ylim=c(0,1),xlim=c(-5,5),type="l",col=2,ylab = "",xlab ="",main="")
  }
  
  return(list(w.cut = w.cut[K,],diff= diff[K], coef=coef[K,] ))
}
####---------------------------------------------------------------------####
VALID.CONT     <- 0
case.store     <- NULL # stores coefs
cutbeta.store  <- NULL
cuttheta.store <- NULL
valid.store    <- NULL # "cn_cv_bv"
invalid.store  <- NULL

if(file.exists(meta.file)){
        load(meta.file)
}

purity         <- unlist(read.table(purity.file))
tmp.vcf        <- read.table(snv.file, sep='\t', stringsAsFactors = FALSE)
mutation.chrom <- as.numeric(tmp.vcf[, 1])
mutation.pos   <- as.numeric(tmp.vcf[, 2])
valid.ind      <- which(!is.na(mutation.chrom))
if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that are on non-sex chromosomes.',snv.file,VALID.CONT))
}

mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
tmp.vcf        <- tmp.vcf[valid.ind,]
minor.read     <- tmp.vcf[,3]
total.read     <- tmp.vcf[,3] + tmp.vcf[,4]
valid.ind      <- intersect(which(minor.read >= 0),which(total.read>=0))
if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that have non-negative reads.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
No.mutations   <- length(valid.ind)

cn.tmp         <- read.table(cnv.file, sep = '\t', stringsAsFactors = FALSE)
No.cnLines     <- nrow(cn.tmp)
if(No.cnLines == 0){
    stop(sprintf('The sample with SNV %s does not have valid copy number status.',snv.file))
}
mut.cna.id     <- unlist( lapply(1:No.mutations, 
						function(x){
                        ret.val <- -1
                        for(i in 1:No.cnLines){
                          if( mutation.chrom[x] == cn.tmp[i, 1] 
                              && mutation.pos[x] >= cn.tmp[i, 2] 
                              && mutation.pos[x] <= cn.tmp[i, 3] ){
                            ret.val <- i
                            break
                          }
                          
                        }
                        return(ret.val)
                      } ) )
valid.ind      <- which(mut.cna.id > 0)
if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that have valid copy number status.',snv.file,VALID.CONT))
}

mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
No.mutations   <- length(valid.ind)
mut.cna.id     <- mut.cna.id[valid.ind]

if(ncol(cn.tmp) == 5){
	minor.copy.lim <- apply(cn.tmp[mut.cna.id,4:5],1,max)
	total.count    <- cn.tmp[mut.cna.id,4] + cn.tmp[mut.cna.id,5]
} else if(ncol(cn.tmp) == 4){
	minor.copy.lim <- cn.tmp[mut.cna.id,4]
	total.count    <- cn.tmp[mut.cna.id,4]
}

# calculate multiplicity
multiplicity   <- round(minor.read/total.read/purity*(total.count*purity+(1-purity)*2))
# multiplicity should not exceed larger copy number of the two allels
minor.count    <- apply(cbind(minor.copy.lim, multiplicity), 1, min)
minor.count[minor.count == 0] <- 1
valid.ind <- intersect(which(minor.count>0),which(total.count>0))
if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that have valid copy number status for SNVs.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
No.mutations   <- length(valid.ind)
total.count    <- total.count[valid.ind]
minor.count    <- minor.count[valid.ind]

sample.coef    <- matrix(0,nrow = No.mutations, ncol = 6)
sample.cutbeta <- matrix(0,nrow = No.mutations, ncol = 2)
sample.diff    <- rep(0,No.mutations)
for(mutation in 1:No.mutations){
	approx.ind <- which(valid.store == paste(c(2,total.count[mutation],minor.count[mutation]),collapse = "_"))
	if(length(approx.ind)>0){
	  sample.coef[mutation,] <- case.store[approx.ind,]
	  sample.cutbeta[mutation,] <- cutbeta.store[approx.ind,]
	} else if(length(which(invalid.store == paste(c(2,total.count[mutation],minor.count[mutation]),collapse = "_") )) != 0 ) {
	  sample.diff[mutation] <- 1
	} else {
	  res <- LinearApproximate(minor.count[mutation],total.count[mutation],2)
	      if(res$diff <= 0.1 ){
	            sample.coef[mutation,]    <- res$coef
	    sample.cutbeta[mutation,]  <- res$w.cut
	    sample.diff[mutation] <- res$diff
	    valid.store <- c(valid.store,paste(c(2,total.count[mutation],minor.count[mutation]),collapse = "_"))
	    case.store <- rbind(case.store,res$coef)
	    cutbeta.store <- rbind(cutbeta.store,res$w.cut)
	    #cuttheta.store <- c(cuttheta.store,res$cut.theta)
	  } else {
	    sample.diff[mutation] <- 1
	    invalid.store <- c(invalid.store, paste(c(2,total.count[mutation],minor.count[mutation]),collapse = "_"))
	  }
	}
}
valid.ind <- which(sample.diff <= 0.1)
if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that have valid approximated theta.',snv.file,VALID.CONT))
}

mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
minor.count    <- minor.count[valid.ind]
total.count    <- total.count[valid.ind]
sample.coef    <- sample.coef[valid.ind,]
sample.cutbeta <- sample.cutbeta[valid.ind,]
No.mutations   <- length(valid.ind)

phi <- 2/(minor.count/(minor.read/total.read) - total.count + 2)
valid.ind <- intersect(which(phi <= 1.5), which(phi > 0 ))
clonal.ind <- which(phi > 1.5)
if(length(clonal.ind) > 0){
	outlier.higherEnd <- cbind(mutation.chrom[clonal.ind], mutation.pos[clonal.ind], total.count[clonal.ind], minor.count[clonal.ind])
	write.table(outlier.higherEnd, file=sprintf("%s/outPosition.txt", output.prefix), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

png(file=sprintf("%s/phi.png", output.prefix))
hist(phi[valid.ind])
dev.off()


if(length(valid.ind) < VALID.CONT ){
	stop(sprintf('The sample with SNV %s has less than %d SNVs that have valid empirical phi.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
minor.count    <- minor.count[valid.ind]
total.count    <- total.count[valid.ind]
sample.coef    <- sample.coef[valid.ind,]
sample.cutbeta <- sample.cutbeta[valid.ind,]
No.mutations   <- length(valid.ind)

index          <- cbind(mutation.chrom,mutation.pos,total.count, minor.count)
# output preprocessed results
output.r       <- sprintf("%s/r.txt",output.prefix)
output.n       <- sprintf("%s/n.txt",output.prefix)
output.minor   <- sprintf("%s/minor.txt",output.prefix)
output.total   <- sprintf("%s/total.txt",output.prefix)
output.index   <- sprintf("%s/multiplicity.txt",output.prefix)
output.pp      <- sprintf("%s/purity_ploidy.txt",output.prefix)
output.coef    <- sprintf("%s/coef.txt",output.prefix)
output.cutbeta <- sprintf("%s/cutbeta.txt",output.prefix)

write.table(minor.read, output.r,  quote=F, col.names = F, row.names = F)
write.table(total.read, output.n,  quote=F, col.names = F, row.names = F)
write.table(minor.count, output.minor,  quote=F, col.names = F, row.names = F)
write.table(total.count, output.total,  quote=F, col.names = F, row.names = F)
write.table(index, output.index,  quote=F, col.names = F, row.names = F)
write.table(purity, output.pp,  quote=F, col.names = F, row.names = F)
write.table(sample.coef, output.coef,  quote=F, col.names = F, row.names = F, sep="\t")
write.table(sample.cutbeta, output.cutbeta,  quote=F, col.names = F, row.names = F)
save(case.store, cutbeta.store, cuttheta.store, valid.store, invalid.store, file = meta.file)

















