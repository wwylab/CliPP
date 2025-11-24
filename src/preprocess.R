options(warn=-1)
# The script takes commandline arguments: input_SNV input_CNV purity_file sample_id Output_dir
args = commandArgs(trailingOnly=TRUE)
# debug use
snv.file      <- args[1]
cn.file       <- args[2]
purity.file   <- args[3]
sample.id     <- args[4]
output.prefix <- args[5]

print(snv.file)
print(cn.file)
print(purity.file)
print(sample.id)
print(output.prefix)

# if the input files do not exist then quit
if(!file.exists(snv.file)){
    stop(sprintf('The Input SNV file: %s Does not exist.', snv.file))
}
if(!file.exists(cn.file)){
    stop(sprintf('The Input CNV file: %s Does not exist.', cn.file))
}
if(!file.exists(purity.file)){
    stop(sprintf('The Input Purity file: %s Does not exist.', purity.file))
}


# if the output directory does not exist, then try to create the full path
if(!dir.exists(output.prefix)){
    dir.create(output.prefix, recursive = TRUE)
}
# utility function used to conduct the linear approximation

LinearApproximate <- function(bv, cv, cn, purity, diag.plot = FALSE) {
  # Local theta() using the parameters of this call
  theta <- function(w) {
    num <- exp(w) * bv
    den <- (1 + exp(w)) * (cn * (1 - purity) + cv * purity)
    num / den
  }
  
  # Simple local linear evaluator
  eval_line <- function(x, a, b) {
    a * x + b
  }
  
  # Fixed grid 
  w    <- -40:40 / 10
  No.w <- length(w)
  
  # True theta values on the grid
  actual.theta <- theta(w)
  
  # In this discrete setting, optimal cut is always -1.8 and 1.8
  i <- which.min(abs(w + 1.8))  # w[i] ≈ -1.8
  j <- which.min(abs(w - 1.8))  # w[j] ≈  1.8
  
  # Segment 1: [w[1], w[i]]
  a1 <- (actual.theta[i] - actual.theta[1]) / (w[i] - w[1])
  b1 <- actual.theta[1] - w[1] * a1
  
  # Segment 2: [w[i], w[j]]
  a2 <- (actual.theta[j] - actual.theta[i]) / (w[j] - w[i])
  b2 <- actual.theta[i] - w[i] * a2
  
  # Segment 3: [w[j], w[No.w]]
  a3 <- (actual.theta[No.w] - actual.theta[j]) / (w[No.w] - w[j])
  b3 <- actual.theta[No.w] - w[No.w] * a3
  
  # Piecewise approximation
  approx.theta <- c(
    eval_line(w[1:i],        a1, b1),
    eval_line(w[(i+1):j],    a2, b2),
    eval_line(w[(j+1):No.w], a3, b3)
  )
  
  # Max absolute error
  diff <- max(abs(actual.theta - approx.theta))
  
  # Optional diagnostic plot
  if (diag.plot) {
    plot(w, actual.theta, type = "l",
         ylim = c(0, 1), xlim = c(-5, 5),
         xlab = "w", ylab = "theta")
    lines(w, approx.theta, col = 2)
    abline(v = c(-1.8, 1.8), lty = 2)
  }
  
  list(
    w.cut = c(w[i], w[j]),           # should be -1.8, 1.8
    diff  = diff,
    coef  = c(a1, b1, a2, b2, a3, b3)
  )
}
###
CombineReasons <- function(chrom, pos, ind, reason){
    res <- NULL
    if(length(ind) > 0){
        for (i in (1:length(ind))){
            res[i] <- sprintf("%d\t%d\t%s",chrom[ind[i]],pos[ind[i]],reason)
        }
    }
    return(res)
}
###

VALID.CONT     <- 0
case.store     <- NULL # stores coefs
cutbeta.store  <- NULL
cuttheta.store <- NULL
valid.store    <- NULL # "cn_cv_bv"
invalid.store  <- NULL
dropped.SNV    <- NULL
# problem.sample <- NULL

pp.table  <- read.table(purity.file)
purity <- pp.table$V1[1]

tmp.vcf        <- read.table(snv.file, header=T, stringsAsFactors = F)
mutation.chrom <- as.numeric(tmp.vcf$chromosome_index)
mutation.pos   <- as.numeric(tmp.vcf$position)
# take off the sex chromosome
valid.ind      <- which(!is.na(mutation.chrom))
drop.ind       <- which(is.na(mutation.chrom))
dropped.SNV    <- CombineReasons(mutation.chrom, mutation.pos, drop.ind, "The SNV is on sex chromosomes.")
if(length(valid.ind) < VALID.CONT ){
    # problem.sample <- rbind(problem.sample,c(sample.id, sprintf("Only %d mutations with more than 2 callers and not sex", length(valid.ind))))
    stop(sprintf('The sample with SNV %s has less than %d SNVs that are on non-sex chromosomes.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- tmp.vcf$alt_count[valid.ind]
total.read     <- tmp.vcf$alt_count[valid.ind] + tmp.vcf$ref_count[valid.ind]
# take only non-negative counts
valid.ind      <- intersect(which(minor.read >= 0),which(total.read>=0))
drop.ind       <- setdiff(1:length(minor.read),valid.ind)
dropped.SNV    <- append(dropped.SNV,CombineReasons(mutation.chrom, mutation.pos, drop.ind, "The SNV has negative reads."))
if(length(valid.ind) < VALID.CONT ){
    # problem.sample <- rbind(problem.sample,c(sample.id,sprintf("Only %d mutations with non-negative reads",length(valid.ind)) ))
    stop(sprintf('The sample with SNV %s has less than %d SNVs that have non-negative reads.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
No.mutations   <- length(valid.ind)

# process copy number

cn.tmp         <- read.table(cn.file,header=T, stringsAsFactors = F)
cn.tmp         <- cn.tmp[which(!is.na(cn.tmp[,"minor_cn"])),]
No.cnLines     <- nrow(cn.tmp)
if(No.cnLines == 0){
    stop(sprintf('The sample with SNV %s does not have valid copy number status.',snv.file))
}
mut.cna.id     <- unlist(lapply(1:No.mutations, function(x){
    ret.val <- -1
    for(i in 1:No.cnLines){
        if( mutation.chrom[x] == cn.tmp[i, "chromosome_index"]
            && mutation.pos[x] >= cn.tmp[i, "start_position"]
            && mutation.pos[x] <= cn.tmp[i, "end_position"]){
            ret.val <- i
            break
        }
    }
    return(ret.val)
}))
valid.ind      <- which(mut.cna.id > 0)
drop.ind       <- setdiff(1:length(minor.read),valid.ind)
dropped.SNV    <- append(dropped.SNV,CombineReasons(mutation.chrom, mutation.pos, drop.ind, "The SNV does not have valid copy number."))
if(length(valid.ind) < VALID.CONT ){
    stop(sprintf('The sample with SNV %s has less than %d SNVs that have valid copy number status.',snv.file,VALID.CONT))
}
mutation.chrom <- mutation.chrom[valid.ind]
mutation.pos   <- mutation.pos[valid.ind]
minor.read     <- minor.read[valid.ind]
total.read     <- total.read[valid.ind]
No.mutations   <- length(valid.ind)
mut.cna.id     <- mut.cna.id[valid.ind]
minor.copy.lim <- apply(cn.tmp[mut.cna.id,c("minor_cn","major_cn")],1,max)
total.count    <- cn.tmp[mut.cna.id,"total_cn"]
# calculate multiplicity
multiplicity   <- round(minor.read/total.read/purity*(total.count*purity+(1-purity)*2))
# multiplicity should not exceed larger copy number of the two allels
minor.count    <- apply(cbind(minor.copy.lim, multiplicity), 1, min)
minor.count[minor.count == 0] <- 1
valid.ind      <- intersect(which(minor.count>0),which(total.count>0))
drop.ind       <- setdiff(1:length(minor.read),valid.ind)
dropped.SNV    <- append(dropped.SNV,CombineReasons(mutation.chrom, mutation.pos, drop.ind, "The SNV has negative multiplicities."))

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
        res <- LinearApproximate(minor.count[mutation],total.count[mutation],2, purity, diag.plot = FALSE)
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
drop.ind       <- setdiff(1:length(minor.read),valid.ind)
dropped.SNV    <- append(dropped.SNV,CombineReasons(mutation.chrom, mutation.pos, drop.ind, 
                                                    "The copy numbers for the SNV is not stable enough to calculate the approximated line."))
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
#sample.cuttheta <- sample.cuttheta[valid.ind,]
No.mutations   <- length(valid.ind)

phi <- 2/(minor.count/(minor.read/total.read) - total.count + 2)
valid.ind <- intersect(which(phi <= 1.5), which(phi > 0 ))
clonal.ind <- which(phi > 1.5)
if(length(clonal.ind) > 0){
    outlier.higherEnd <- cbind(mutation.chrom[clonal.ind],mutation.pos[clonal.ind],
                               total.count[clonal.ind], minor.count[clonal.ind])
    write.table(outlier.higherEnd, file=sprintf("%s/outPosition.txt", output.prefix),
                quote=F, col.names = F, row.names = F)
}


drop.ind       <- setdiff(1:length(minor.read),valid.ind)
dropped.SNV    <- append(dropped.SNV,CombineReasons(mutation.chrom, mutation.pos, drop.ind, 
                                                    "The empirical CP is off the chart, which may be caused by incorrect copy number or existence of super cluster(s)"))
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

output.dropped <- sprintf("%s/excluded_SNVs.txt",output.prefix)
write.table(minor.read, output.r,  quote=F, col.names = F, row.names = F)
write.table(total.read, output.n,  quote=F, col.names = F, row.names = F)
write.table(minor.count, output.minor,  quote=F, col.names = F, row.names = F)
write.table(total.count, output.total,  quote=F, col.names = F, row.names = F)
write.table(index, output.index,  quote=F, col.names = F, row.names = F)
write.table(purity, output.pp,  quote=F, col.names = F, row.names = F)
write.table(sample.coef, output.coef,  quote=F, col.names = F, row.names = F, sep="\t")
write.table(sample.cutbeta, output.cutbeta,  quote=F, col.names = F, row.names = F)
write.table(dropped.SNV, output.dropped,  quote=F, col.names = F, row.names = F)
# save(case.store, cutbeta.store, cuttheta.store, valid.store, invalid.store, sample.id, problem.sample, file="approximate_meta_data.Rdata")
