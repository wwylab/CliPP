#args[1]: snv input file
#args[2]: cnv input file
#args[3]: purity file
#args[4]: output folder
#args[5]: mount point
args <- commandArgs(trailingOnly = TRUE)
setwd(args[5])
if(!dir.exists('./intermediate')){
	dir.create('./intermediate')
}
if(!dir.exists(args[4])){
	dir.create(args[4],recursive = TRUE)
}
if(!dir.exists('./semifinal')){
	dir.create('./semifinal')
}
#setwd('/CliP/')
cmd <- sprintf('python3.6 /CliP/preformat.py %s %s %s', args[1], args[2], args[3])
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- 'Rscript /CliP/preprocess.R ./intermediate/snv_input.txt ./intermediate/cnv_input.txt ./intermediate/purity.txt ./intermediate/middle/ ./intermediate/meta.Rdata'
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('python3.6 /CliP/run_CliP.py ./intermediate/middle/ ./intermediate/results1/ /CliP/ 0.1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('Rscript /CliP/post_analysis_run.R ./intermediate/results1/ ./intermediate/middle/ ./semifinal/res1/ 0.1 1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)


cmd <- sprintf('python3.6 /CliP/run_CliP.py ./intermediate/middle/ ./intermediate/results1/ /CliP/ 0.5')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('Rscript /CliP/post_analysis_run.R ./intermediate/results1/ ./intermediate/middle/ ./semifinal/res2/ 0.5 1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)


cmd <- sprintf('python3.6 /CliP/run_CliP.py ./intermediate/middle/ ./intermediate/results1/ /CliP/ 1.1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('Rscript /CliP/post_analysis_run.R ./intermediate/results1/ ./intermediate/middle/ ./semifinal/res3/ 1.1 1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)


cmd <- sprintf('python3.6 /CliP/run_CliP.py ./intermediate/middle/ ./intermediate/results1/ /CliP/ 1.5')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('Rscript /CliP/post_analysis_run.R ./intermediate/results1/ ./intermediate/middle/ ./semifinal/res4/ 1.5 1')
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

cmd <- sprintf('Rscript /CliP/output_results.R ./intermediate/purity.txt ./intermediate/snv_input.txt ./intermediate/middle/multiplicity.txt %s %s', args[4], args[5])
a <- system(cmd,intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE)

a <- system('rm -rf intermediate', intern=TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE) 

a <- system('rm -rf semiginal', intern=TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE) 

