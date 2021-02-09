#accept command line arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

#input is path and tf
setwd(path)
d = read.delim(paste(tf, "_thresholds.txt", sep=""), header=T)
fit = lm(-log10(d$P.value) ~ d$thres)
outfile = paste(tf, "_plinear", sep="")
write.table(fit$coefficients, file=outfile, quote=F, row.names=F, col.names=F, sep="\t")
