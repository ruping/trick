#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for calculating the LOD of cn segments between a pair of samples

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 6) stop("Wrong number of input parameters: 'sample1 sample2 sn1 sn2 cnvres path'")


sample1 <- inputpar[1]
sample2 <- inputpar[2]
sn1     <- inputpar[3]
sn2     <- inputpar[4]
cnvres  <- inputpar[5]
path    <- inputpar[6]

library(matrixStats)
setwd(path)


cnlod <- function(tt1, tt2) {
    #a = b
    rowsum1 = rowSums(log(tt1[[1]]))+rowSums(log(tt1[[2]]))
    rowsum2 = rowSums(log(tt2[[1]]))+rowSums(log(tt2[[2]]))
    rowsumt = rowsum1 + rowsum2
    anyvalue = max(rowsumt)
    logEQ = anyvalue + log(sum(exp(rowsumt - anyvalue)))
    
    #a != b
    combs = combn(25,2)
    diffsumt = as.vector(apply(combs, 2, function(x, rowsum1, rowsum2){
              c(rowsum1[x[1]] + rowsum2[x[2]], rowsum2[x[1]] + rowsum1[x[2]])
          }, rowsum1=rowsum1, rowsum2=rowsum2))
    anyvalue = max(diffsumt)
    logDF = anyvalue + log(sum(exp(diffsumt - anyvalue)))

    #compute lod
    lod = logDF - logEQ
    return(lod)
}

cnlrbafprob <- function(data, params,
                        major_cn_code = c(0,1,2,1,3,2,4,3,2,5,4,3,6,5,4,3,7,6,5,4,8,7,6,5,4),
                        minor_cn_code = c(0,0,0,1,0,1,0,1,2,0,1,2,0,1,2,3,0,1,2,3,0,1,2,3,4)) {

    norcon = tail(params$n,1)
    #message(paste("normal contamination is: ", norcon, sep=""))
    ploidy = tail(params$phi,1)
    #message(paste("ploidy is: ", ploidy, sep=""))
    mu = params$muC[1:25,1,5]
    var = params$var[,5]
    lr = log(2^as.numeric(data$LogRatio))
    #message(paste(lr, collapse="\t"))
    nprob = sapply(lr, function(x) {pnormTwoTails(rep(x,25), mu, sqrt(var))})
    
    prev1= sapply(as.numeric(data$CellularPrevalence),function(x){
                      if (is.na(x)){1} else {x}
                      })
    N = as.numeric(data$RefCount) + as.numeric(data$NRefCount)
    k = pmax(as.numeric(data$RefCount), as.numeric(data$NRefCount))
    bprob = sapply(1:length(k), function(x, k, N, prev) {
                       omega = balleleRatio(n=rep(norcon, 25), prev=rep(prev1[x],25),
                           Ct=major_cn_code+minor_cn_code, Cb=minor_cn_code)
                       pbinomTwoTails(rep(k[x],25), rep(N[x],25), omega)
                   }, k=k, N=N, prev=prev1)
    return(list(nprob,bprob))
}

pnormTwoTails <- function(q, mean, sd) {
    p = sapply(1:length(q), function(x, q, mean, sd){
                   if (q[x] < mean[x]) {
                       2*pnorm(q[x], mean[x], sd[x])
                   } else {
                       2*pnorm(q[x], mean[x], sd[x], lower.tail = FALSE)
                   }
               }, q=q, mean=mean, sd=sd)
    p = as.vector(p)
    return(p)
}

pbinomTwoTails <- function(q, size, prob) {
    p = sapply(1:length(q), function(x, q, size, prob) {
                 binom.test(q[x], size[x], prob[x], alternative="two.sided")$p.value
               }, q=q, size=size, prob=prob)
    p = as.vector(p)
    return(p)
}

logRatio <- function(n, phi, Ct, prev, Cn=2) {
    lr = (n*Cn + (1-n)*(1-prev)*Cn + (1-n)*prev*Ct)/(n*Cn + (1-n)*phi)
    return(lr)
}

balleleRatio <- function(n, rn=0.5, Cn=2, prev, Ct, Cb) {
    numRefAlleles <- n*rn*Cn + (1-n)*(1-prev)*rn*Cn + (1-n)*prev*(Ct-Cb)
    totalAlleles <- n*Cn + (1-n)*(1-prev)*Cn + (1-n)*prev*Ct
    br = numRefAlleles/totalAlleles
    return(br)
}


#start processing
cnvA2 = read.delim(cnvres)
load(sample1)
sample1.titanres = titancnaresults
load(sample2)
sample2.titanres = titancnaresults
cnlodres = vector()
nprobes = vector()
for (i in 1:dim(cnvA2)[1]) {
    chr = cnvA2[i,"chrom"]
    start = as.numeric(cnvA2[i,"loc.start"])
    end = as.numeric(cnvA2[i, "loc.end"])
    cn1 = cnvA2[i, sn1]
    cn2 = cnvA2[i, sn2]
    data = sample1.titanres[[1]]$results[with(sample1.titanres[[1]]$results, Chr == chr &
                                           as.numeric(Position) >= start & as.numeric(Position) <= end),]
    data2 = sample2.titanres[[1]]$results[with(sample2.titanres[[1]]$results, Chr == chr &
                                            as.numeric(Position) >= start & as.numeric(Position) <= end),]
    allpos = intersect(data$Position[which(data$TITANcall == names(which.max(table(data$TITANcall))))],
        data2$Position[which(data2$TITANcall == names(which.max(table(data2$TITANcall))))])
    npos = length(allpos)
    data = data[match(allpos, data$Position),]
    data2 = data2[match(allpos, data2$Position),]
    
    cnlodc = 0
    if (npos >= 2) {
        tt1 = cnlrbafprob(data, sample1.titanres[[1]]$convergeParams)
        tt2 = cnlrbafprob(data2, sample2.titanres[[1]]$convergeParams)
        cnlodc = cnlod(tt1, tt2)/npos
    }
    message(paste(c(i,chr,start,end,cn1,cn2,npos,cnlodc), collapse="\t"))
    cnlodres = append(cnlodres,  cnlodc)
    nprobes = append(nprobes, npos)
}

cnvA2 = data.frame(cnvA2, nprobes, cnlodres)
colnames(cnvA2)[(dim(cnvA2)[2]-1):dim(cnvA2)[2]] = c("npos", paste("cnlod",sn1,sn2,sep="."))
write.table(cnvA2, file=paste("cnlod",sn1,sn2,"tsv",sep="."), quote=F, sep="\t", row.names=F)
