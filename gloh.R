#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for germline LOH detection

inputpar <- commandArgs(TRUE)
if (length(inputpar) < ) stop("Wrong number of input parameters: 'path sampleName numMark'")

path <- inputpar[1]
sampleName <- inputpar[2]
numMark <- inputpar[3]

options(width=130)
library(DNAcopy)
library(HMMcopy)


#LOH
#list.files("./","gloh",full=T)
file = paste(path, "/", sampleName, ".gloh", sep="")
d = read.delim(file)

d = d[which(d$chr != "X" & d$chr != "Y" & d$chr != "MT"),]
chrs = names(table(d$chr))
mvalue = vector()
for (i in 1:22) {
  message(chrs[i])
  mvalue = append(mvalue, runmed(d$value[which(d$chr == chrs[i])], 11))
}

loh = data.frame(chrom=d$chr, maploc=d$pos, mvalue=mvalue)
pdf(file=paste(path,"/",sampleName,".pdf",sep=""), width=16,height=10)
par(mfrow=c(4,6))
for (i in 1:22) {
    plot(loh$mvalue[which(loh$chrom == i)], ylim=c(0,1), main=paste("chr", i, sep=""),
         col=as.vector(sapply(loh$mvalue[which(loh$chrom == i)],function(x){if (x == 0){"blue"} else {"red"}})),
         axes=F,ylab="",xlab="SNP coor index",cex.lab=1.5,cex.main=2)
    #smoothScatter(loh$maploc[which(loh$chrom == i)],loh$mvalue[which(loh$chrom == i)], ylim=c(0,1), main=paste("chr", i, sep=""),
         #col=as.vector(sapply(loh$mvalue[which(loh$chrom == i)],function(x){if (x == 0){"blue"} else {"red"}})),
    #     axes=F,ylab="",xlab="SNP coor index",cex.lab=1.5,cex.main=2)
    axis(side=2, at=c(0,1), labels=c("hetero","homo"),cex.axis=1.5)
}
dev.off()


#par(mfrow=c(1,1))
sloh = CNA(loh$mvalue,loh$chrom,loh$maploc,
   data.type="binary",sampleid=sampleName)
segloh = segment(sloh)
#plot(density(log2(as.numeric(segloh$output$num.mark))))
#plot(segloh, plot.type="w")

seg = segloh$output
seg = seg[which(seg$num.mark > numMark),]
write.table(seg, file=paste(path,"/",sampleName,"gloh.seg",sep=""), quote=F, sep="\t", row.names=F)
