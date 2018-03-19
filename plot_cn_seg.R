#plot CNA segments
#Ruping Sun
#regularhand@gmail.com

#prepare cn data foa a sample
prepareCNAinput <- function(ff, subpre, subsuf, minnMark=40) {
    d = data.frame()
    for (i in 1:length(ff)) {
        name = gsub(subpre, "", ff[i])
        name = gsub(subsuf, "", name)
        cnvA = read.delim(ff[i])
        cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]
        cnvA$nt = cnvA$copynumber
        cnvA$nb = cnvA$minor_cn
        cnvA = cnvA[which(cnvA$num.mark > minnMark & cnvA$copynumber > 0),]
        
        cnv = data.frame(FID=paste(name,"Nt",sep="_"),CHR=cnvA$chrom,BP1=cnvA$loc.start,BP2=cnvA$loc.end,TYPE=cnvA$nt)
        cnv = rbind(cnv, data.frame(FID=paste(name,"Nb",sep="_"), CHR=cnvA$chrom,BP1=cnvA$loc.start,BP2=cnvA$loc.end,TYPE=cnvA$nb))
        if (i == 1) {
            d = cnv
        } else {
            d = rbind(d,cnv)
        }
    }
    return(d)
}


plotCNA <- function (d, ofile, w=12.5, h=2.35, layouth=c(2,3,1.3), cex.label=1.5, samples = vector(), cy_file="../OS/timing/cytoband_hg38.txt", sampleReOrder = FALSE) {

    #samples order
    if (length(samples) == 0) {
        samples = unique(d[,1])
    }
    ns = length(samples)
    s.sort = 1:ns
    if (sampleReOrder) {
        sorder = order(as.numeric(sapply(samples, function(x){sum(d$BP2[which(d$FID == x)] - d$BP1[which(d$FID == x)])})),
            decreasing=T)
        samples = samples[sorder]
    }
    ind = (d[,2] == "X" | d[,2] == "Y")
    d = d[!ind,]
    d[,2] = as.numeric(sub("chr","",d[,2]))     

    ## Cytoband stuff and Chromosome boundaries
    turn = 0
    cytoband=read.delim(cy_file,header=F)
    ind = (cytoband[,1] == "chrX" | cytoband[,1] == "chrY" | cytoband[,1] == "chrM")
    cytoband = cytoband[!ind,]
    cytoband[,1] = as.numeric(sub("chr","",cytoband[,1]))      #how to change chromosome name
    
    nchr=22
    chrLim=matrix(0,nchr,2)      #chromosome limits
    centromere=numeric()         #centromere limits
    cum=0
    for(i in 1:nchr)
        {
            ind=cytoband[,1]==i
            indP=regexpr("p",cytoband[ind,4])>0
            indQ=regexpr("q",cytoband[ind,4])>0
            
            chrLim[i,]=cbind(min(cytoband[ind,2]),max(cytoband[ind,3]))
            if(nchr==22)
                centromere=c(centromere, max(cytoband[ind,3][indP]))
            cum=c(cum,cum[i]+chrLim[i,2])
        }
    
    pdf(file=ofile, width = w, height=h)
    par(mar=c(0,0,0,0), cex.lab=3)
    par(oma=c(0,0,0,0))
    layout(matrix(c(1,2,3),3,1),height=layouth)

    #plot genome annotation                                                                         #200000000
    plot(NULL,NULL,axes=FALSE,xlim=c(cum[1]-200000000,cum[nchr+1]),ylim=c(-1.1,4),xlab="",ylab="")  #150000000
    sig=-1
    for(i in 1:nchr)
        {
            off=cum[i]
            if(sig == -1)
                ccol = "black"
            else
                ccol = "white"
            rect(cum[i],0,cum[i+1],-1,col=ccol)  #lwd=3
            if(nchr==22)
                lines(c(off+centromere[i],off+centromere[i]),c(0,-1),col="darkgray")  #lwd=2
            if(turn == 0)
                text(off+mean(chrLim[i,]),1.3-sig*0.6,paste(i),cex=2.5,font=2.5)
            if(turn == 1)
                text(off+mean(chrLim[i,]),1.3+sig*0.6,paste(i),cex=2.5,font=2.5,srt=90)
            sig=sig*(-1)
        }

    ##plot copy number: all
    plot(NULL,NULL,axes=FALSE,xlim=c(cum[1]-200000000,cum[nchr+1]),ylim=c(0,1),xlab="",ylab="")   #150000000
    for(k in 1:ns)
        {
            dd = d[d[,1] == samples[s.sort[k]],]
            y.down = (ns-k)/ns
            y.up = (ns-k+1)/ns
            trackn = "Total CN"
            if (grepl("_Nt", samples[s.sort[k]])) {
                y.down = y.down - 1/(2*ns)
            } else if (grepl("_Nb", samples[s.sort[k]])) {
                y.up = y.up - 1/(2*ns)
                trackn = "Minor CN"
            }

            #text(cum[1], (y.down+y.up)/2, labels=trackn, pos=2, cex=cex.label, font=2)
            if (k %% 2 == 1) {
                printSN = gsub("_Nt|_Nb", "", samples[s.sort[k]])
                printSN = gsub("^OS", "", printSN)
                text(cum[1]-20000000, (ns-k+1)/ns-1/(2*ns), labels=printSN, pos=2, cex=cex.label, font=2)
            }
            
            for(i in 1:nchr) {
                segments(cum[i], y.down, cum[i], y.up, lty=2, col=rgb(10/255,10/255,10/255,1/2))   #draw lattice #lwd=2
                ind = (dd[,2] == i)
                d.chr = dd[ind,]
                if(dim(d.chr)[1]!=0)
                    {
                        for(j in 1:length(d.chr[,2]))
                            {
                                currentCol = colors.cn(d.chr[j,5])
                                if (grepl("Nb",samples[s.sort[k]])){
                                    currentCol = colors.cn.haploid(d.chr[j,5])
                                }
                                rect(cum[i]+d.chr[j,3], y.down, cum[i]+d.chr[j,4], y.up, col = currentCol, border = NA, lwd =0)
                            }
                    }
            }
            if ( k %% 2 == 1 & k != 1 ) {
                segments(0, y.up, cum[length(cum)], y.up, lty=1, col=rgb(10/255,10/255,10/255,1/2)) #draw lattice for samples
            }
        }
    rect(cum[1],0,cum[nchr+1],1,border=T)  #lwd=3


    
    par(mar=c(2.5,0,0.5,0),cex.lab=3)
    plot(NULL,NULL,axes=FALSE, xlim=c(-3,43), ylim=c(0,1),xlab="",ylab="")
    off = 0
    for(i in seq(0.5,6,by=0.1)) {
        j = i+10
        rect(j,0,j+.1,1,col=colors.cn(i),border=NA,lwd=0)
    }
    segments(10.5,0,10.5,1)  #lwd=2
    segments(10.5,1,16,1)  #lwd=2
    segments(16,0,16,1)  #lwd=2
    segments(10.5,0,16,0)  #lwd=2
    axis(side=1, at=seq(11,16,by=1), labels=seq(1,6,by=1), cex.axis=1.6, font=2) #lwd=1.5
    text(6.5, 0.1, labels="Copy Number:", pos=2, cex=2.2, font=2)   
    text(10.5, 0.1, labels="Nt (total)", pos=2, cex=2.2, font=2)
    #additional one, may be removed
    off = 0
    for(i in seq(0,3,by=0.1)) {
        j = i+25
        rect(j,0,j+.1,1,col=colors.cn.haploid(i),border=NA,lwd=0)
    }
    segments(25,0,25,1)  #lwd=2
    segments(25,1,28,1) #lwd=2
    segments(28,0,28,1) #lwd=2
    segments(25,0,28,0)  #lwd=2
    axis(side=1, at=seq(25,28,by=1), labels=seq(0,3,by=1), cex.axis=1.6, font=2) #lwd=1.5
    text(25, 0.1, labels="Nb (minor allele)", pos=2, cex=2.2, font=2)

    dev.off()
}


colors.cn = function(cn,amp.start = 2.3,amp.end = 5,del.start = 1.7,del.end = 0.7)
  {
    if(cn > del.start  & cn < amp.start)
      return("white")
    else if(cn >= amp.start)
      {
        a = min(1,cn/(amp.end-amp.start)-amp.start/(amp.end-amp.start))
        return(hsv(0,alpha = a))
      }
        else if(cn <= del.start)
          {
            a = min(1,cn/(del.end-del.start)-del.start/(del.end-del.start))
            return(hsv(0.666666667,alpha = a))
          }
  }


colors.cn.haploid = function(cn,amp.start = 1.3,amp.end = 3,del.start = 0.7,del.end = 0)
    {
        if(cn > del.start  & cn < amp.start)
            return("white")
        else if(cn >= amp.start)
            {
                a = min(1,cn/(amp.end-amp.start)-amp.start/(amp.end-amp.start))
                return(hsv(0,alpha = a))
            }
        else if(cn <= del.start)
            {
                a = min(1,cn/(del.end-del.start)-del.start/(del.end-del.start))
                return(hsv(0.4908,0.8202,0.6980,alpha = a))
            }
    }
