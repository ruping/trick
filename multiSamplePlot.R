# plot multi sample mutation Table
script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir, "smkey.R", sep="/"))

getSampMutMulti <- function(samples, normal, d, cmedianTh, original) {
    rindexTF = vector()
    cindex = vector()
    for (si in 1:length(samples)) {
        if (length(rindexTF) > 0) {
            rindexTF = rindexTF | grepl(paste(samples[si],"\\[",sep=""), d$somatic)
        } else {
            rindexTF = grepl(paste(samples[si],"\\[",sep=""), d$somatic)
        }
        cindex = append(cindex, match(paste(samples[si], "maf", sep=""), colnames(d)))
    }
    rindex = which(rindexTF)
    cindex = as.vector(sapply(cindex, function(x){c(x-1,x,x+1)}))
    
    ncindex = match(paste(normal, "maf", sep=""), colnames(d))
    ncindex = as.vector(sapply(ncindex, function(x){c(x,x+1)}))  #normal samples maf and depth indexes
    
    dronindex = match("dron", colnames(d))
    gnindex = match("geneName", colnames(d))
    glindex = match("geneLoc", colnames(d))
    gfindex = match("functionalClass", colnames(d))
    caddindex = match("CADD_phred", colnames(d))
    somindex = match("somatic", colnames(d))

    if (! is.na(caddindex)) {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,caddindex,dronindex,somindex,ncindex)]
    } else {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,dronindex,somindex,ncindex)]
    }
    
    resColnames = colnames(res)
    resAdd = t(data.frame(apply(res, 1, function(x, cindex, original, resColnames) {     #x is every row
              maxMaf = 0
              maxTlod = 0
              ssb = 0                    #combined strand bias check
              ssbc = 0
              for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample get maxMaf
                  ss = strsplit(as.character(x[j+1]), "\\|")
                  mafTmp = as.numeric(ss[[1]][1])
                  oriTlodTmp = 0
                  if ( grepl("\\|", as.character(x[j])) ) {
                      ori = strsplit(as.character(x[j]), "\\|")
                      oriLod = strsplit(ori[[1]][3], ",")            
                      oriTlodTmp = as.numeric(oriLod[[1]][1])
                  }
                  if (mafTmp > maxMaf) {
                      maxMaf = mafTmp    #define max maf
                  }
                  if (original > 0 & oriTlodTmp > maxTlod) {
                      maxTlod = oriTlodTmp
                  }
              }
              resVector = vector()
              for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample
                  sampleName = resColnames[j]
                  sampleNameMaf = paste(sampleName, "mafc", sep="")
                  sampleNameMafa = paste(sampleName, "mafa", sep="")
                  sampleNameCcf = paste(sampleName, "ccf", sep="")
                  sampleNameCcfSD = paste(sampleName, "ccfSD", sep="")
                  sampleNameAlt = paste(sampleName, "altc", sep="")
                  sampleNameRef = paste(sampleName, "refc", sep="")
                  oriTlod = 0
                  if ( grepl("\\|", as.character(x[j])) ) {
                      ori = strsplit(as.character(x[j]), "\\|")
                      oriLod = strsplit(ori[[1]][3], ",")
                      oriTlod = as.numeric(oriLod[[1]][1])
                  }
                  ss = strsplit(as.character(x[j+1]), "\\|")
                  mafTmp = as.numeric(ss[[1]][1])
                  cmeme = strsplit(ss[[1]][3], ",")
                  endBias = as.numeric(ss[[1]][2])
                  strandBias = as.numeric(ss[[1]][4])
                  mappingBias = as.numeric(ss[[1]][5])
                  cmedianav = as.numeric(cmeme[[1]][2])
                  cmemeSum = sum(as.numeric(cmeme[[1]]))

                  
                  if ( mafTmp > 0 ) {
                      ssb = ssb+strandBias
                      ssbc = ssbc + 1
                  }
                  
                  #decide a b allele count
                  refnow = round(as.numeric(x[j+2])*(1-mafTmp))
                  altnow = round(as.numeric(x[j+2])*mafTmp)

                  #decide mafNow
                  if ( mafTmp == 0 ) {
                      mafNow = mafTmp
                  } else if (endBias < 0.9 & strandBias != 1 & strandBias != 0 & mappingBias < 0.8 & cmemeSum < 5.2 & cmedianav < cmedianTh) {
                  #} else if (endBias < 0.9 & mappingBias < 0.99 & cmemeSum < 7.5 & cmedianav < cmedianTh) {    
                      mafNow = mafTmp
                      if (original > 0) {   #if oriLod
                          if (oriTlod < original) {
                              if (!(maxMaf > 0.1 & (original == 0 | (original > 0 & maxTlod >= original)))) {   #not called
                                  mafNow = 0
                              }
                          }
                      }
                  } else {
                      if (maxMaf > 0.2 & (original == 0 | (original > 0 & maxTlod >= original)) & mafTmp >= 0.01) {
                          if (cmedianav < 4) {
                              mafNow = mafTmp
                          } else {
                              mafNow = 0
                          }
                      } else if (maxMaf > 0.05 & (original == 0 | (original > 0 & maxTlod < original)) & mafTmp >= 0.04) {
                          if (cmedianav == 1 & mappingBias == 0) {
                              mafNow = mafTmp
                          } else {
                              mafNow = 0
                          }
                      } else {
                          mafNow = 0
                      }
                  }
                  resVector = c(resVector, c(mafNow, 0, 0, 0, refnow, altnow))
                  names(resVector)[(length(resVector)-5):length(resVector)] = c(sampleNameMaf,sampleNameMafa,sampleNameCcf,sampleNameCcfSD,sampleNameRef,sampleNameAlt)
              } #for each sample
              
              ssb = ssb/ssbc
              if (ssb >= 0.95 | ssb <= 0.05) {
                  for (j in seq(6,(3+length(cindex)), by=3)) {
                      sampleName = resColnames[j]
                      sampleNameMaf = paste(sampleName, "mafc", sep="")
                      resVector[sampleNameMaf] = 0
                  }
              }

              totalMaf = sum(as.numeric(resVector[grepl("mafc", names(resVector))]))
              totalAlt = sum(as.numeric(resVector[grepl("altc", names(resVector))]))
              totalRef = sum(as.numeric(resVector[grepl("refc", names(resVector))]))
              mergeMAFC = round(totalAlt/(totalAlt+totalRef), 4)
              mergeMAFA = mergeMAFC
              mergedCCF = mergeMAFC
              mergedCCFsd = mergeMAFC
              resVector = c(resVector, totalMaf, totalAlt, totalRef, mergeMAFC, mergeMAFA, mergedCCF, mergedCCFsd)
              names(resVector)[(length(resVector)-6):length(resVector)] = c("totalMaf", "totalAlt", "totalRef", "mergeMAFC", "mergeMAFA", "mergeCCF", "mergeCCFsd")

              resVector                #return the result

          }, cindex=cindex, original=original, resColnames=resColnames)))
    
    

    res = cbind(res, resAdd)
    res = res[which(res$totalMaf != 0),]
    return(res)
}



#adjust CCF titan for multi samples

adjust.ccf.titan.multi <- function(sampAB, samples, t, titanPath="./titan/", correctColname=FALSE) {

    purities = vector()
    
    for (i in 1:length(samples)) {
        sn = samples[i]
        message(sn)
        cnv.inputA = paste(titanPath, sn, "_nclones1.TitanCNA.segments.txt", sep="")
        #message(cnv.inputA)
        cnvA = read.delim(cnv.inputA)
        cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]
        cnvA$nt = cnvA$copynumber
        if ("minor_cn" %in% colnames(cnv.inputA)) {
            cnvA$nb = cnvA$minor_cn
        } else {
            #cnvA$nb = partialRound(cnvA$copynumber*(1-(cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion)))
            cnvA$nb = partialRound(cnvA$copynumber*(
                1-((cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion) - (1-cnvA$cellularprevalence)*0.5)
                /cnvA$cellularprevalence))
        }
        cnvA$cellularprevalence = sapply(cnvA$cellularprevalence, function(x){if (x == 1){0.99} else {x}})
        
        cnvSeqNames = cnvA$chrom
        if (!grepl("chr",cnvA$chrom[1])){
            cnvSeqNames = paste("chr",cnvA$chrom,sep="")
        }
        snvSeqNames = sampAB$chr
        if (!grepl("chr",sampAB$chr[1])){
            snvSeqNames = paste("chr",sampAB$chr,sep="")
        }
        cnvRangeA = GRanges(seqnames = cnvSeqNames, ranges = IRanges(cnvA$loc.start, end=cnvA$loc.end), strand=rep('+',dim(cnvA)[1]))
        snvRange = GRanges(seqnames = snvSeqNames, ranges = IRanges(sampAB$pos, end=sampAB$pos), strand=rep('+',dim(sampAB)[1]))
        foA = findOverlaps(snvRange, cnvRangeA)
        
        pa1 = sapply(1:dim(sampAB)[1], function(x) {
                         if (x %in% queryHits(foA)){cnvA$cellularprevalence[subjectHits(foA)[match(x, queryHits(foA))]]} else {0}}) 
        #pu1 = 1-cnvA$normalproportion[1]
        con1 = as.numeric(cnvA$normalproportion[1])
        pu1 = 1 - (con1 + (1-max(as.numeric(pa1)))*(1-con1))
        message(pu1)
        purities = append(purities, pu1)
        names(purities)[length(purities)] = sn
        nt1 = sapply(1:dim(sampAB)[1], function(x) {
                         if (x %in% queryHits(foA)){cnvA$nt[subjectHits(foA)[match(x, queryHits(foA))]]} else {2}})
        nb1 = sapply(1:dim(sampAB)[1], function(x) {
                         if (x %in% queryHits(foA)){cnvA$nb[subjectHits(foA)[match(x, queryHits(foA))]]} else {1}})
        seg1 = sapply(1:dim(sampAB)[1], function(x) {
                         if (x %in% queryHits(foA)){subjectHits(foA)[match(x, queryHits(foA))]} else {0}})   #CNA segments
        sampAB = data.frame(sampAB, pu=pu1, pa=pa1, nt=nt1, nb=nb1, seg=seg1)
        colnames(sampAB)[(dim(sampAB)[2]-4):dim(sampAB)[2]] = paste(sn,colnames(sampAB)[(dim(sampAB)[2]-4):dim(sampAB)[2]], sep="")
        if ( correctColname == TRUE ) {
            colnames(sampAB) = gsub("\\.","-",colnames(sampAB))
        }
    }

    
    for( i in 1:dim(sampAB)[1]) {  # rescale the maf and calculate CCF
        
        foundSites = 0 #count how many sites found
        #message(i)
        for (j in 1:length(samples)) {
            sn = samples[j]
            #message(sn)
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            #message(maf1)
            if (maf1 > t) {
                foundSites = foundSites+1
            }
        }
        #message(foundSites)
        for (j in 1:length(samples)) {
            sn = samples[j]

            pa1 = as.numeric(sampAB[i, match(paste(sn, "pa", sep=""), colnames(sampAB))])
            nt1 = as.numeric(sampAB[i, match(paste(sn, "nt", sep=""), colnames(sampAB))])
            nb1 = as.numeric(sampAB[i, match(paste(sn, "nb", sep=""), colnames(sampAB))])
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            refc1 = as.numeric(sampAB[i, match(paste(sn, "refc", sep=""), colnames(sampAB))])
            altc1 = as.numeric(sampAB[i, match(paste(sn, "altc", sep=""), colnames(sampAB))])
            pu1 = as.numeric(sampAB[i, match(paste(sn, "pu", sep=""), colnames(sampAB))])       #cell purity
            pu1 = nt1*pu1/(nt1*pu1+2*(1-pu1))                                                #effective purity
            
            if (maf1 > 0) {
                if (maf1 > t & foundSites >= 2) {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,nt1,nb1,"unknown")
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                } else {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,nt1,nb1,"late")
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    #sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = maf1/pu1
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                }
            }

            if ( sn == names(purities)[which.max(purities)] ) {       #rescale for merged MAF
                maf1 = as.numeric(sampAB[i, match("mergeMAFC", colnames(sampAB))])
                if (maf1 > 0) {
                    if (maf1 > t & foundSites >= 2) {
                        CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,nt1,nb1,"unknown")
                        sampAB[i, match("mergeCCF", colnames(sampAB))] = as.numeric(CCF1[3])
                        sampAB[i, match("mergeCCFsd", colnames(sampAB))] = as.numeric(CCF1[4])
                        sampAB[i, match("mergeMAFA", colnames(sampAB))] = as.numeric(CCF1[1])/2
                    } else {
                        CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,nt1,nb1,"late")
                        sampAB[i, match("mergeCCF", colnames(sampAB))] = as.numeric(CCF1[3])
                        sampAB[i, match("mergeCCFsd", colnames(sampAB))] = as.numeric(CCF1[4])
                        sampAB[i, match("mergeMAFA", colnames(sampAB))] = as.numeric(CCF1[1])/2
                    }
                }
            } #for merged MAF
        }     #for each sample
        
    }
    return(sampAB)
}



adjust.ccf.titan.single <- function(sampAB, cnv.inputA, minAF=0.05) {    #may not be working, needs revision
    cnvA = read.delim(cnv.inputA)
    cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]
    cnvA$nt = cnvA$copynumber
    if ("minor_cn" %in% colnames(cnv.inputA)) {
        cnvA$nb = cnvA$minor_cn
    } else {
        cnvA$nb = partialRound(cnvA$copynumber*(
            1-((cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion) - (1-cnvA$cellularprevalence)*0.5)
            /cnvA$cellularprevalence))
    }
    cnvA$cellularprevalence = sapply(cnvA$cellularprevalence, function(x){if (x == 1){0.99} else {x}})
    
    cnvRangeA = GRanges(seqnames = paste("chr",cnvA$chrom,sep=""), ranges = IRanges(cnvA$loc.start, end=cnvA$loc.end), strand=rep('+',dim(cnvA)[1]))
    snvRange = GRanges(seqnames = paste("chr",sampAB$chr,sep=""), ranges = IRanges(sampAB$pos, end=sampAB$pos), strand=rep('+',dim(sampAB)[1]))
    foA = findOverlaps(snvRange, cnvRangeA)
    
    pa1 = sapply(1:dim(sampAB)[1], function(x) {
                     if (x %in% queryHits(foA)){cnvA$cellularprevalence[subjectHits(foA)[match(x, queryHits(foA))]]} else {0}}) 
    pu1 = 1-cnvA$normalproportion[1]
    nt1 = sapply(1:dim(sampAB)[1], function(x) {
                     if (x %in% queryHits(foA)){cnvA$nt[subjectHits(foA)[match(x, queryHits(foA))]]} else {2}})
    nb1 = sapply(1:dim(sampAB)[1], function(x) {
                     if (x %in% queryHits(foA)){cnvA$nb[subjectHits(foA)[match(x, queryHits(foA))]]} else {1}})
    
    sampAB = data.frame(sampAB, pu1=pu1,pa1=pa1,nt1=nt1,nb1=nb1)

    for( i in 1:dim(sampAB)[1]) {  # rescale the maf
        pu1 = as.numeric(sampAB$pu1[i])
        pa1 = as.numeric(sampAB$pa1[i])
        nt1 = as.numeric(sampAB$nt1[i])
        nb1 = as.numeric(sampAB$nb1[i])
        maf1 = sampAB$maf1[i]
        refc1 = sampAB$refc1[i]
        altc1 = sampAB$altc1[i]

        CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,nt1,nb1)
        sampAB$maf1[i] = as.numeric(CCF1[1])/2
    }
    return(sampAB)
}



# plot pair-wise multi sample plot

plotRes.multi.matrix.pdf <- function(sampAB, sroot, samples, pdfsize = 16) {
    combinations = combn(length(samples),2)
    nplots = dim(combinations)[2]
    ndims = length(samples)-1

    resStats = list()
    pdf(file = paste(sroot,"multi_hist.pdf",sep=""), width=pdfsize, height=pdfsize)
    par(mfcol=c(ndims,ndims))
    pindex = 0
    for (ci in 1:ndims) {                          #for each column
        noplotRow = vector()
        if (ci > 1) {
            noplotRow = (1:ndims)[1:(ci-1)]         #which row do not plot?
        }
        for (ri in 1:ndims) {                       #for each row
            if (ri %in% noplotRow) {               #no plot
                plot(NULL,NULL,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
            } else {                               #plot
                pindex = pindex + 1
                pair = combinations[,pindex]
                sn1 = samples[pair[1]]             #samplename1
                sn1mafa = paste(sn1, "mafa",sep="")
                sn1d = paste(sn1, "d",sep="")
                sn1short = gsub("CRCTumor","",sn1)
                sn1short = gsub("Recurrence","Rec",sn1short)
                sn1short = gsub(paste(sroot,"_",sep=""),"",sn1short)
                sn1short = gsub(paste(sroot,"-",sep=""),"",sn1short)
                sn2 = samples[pair[2]]             #samplename2
                sn2mafa = paste(sn2, "mafa",sep="")
                sn2d = paste(sn2, "d",sep="")
                sn2short = gsub("CRCTumor","",sn2)
                sn2short = gsub("Recurrence","Rec",sn2short)
                sn2short = gsub(paste(sroot,"_",sep=""),"",sn2short)
                sn2short = gsub(paste(sroot,"-",sep=""),"",sn2short)
                main.title = paste(sroot,sn1short,sn2short,"a.", sep=" ")
                statName = paste(sroot, sn1short, sn2short, "stats", sep="_")
                resStats[[statName]] = plotRes.multi.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, sn1mafa, sn2mafa, 0.05, plotDepth=F, plotDensity=F, pdf=F)
                #names(resStats)[length(resStats)] = statName
            }
        }
    }
    dev.off()
    return(resStats)
}



plotRes.multi.pdf <- function(sampAB, sampName, sn1, sn2, minAF, ratio=1, plotAF=TRUE, plotDensity=TRUE, plotDepth=TRUE, pdf=TRUE) {

  sn1s = gsub("mafc","",sn1)
  sn1s = gsub("mafa","",sn1s)
  sn1s = gsub("ccf","",sn1s)
  sn2s = gsub("mafc","",sn2)
  sn2s = gsub("mafa","",sn2s)
  sn2s = gsub("ccf","",sn2s)
  nb1i = match(paste(sn1s, "nb", sep=""),colnames(sampAB))
  nb2i = match(paste(sn2s, "nb", sep=""),colnames(sampAB))
  dp1i = match(paste(sn1s, "d", sep=""),colnames(sampAB))
  dp2i = match(paste(sn2s, "d", sep=""),colnames(sampAB))
  sampAB = sampAB[which(grepl(paste(sn1s,"\\[",sep=""), sampAB$somatic) | grepl(paste(sn2s,"\\[",sep=""), sampAB$somatic)),]
  sampAB = sampAB[which((sampAB[,nb1i] > 0 & sampAB[,nb2i] > 0) | (sampAB[,nb1i] == 0 & sampAB[,nb2i] == 0)),]                 #remove different LOH locations
  
  maf1Index = match(sn1, colnames(sampAB))
  maf2Index = match(sn2, colnames(sampAB))
  mafa1Index = match(paste(sn1s,"mafa",sep=""), colnames(sampAB))    #for mafa
  mafa2Index = match(paste(sn2s,"mafa",sep=""), colnames(sampAB))    #for mafa

  #check depth power to reject a presence of a mutation
  depthPowerKeep <- as.vector(apply(sampAB, 1, function(x,mafa1i,mafa2i,dp1i,dp2i) {
                                        if(as.numeric(x[mafa1i]) == 0){vaf = as.numeric(x[mafa2i])
                                                    if (vaf > 1){FALSE}  else if (pbinom(0,as.numeric(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else if (as.numeric(x[mafa2i]) == 0){vaf = as.numeric(x[mafa1i])
                                                    if (vaf > 1){FALSE}  else if (pbinom(0,as.numeric(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else {TRUE}
                                    }, mafa1i=mafa1Index,mafa2i=mafa2Index,dp1i=dp1i,dp2i=dp2i))
  sampAB = sampAB[depthPowerKeep,]

  subMuts = subclonalMut(sampAB, sn1s, sn2s, minAF)  #subclonal mutations
  
  allA_Rows = which(sampAB[,maf1Index] > minAF & sampAB[,mafa1Index] > minAF & sampAB[,maf1Index] <= 1)
  subA_Rows = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))
  ssA_Rows  = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1 & sampAB[,maf2Index] == 0))
  sampAh = hist(sampAB[allA_Rows, maf1Index]/ratio, breaks=(max(sampAB[allA_Rows, maf1Index]/ratio)-min(sampAB[allA_Rows, maf1Index]/ratio)+0.01)/0.02,  plot=F)
  sampAhsub = hist(sampAB[subA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  sampAhss = hist(sampAB[ssA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  ylimup = max(sampAh$count)
  

  allB_Rows = which(sampAB[,maf2Index] > minAF & sampAB[,mafa2Index] > minAF & sampAB[,maf2Index] <= 1)
  subB_Rows = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))
  ssB_Rows  = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1 & sampAB[,maf1Index] == 0))
  sampBh = hist(sampAB[allB_Rows, maf2Index]/ratio, breaks=(max(sampAB[allB_Rows, maf2Index]/ratio)-min(sampAB[allA_Rows, maf2Index]/ratio)+0.01)/0.02, plot=F)
  sampBh$counts = sampBh$counts*(-1)
  sampBhsub = hist(sampAB[subB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhsub$counts = sampBhsub$count*(-1)
  sampBhss = hist(sampAB[ssB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhss$counts = sampBhss$counts*(-1)
  ylimdown = min(sampBh$count)

  #set ylim up and down equal to the same scale
  if (abs(ylimup) >= abs(ylimdown)) {
      ylimdown = (-1)*ylimup
  } else {
      ylimup = (-1)*ylimdown
  }
      
  
  ssfr = seq(0.05,0.24,by=0.01)
  ssfs = as.vector(sapply(ssfr, function(x, sAB, maf1I, maf2I){ length(which((sAB[,maf1I] >= x & sAB[, maf1I] < x+0.01 & sAB[, maf2I] == 0) | (sAB[,maf2I] >= x & sAB[, maf2I] < x+0.01 & sAB[, maf1I] == 0)))/
                                             length(which((sAB[,maf1I] >= x & sAB[, maf1I] < x+0.01) | (sAB[,maf2I] >= x & sAB[, maf2I] < x+0.01)))}, sAB = sampAB, maf1I=maf1Index, maf2I=maf2Index))
  
  if (plotAF == TRUE){

      if (pdf == TRUE) {
          pdf(file = paste(sampName, "hist.pdf", sep="_"), width = 8, height = 8)
      }
      par(mar=c(4.5,5,4.5,0))

      plot( sampAh, col=rgb(0,0,0,1/4), xlim=c(0, 1), ylim=c(ylimdown,ylimup), border=F, ylab="# of Mutations", xlab="Allele Frequency", axes = F, main = sampName,cex.lab = 2.3, cex.main = 2.3 )  # first histogram
      #plot( sampAhsub, col=rgb(242/255,108/255,108/255,1), add=T, border=F )   #subclonal red
      plot( sampAhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green set border
      #plot( sampAhss, col=rgb(248/255,181/255,53/255,1), add=T, border=F )     #site specific yellow
      plot( sampAhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue

      
      plot( sampBh,col=rgb(0,0,0,1/4), border=F, add=T )  # second histogram
      #plot( sampBhsub, col=rgb(242/255,108/255,108/255,1), add=T, border=F )   #subclonal red
      plot( sampBhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green
      #plot( sampBhss, col=rgb(43/255,98/255,223/255,1), add=T, border=F)       #site specific blue old
      plot( sampBhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
      sampAhss$counts = 0                                  #make black line
      sampBhss$counts = 0                                  #make black line
      plot( sampAhss, col="black", add=T, border=F)        #make black line
      plot( sampBhss, col="black", add=T, border=F)        #make black line
      fitxBstart = sampBh$breaks[match(max(sampBh$density[1:10]),sampBh$density)]
      fitxB = seq(fitxBstart,0.25,by=0.01)


      sn1s = gsub(".+(Core\\d+)","\\1",sn1, perl=T, ignore.case=T)
      sn1s = gsub(".+(Cor\\d+)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Sec(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Primary\\_(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Primary\\_(\\w+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Met\\_(\\d+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub(".+(Met\\_(\\w+)?)","\\1",sn1s, perl=T, ignore.case=T)
      sn1s = gsub("Primary","Pri",sn1s, perl=T, ignore.case=T)
      sn1s = gsub("mafc","",sn1s)
      sn1s = gsub("mafa","",sn1s)
      sn1s = gsub("CRCTumor","",sn1s)
      sn1s = gsub("HCT116_","",sn1s)
      sn2s = gsub(".+(Core\\d+)","\\1",sn2, perl=T, ignore.case=T)
      sn2s = gsub(".+(Cor\\d+)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Sec(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Primary\\_(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Primary\\_(\\w+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Met\\_(\\d+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub(".+(Met\\_(\\w+)?)","\\1",sn2s, perl=T, ignore.case=T)
      sn2s = gsub("Primary","Pri",sn2s, perl=T, ignore.case=T)
      sn2s = gsub("mafc","",sn2s)
      sn2s = gsub("mafa","",sn2s)
      sn2s = gsub("CRCTumor","",sn2s)
      sn2s = gsub("HCT116_","",sn2s)
      message(sn1s)
      message(sn2s)
      
      axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=1.7)
      axis(side=2,at=decideTickAt(ylimdown, ylimup),labels=allAbs(decideTickAt(ylimdown, ylimup)),cex.axis=1.7)
      text(x=0.5,y=5*ylimdown/6,labels=paste(sn2s,",",length(which(sampAB[,mafa2Index] > minAF)), "SSNVs", sep = " "), cex=2.2)
      text(x=0.5,y=ylimdown,labels=paste("Site Specific,",length(which(sampAB[,mafa2Index] > minAF & sampAB[,mafa1Index]==0)), "SSNVs", sep = " "), cex=2)
      text(x=0.5,y=ylimup,labels=paste(sn1s,",",length(which(sampAB[,mafa1Index] > minAF)), "SSNVs", sep = " "), cex=2.2)
      text(x=0.5,y=5*ylimup/6,labels=paste("Site Specific,",length(which(sampAB[,mafa1Index] > minAF & sampAB[,mafa2Index]==0)), "SSNVs", sep = " "), cex=2)
      ssfrac = length(which((sampAB[,maf2Index] > 0 & sampAB[,maf1Index] == 0 & sampAB[,maf2Index] < 0.25) | (sampAB[,maf1Index] > 0 & sampAB[,maf2Index] == 0 & sampAB[,maf1Index] < 0.25)))/length(which((sampAB[,maf2Index] > 0 & sampAB[,maf2Index] < 0.25) | (sampAB[,maf1Index] > 0 & sampAB[, maf1Index] < 0.25)))
      message(ssfrac)

      #stats starting from here!
      jsd = round(subMuts$JSD,4)
      fst = round(subMuts$FST,4)
      text(x=0.6,y=2*ylimup/3, labels=bquote(paste("JSD"["sub"], " = ", .(jsd))),cex=2.2)
      text(x=0.6,y=1*ylimup/2, labels=bquote(paste("FST"["sub"], " = ", .(fst))),cex=2.2)

      npub = length(which(sampAB[,mafa1Index] > minAF)) - length(subMuts$subAi)
      legend(0.65,ylimdown/4, legend=c(paste("Public ","(",npub,")",sep=""),"Shared Sub","Site Specific"),
             col=c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1)),pch=15,bty="n",cex=1.7)
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotDensity == TRUE){
      #two way density plot
      pdf(file = paste(sampName, "density.pdf", sep="_"), width = 8, height = 8)
      par(mar=c(4.5,5,4.5,2))
      smkey(sampAB[,maf1Index],sampAB[,maf2Index],xlab=paste("MAF",sn1s,sep=" "), ylab=paste("MAF",sn2s,sep=" "),
            main = sampName, xlim=c(0,1), ylim=c(0,1),cex.lab = 2.3, cex.main = 2.3,cex.axis=1.7)
      if (length(which(sampAB$dron == 1)) > 0) {
          dronindex = which(sampAB$dron == 1)
          pointLabel(sampAB[dronindex,maf1Index],sampAB[dronindex,maf2Index],labels=as.character(sampAB$geneName[dronindex]),
                     col=decideDriverColor(sampAB[dronindex,maf1Index],sampAB[dronindex,maf2Index]),cex=1.7,font=2)
      }
      segments(0,0,0.8,0.8)
      dev.off()
  }


   if (plotDepth == TRUE) {
      pdf(file = paste(sampName, "depth.pdf", sep="_"), width = 10, height = 6.6)
      par(mfrow=c(2,3))
      ymaxdepth = max(quantile(sampAB[which(sampAB[,mafa1Index] > minAF),dp1i], prob=seq(0,1,0.05))["95%"],
          quantile(sampAB[which(sampAB[,mafa1Index] > minAF),dp2i], prob=seq(0,1,0.05))["95%"])
      calledA = which(sampAB[,mafa1Index] > minAF & sampAB[,mafa1Index] <= 1)
      calledB = which(sampAB[,mafa2Index] > minAF & sampAB[,mafa2Index] <= 1)
      plot(sampAB[calledA,mafa1Index],sampAB[calledA,dp1i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("VAF",sn1s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste("mutations called in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      legend("topright",legend=c("Site Specific","Shared Sub","Public"),col=c(rgb(31/255,120/255,180/255,1),rgb(178/255,223/255,138/255,1),rgb(0,0,0,1/4)),pch=19, cex=1.3)
      plot(sampAB[calledA,mafa1Index],sampAB[calledA,dp2i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("VAF",sn1s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste("mutations called in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledA,dp1i],sampAB[calledA,dp2i],xlim=c(0,ymaxdepth),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("depth",sn1s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste("mutations called in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,mafa2Index],sampAB[calledB,dp2i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("VAF",sn2s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste("mutations called in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,mafa2Index],sampAB[calledB,dp1i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("VAF",sn2s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste("mutations called in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,dp2i],sampAB[calledB,dp1i],xlim=c(0,ymaxdepth),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("depth",sn2s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste("mutations called in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      dev.off()
  }

  return(subMuts)
  #return(ssfs)
}


subSSColor <- function(allindex, subindex, ssindex){
    subsscolor = as.vector(sapply(allindex, function(x, subindex, ssindex){
                            if (x %in% ssindex) {
                                rgb(31/255,120/255,180/255,1)
                            } else if (x %in% subindex) {
                                rgb(178/255,223/255,138/255,1)
                            } else {
                                rgb(0,0,0,1/4)
                            }},subindex=subindex, ssindex=ssindex))
    return(subsscolor)
}


decideDriverColor <- function(x,y){
    dc = rep(rgb(0,0,0,3/4), length(x))
    for (i in 1:length(x)){
        if (x[i] > 0.75 & y[i] > 0.75) {
            dc[i] = rgb(1,1,1,3/4)
        }
    }
    return(dc)
}

nearestDecimal <- function(x) {
    r = x %% 10
    nearD = 0
    if (r > 5) {
        nearD = x + (10-r)
    } else {
        nearD = x - r
    }
    return(nearD)
}

decideTickAt <- function(ylimdown, ylimup) {  #assuming they are equal in abs

    if (abs(ylimdown) != abs(ylimup)){
        stop("ylimdown and up should be equal in abs!")
    }

    realmax = nearestDecimal(ylimup)
    incre = nearestDecimal(realmax/5)
    if (incre == 0) {
        incre = 5
    }
    njump = floor(realmax/incre)
    tailgap = realmax %% incre

    ylimup = realmax
    ylimdown = (-1)*realmax
    tickAt = c(ylimdown, seq(ylimdown + tailgap, 0, by=incre), seq(0, ylimup - tailgap, by=incre), ylimup)
    return(tickAt)

}


allAbs <- function(x) {
    for (i in 1:length(x)){
        x[i] = abs(x[i])
    }
    return(x)
}


computeCCF <- function(f, A, S, pu, pa, nt, nb, prior="unknown") {
    ccf = 0
    ccf2 = 0
    sd = 0
    cc <- seq(0.02, 1, by = 0.01)
    evoType = "A1/A2/B/C"
    sigTh = 0.90
    N = A + S
    nc = nt * pa + 2 * (1 - pa)
    #message(paste(f,A,S,pu,pa,nt,nb,nc,sep="  "))
    if (nb == 1 & nt == 2) {   #normal diploid
        ccf = 2*(f/pu)
        ff = pu*cc/2
        Ms = computeSD(N, S, ff)
        ccf2 <- Ms$M1
        sd <- Ms$SD
    }
    else if (nt == 1) {        #het deletion
        ccf = (f/pu)*nc
        ff.C <- pu*cc/nc                       #dbinom
        Ms.C <- computeSD(N, S, ff.C)          #dbinom
        ccf2 <- Ms.C$M1                        #dbinom
        sd <- Ms.C$SD
        fh.ea <- (pa * (nt - nb) + 1 - pa)/nc
        fl.ea <- (pa * (nt - nb))/nc
        fh.t <- pa/nc
        fh.e <- (1 - pa)/nc
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cp.A <- pEarly.a/Ptot
        cp.CD <- 1 - cp.A
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        cp.AC <- 1 - cp.D
        cp.AD <- 1 - cp.C
        if (cp.A >= sigTh) {
            evoType <- "A1"
        } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh) {
            evoType <- "B/C"
        } else if (cp.C >= sigTh) {
            evoType <- "B"
        } else if (cp.D >= sigTh) {
            evoType <- "C"
        } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh) {
            evoType <- "A1/B"
        } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh) {
            evoType <- "A1/C"
        }
    } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
        fh.ea <- (pa * (nt - nb) + 1 - pa)/nc
        fl.ea <- (pa * (nt - nb))/nc
        fh.t <- pa/nc
        fh.e <- (1 - pa)/nc
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cpEarly.a <- pEarly.a/Ptot
        cpLate.eup <- 1 - cpEarly.a
        cpLate <- pLate/Ptot
        cpEup <- pEuploid/Ptot
        #message(paste(pEarly.a, pLate, pEuploid, Ptot, sep="  "))
        if (Ptot > 0) {
            if (cpEarly.a >= sigTh){
                evoType <- "A1"
            } else if (cpLate.eup >= sigTh){
                evoType <- "B/C"
            } else if (cpLate >= sigTh){
                evoType <- "B"
            } else if (cpEup >= sigTh){
                evoType <- "C"
            }
        }
        allprobs = c(pEarly.a, pLate, pEuploid)
        names(allprobs) = c("pEarly.a", "pLate", "pEuploid")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        if (maxType == "pEarly.a" & prior != "late") {
            ccf = (f/pu)*nc - (nt - nb - 1)*pa
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            sd <- Ms.A$SD
        } else {
            ccf = (f/pu)*nc
            ff.C <- pu*cc/nc                        #dbinom
            Ms.C <- computeSD(N, S, ff.C)           #dbinom
            ccf2 <- Ms.C$M1                         #dbinom
            sd <- Ms.C$SD
        }
    } else if (nb >= 1 & nt > 2) {
        fh.ea <- (pa * (nt - nb) + 1 - pa)/nc
        fl.ea <- (pa * (nt - nb))/nc
        fh.eb <- (nb * pa + 1 - pa)/nc
        fl.eb <- nb * pa/nc
        fh.t <- pa/nc
        fh.e <- (1 - pa)/nc
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
        cp.A <- pEarly.a/Ptot
        cp.B <- pEarly.b/Ptot
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        cp.AB <- 1 - cp.C - cp.D
        cp.AC <- 1 - cp.B - cp.D
        cp.AD <- 1 - cp.B - cp.D
        cp.BC <- 1 - cp.A - cp.D
        cp.BD <- 1 - cp.A - cp.C
        cp.CD <- 1 - cp.A - cp.B
        cp.ABC <- 1 - cp.D
        cp.ABD <- 1 - cp.C
        cp.ACD <- 1 - cp.B
        cp.BCD <- 1 - cp.A
        if (cp.A >= sigTh){
            evoType = "A1"
        } else if (cp.B >= sigTh){
            evoType <- "A2"
        } else if (cp.C >= sigTh){
            evoType <- "B"
        } else if (cp.D >= sigTh){
            evoType <- "C"
        } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh){
            evoType <- "B/C"
        } else if (cp.AB >= sigTh & cp.A < sigTh & cp.B < sigTh){
            evoType <- "A1/A2"
        } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh){
            evoType <- "A1/B"
        } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh){
            evoType <- "A1/C"
        } else if (cp.BC >= sigTh & cp.B < sigTh & cp.C < sigTh){
            evoType <- "A2/B"
        } else if (cp.BD >= sigTh & cp.B < sigTh & cp.D < sigTh){
            evoType <- "A2/C"
        } else if (cp.BCD >= sigTh & cp.BC < sigTh & cp.BD < sigTh & cp.CD < sigTh & cp.B < sigTh & cp.C < sigTh & cp.D < sigTh){
            evoType <- "A2/B/C"
        } else if (cp.ABC >= sigTh & cp.BC < sigTh & cp.AB < sigTh & cp.AC < sigTh & cp.B < sigTh & cp.C < sigTh & cp.A < sigTh){
            evoType <- "A1/A2/B"
        } else if (cp.ABD >= sigTh & cp.AB < sigTh & cp.AD < sigTh & cp.BD < sigTh & cp.B < sigTh & cp.D < sigTh & cp.A < sigTh){
            evoType <- "A1/A2/C"
        } else if (cp.ACD >= sigTh & cp.AC < sigTh & cp.AD < sigTh & cp.CD < sigTh & cp.A < sigTh & cp.D < sigTh & cp.C < sigTh){
            evoType <- "A1/B/C"
        }
        
        allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid)
        names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        if (maxType == "pEarly.a" & prior != "late") {
            ccf = (f/pu)*nc - (nt - nb - 1)*pa  #early A1
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            sd <- Ms.A$SD
        } else if (maxType == "pEarly.b" & prior != "late") {
            ccf = (f/pu)*nc - (nb - 1)* pa      #early A2
            ff.B <- pu*(cc - pa + nb * pa)/nc      #dbinom
            Ms.B <- computeSD(N, S, ff.B)          #dbinom
            ccf2 <- Ms.B$M1                         #dbinom
            sd <- Ms.B$SD
        } else {
            ccf = (f/pu)*nc                     #other
            ff.C <- pu*cc/nc                       #dbinom
            Ms.C <- computeSD(N, S, ff.C)          #dbinom
            ccf2 <- Ms.C$M1                         #dbinom
            sd <- Ms.C$SD
        }
    }
    if ( f > 0.1 & ccf >= 1.6 ) {    #correct for over-adjustment
        if (evoType != "A1") {
            if ( (nt-nb) >= 3 ){
                ccf = (f/pu)*2
            } else if ( nt >= 2 & (nt-nb) < 3 ) {
                ccf = (f/pu)*nc - (nt - nb - 1)*pa
            } else {
                ccf = (f/pu)*nc
            }
        } else {
            ccf = f*2
        }
    }
    if ( f > 0.45 & ccf > 2 ) {
        ccf = f*2
    }
    return(c(ccf, evoType, ccf2, sd))
}


computeSD <- function(N, S, f) {
    M1list <- c()
    M2list <- c()
    MLElist <- c()
    cc <- seq(0.02, 1, by = 0.01)
    for (ii in 1:length(N)) {
        PF <- sum(dbinom(S[ii], N[ii], f), na.rm = TRUE)
        M1 <- sum(dbinom(S[ii], N[ii], f) * cc, na.rm = TRUE)/PF
        M2 <- sum(dbinom(S[ii], N[ii], f) * cc^2, na.rm = TRUE)/PF
        M1list <- c(M1list, M1)
        M2list <- c(M2list, M2)
        MLElist <- c(MLElist, cc[which.max(dbinom(S[ii], N[ii], f))])
    }
    return(list(M1 = MLElist, SD = sqrt(M2list - M1list^2)))
}


partialRound <- function(x) {
    r = x
    for (i in 1:length(x)) {
        r[i] = round(x[i])
        if (abs(r[i]-x[i]) > 0.3 & x[i] > 1) {
            r[i] = round(x[i],1)
        }
        if (r[i] < 0) {
            r[i] = 0
        }
    }
    return(r)
}


effectivePurity <- function(p, Nt) {
    pu = Nt*p/(Nt*p + 2*(1-p))
    return(pu)
}


pmf <- function(muts, minAF=0.04) {                          # calculate the probability distribution from AF data
    bin = seq(minAF,0.8,0.02)                     # the bins I used
    mhist = hist(muts, breaks = bin, plot=F)
    mfreq = mhist$counts/sum(mhist$counts)

    for(k in 1:length(mfreq)) {
        if(mfreq[k]==0) {
            mfreq[k] = 0.000001
        }
    }
    return(mfreq)
}


JS.divergence <- function(sub1,sub2, minAF=0.04) {
    sub1 = sub1[which(sub1 <= 0.8 & sub1 >= minAF)]
    sub2 = sub2[which(sub2 <= 0.8 & sub2 >= minAF)]
    freqs1 = pmf(sub1, minAF)
    #message(paste(freqs1,"",sep=" "))
    freqs2 = pmf(sub2, minAF)
    #message(paste(freqs2,"",sep=" "))
    m<-0.5*(freqs1 + freqs2)
    JS<-0.5*(sum(freqs1*log(freqs1/m)) + sum(freqs2*log(freqs2/m)))
    return(JS)
}


subclonalMut <- function(sampAB, snA, snB, minAF=0.08, statsAF=0.08, highAF=0.2, ratio=1)   {                  #determinine subclonal mutations
    ccfAi = match(paste(snA, "ccf", sep=""), colnames(sampAB))
    ccfBi = match(paste(snB, "ccf", sep=""), colnames(sampAB))
    ccfsdAi = match(paste(snA, "ccfSD", sep=""), colnames(sampAB))
    ccfsdBi = match(paste(snB, "ccfSD", sep=""), colnames(sampAB))
    mafaAi = match(paste(snA, "mafa", sep=""), colnames(sampAB))
    mafaBi = match(paste(snB, "mafa", sep=""), colnames(sampAB))
    nbAi = match(paste(snA, "nb", sep=""), colnames(sampAB))
    nbBi = match(paste(snB, "nb", sep=""), colnames(sampAB))
    depthAi = match(paste(snA, "d", sep=""), colnames(sampAB))
    depthBi = match(paste(snB, "d", sep=""), colnames(sampAB))

    gNIndex = match("geneName",colnames(sampAB))
    gLIndex = match("geneLoc",colnames(sampAB))
    fCIndex = match("functionalClass",colnames(sampAB))
    
    # for JSD
    subAi = which( sampAB[,mafaAi] > minAF &
            (((sampAB[,ccfAi]+2.58*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+2.58*sampAB[,ccfsdBi]) < 1) &       #either one side is below CCF+sd 1
            (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)) &                                               #or one side is below VAF 0.25
                ((sampAB[,mafaBi] == 0 & (sampAB[,nbBi] != 0 | sampAB[,nbAi] == 0)) | sampAB[,mafaBi] != 0) )  #and the other side VAF > 0 or VAF == 0 (either not LOH or the other side is the same LOH)
    mutsA = sampAB[subAi,mafaAi]/ratio

    subBi = which( sampAB[,mafaBi] > minAF &
            (((sampAB[,ccfAi]+2.58*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+2.58*sampAB[,ccfsdBi]) < 1) &       #either one side is below CCF+sd 1
            (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)) &                                               #or one side is below VAF 0.25
            ((sampAB[,mafaAi] == 0 & (sampAB[,nbAi] != 0 | sampAB[,nbBi] == 0)) | sampAB[,mafaAi] != 0))       #and the other side VAF > 0 or VAF == 0 (either not LOH or the other side is the same LOH)
    mutsB = sampAB[subBi,mafaBi]/ratio
    JSD = JS.divergence(mutsA, mutsB, minAF=statsAF)

    # for mutational function dNdS
    subTi = union(subAi, subBi)
    pubTi = which( ((sampAB[,ccfAi]+2.58*sampAB[,ccfsdAi]) >= 1 & (sampAB[,ccfBi]+2.58*sampAB[,ccfsdBi]) >= 1) |
                      (sampAB[,mafaAi] >= 0.25 & sampAB[,mafaBi] >= 0.25) )
    subMutGeneFunc = sampAB[subTi, c(gNIndex, gLIndex, fCIndex)]
    pubMutGeneFunc = sampAB[pubTi, c(gNIndex, gLIndex, fCIndex)]

    # for FST
    mutsSub = sampAB[which((sampAB[,mafaAi] >= statsAF | sampAB[,mafaBi] >= statsAF) & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 &
                     ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 & (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                     (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),]
    mutsSub = data.frame(maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi])
    FST = mean(fst.wc84(mutsSub, minAF=statsAF))
    
    # for other stats
    mutsA2 = sampAB[which( sampAB[,mafaAi] >= statsAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaAi]
    mutsAh2 = sampAB[which( sampAB[,mafaAi] >= highAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaAi]
    mutsASp2 = sampAB[which( sampAB[,mafaAi] >= statsAF & sampAB[,mafaBi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15),mafaAi]
    mutsASph2 = sampAB[which( sampAB[,mafaAi] >= highAF & sampAB[,mafaBi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15),mafaAi]
    mutsB2 = sampAB[which( sampAB[,mafaBi] >= statsAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaBi]
    mutsBh2 = sampAB[which( sampAB[,mafaBi] >= highAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaBi]
    mutsBSp2 = sampAB[which( sampAB[,mafaBi] >= statsAF & sampAB[,mafaAi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 ),mafaBi]
    mutsBSph2 = sampAB[which( sampAB[,mafaBi] >= highAF & sampAB[,mafaAi] == 0 &
                                 sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 ),mafaBi]
        
    # list for output    
    muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, fstInput=mutsSub,
        lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
        lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
        lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
        lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
        FST=FST, JSD=JSD, subMutGeneFunc = subMutGeneFunc, pubMutGeneFunc = pubMutGeneFunc)
    return(muts)
}


pubOrSub <- function(sampAB, samples, minAF=0.08, statsAF=0.08, highAF=0.2, minDep=9, ratio=1) {

    originalColNames = colnames(sampAB)
    
    CCFbelowOne = vector()          #sub
    VAFbelowQua = vector()          #sub
    CCFaboveOne = vector()          #pub
    VAFaboveQua = vector()          #pub
    aboveContri = vector()
    
    for (i in 1:length(samples)) {                                     #get subclonal ones 1st round: Raw
        sn = samples[i]
        message(sn)
        ccfi = match(paste(sn, "ccf", sep=""), colnames(sampAB))
        ccfsdi = match(paste(sn, "ccfSD", sep=""), colnames(sampAB))
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        depthi = match(paste(sn, "d", sep=""), colnames(sampAB))

        depthQualTmp = sampAB[,depthi] >= minDep

        CCFbelowOneTmp = (sampAB[,ccfi]+2.58*sampAB[,ccfsdi]) < 1      #at least one site is below CCF+sd < 1
        if ( length(CCFbelowOne) > 0 ) {
            CCFbelowOne = CCFbelowOne | (CCFbelowOneTmp & depthQualTmp)
        } else {
            CCFbelowOne = CCFbelowOneTmp & depthQualTmp
        }

        CCFaboveOneND = (sampAB[,ccfi]+2.58*sampAB[,ccfsdi]) >= 1
        CCFaboveOneTmp = CCFaboveOneND | sampAB[,depthi] <= 6
        if ( length(CCFaboveOne) > 0 ) {
            CCFaboveOne = CCFaboveOne & CCFaboveOneTmp
            aboveContri = aboveContri + as.numeric(CCFaboveOneND)
        } else {
            CCFaboveOne = CCFaboveOneTmp
            aboveContri = as.numeric(CCFaboveOneND)
        }
        
        VAFbelowQuaTmp = sampAB[,mafai] < 0.25                         #and one site is below VAF 0.25
        if ( length(VAFbelowQua) > 0 ) {
            VAFbelowQua = VAFbelowQua | (VAFbelowQuaTmp & depthQualTmp)
        } else {
            VAFbelowQua = VAFbelowQuaTmp & depthQualTmp
        }

        VAFaboveQuaTmp = sampAB[,mafai] >= 0.25 | sampAB[,depthi] <= 6
        if ( length(VAFaboveQua) > 0 ) {
            VAFaboveQua = VAFaboveQua & VAFaboveQuaTmp
        } else {
            VAFaboveQua = VAFaboveQuaTmp
        }
    }
    #message(CCFbelowOne[222])
    #message(VAFbelowQua[222])
    #message(CCFaboveOne[222])
    #message(VAFaboveQua[222])
    submutIndex = list()
    submutMafa = list()
    allsubclone = CCFbelowOne & VAFbelowQua
    allpubclone = (CCFaboveOne | VAFaboveQua) & (aboveContri >= 2)
    for (i in 1:length(samples)) {                                     #get subclonal ones 2nd round: Refine
        sn = samples[i]
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        depthi = match(paste(sn, "d", sep=""), colnames(sampAB))

        mafaQualTmp = sampAB[,mafai] >= minAF
        noDifferLOH = vector()
        
        for (j in 1:length(samples)) {
            if (i != j) {
                mafaJi = match(paste(samples[j], "mafa", sep=""), colnames(sampAB))
                nbJi = match(paste(samples[j], "nb", sep=""), colnames(sampAB))
                depthJi = match(paste(samples[j], "d", sep=""), colnames(sampAB))
                noDifferLOHTmp = (sampAB[,mafaJi] == 0 & (sampAB[,nbJi] != 0 | sampAB[,nbi] == 0)) | sampAB[,mafaJi] != 0    #the other site VAF > 0 or VAF == 0 (either not LOH or the other side is the same LOH)
                if ( length(noDifferLOH) > 0 ){
                    noDifferLOH = noDifferLOH & noDifferLOHTmp
                } else {
                    noDifferLOH = noDifferLOHTmp
                }
            }
        }   #another sample
        subi = which(allsubclone & mafaQualTmp & noDifferLOH)
        mafa = sampAB[subi,mafai]/ratio
        submutIndex = append(submutIndex, list(subi))
        names(submutIndex)[length(submutIndex)] = paste(sn, "Subi", sep="")
        submutMafa = append(submutMafa, list(mafa))
        names(submutMafa)[length(submutMafa)] = paste(sn, "mafa", sep="")
    }

    pstype = sapply(1:dim(sampAB)[1], function(x, pub, sub) {
                        pstype = "unknown"
                        if (x %in% pub){
                            pstype = "public"
                        } else {
                            for(j in 1:length(submutIndex)){
                                sampname= gsub("Subi","",names(submutIndex)[j])
                                if (x %in% submutIndex[[j]]) {
                                    if (pstype == "unknown"){
                                        pstype = paste("private=",sampname,sep="")
                                    } else {
                                        pstype = paste(pstype, sampname, sep=",")
                                    }
                                }
                            }
                        }
                        pstype
                    }, pub=which(allpubclone), sub=submutIndex)
    
    sampAB = data.frame(sampAB, pubOrSub=pstype)

    colnames(sampAB) = c(originalColNames,"pubOrSub")
    
    return(sampAB)

}



subclonalMutAll <- function(sampAB, samples, minAF=0.08, statsAF=0.08, highAF=0.2, minDep=10, ratio=1)   {                  #determinine subclonal mutations

    combinations = combn(length(samples),2)
    
    CCFbelowOne = vector()          #sub
    VAFbelowQua = vector()          #sub
    CCFaboveOne = vector()          #pub
    VAFaboveQua = vector()          #pub
    
    for (i in 1:length(samples)) {                                     #get subclonal ones 1st round: Raw
        sn = samples[i]
        ccfi = match(paste(sn, "ccf", sep=""), colnames(sampAB))
        ccfsdi = match(paste(sn, "ccfSD", sep=""), colnames(sampAB))
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        depthi = match(paste(sn, "d", sep=""), colnames(sampAB))

        depthQualTmp = sampAB[,depthi] >= minDep
        CCFbelowOneTmp = (sampAB[,ccfi]+2.58*sampAB[,ccfsdi]) < 1      #at least one site is below CCF+sd 1
        CCFbelowOne = CCFbelowOne | (CCFbelowOneTmp & depthQualTmp)
        CCFaboveOneTmp = ((sampAB[,ccfi]+2.58*sampAB[,ccfsdi]) >= 1) | sampAB[,depthi] <= 3
        CCFaboveOne = CCFaboveOne & CCFaboveOneTmp
        
        VAFbelowQuaTmp = sampAB[,mafai] < 0.25                         #and one site is below VAF 0.25
        VAFbelowQua = VAFbelowQua | (VAFbelowQuaTmp & depthQualTmp)
        VAFaboveQuaTmp = sampAB[,mafai] >= 0.25 | sampAB[,depthi] <= 3
        VAFaboveQua = VAFaboveQua & VAFaboveQuaTmp
    }


    submutIndex = list()
    submutMafa = list()
    allsubclone = CCFbelowOne & VAFbelowQua
    allpubclone = CCFaboveOne | VAFaboveQua
    for (i in 1:length(samples)) {                                     #get subclonal ones 2nd round: Refine
        sn = samples[i]
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        depthi = match(paste(sn, "d", sep=""), colnames(sampAB))

        mafaQualTmp = sampAB[,mafai] >= minAF
        noDifferLOH = vector()
        
        for (j in 1:length(samples)) {
            if (i != j) {
                mafaJi = match(paste(samples[j], "mafa", sep=""), colnames(sampAB))
                nbJi = match(paste(samples[j], "nb", sep=""), colnames(sampAB))
                depthJi = match(paste(samples[j], "d", sep=""), colnames(sampAB))
                noDifferLOHTmp = (sampAB[,mafaJi] == 0 & (sampAB[,nbJi] != 0 | sampAB[,nbi] == 0)) | sampAB[,mafaJi] != 0    #the other site VAF > 0 or VAF == 0 (either not LOH or the other side is the same LOH)
                noDifferLOH = noDifferLOH & noDifferLOHTmp
            }
        }   #another sample
        subi = which(allsubclone & mafaQualTmp & noDifferLOH)
        mafa = sampAB[subi,mafai]/ratio
        submutIndex = append(submutIndex, list(subi))
        names(submutIndex)[length(submutIndex)] = paste(sn, "Subi", sep="")
        submutMafa = append(submutMafa, list(mafa))
        names(submutMafa)[length(submutMafa)] = paste(sn, "mafa", sep="")
    }
    
    # for JSD
    JSD = vector()
    for (i in 1:dim(combinations)[2]) {
      JSD = c(JSD, JS.divergence(submutMafa[[combinations[1,i]]], submutMafa[[combinations[2,i]]], minAF=statsAF))
      names(JSD)[length(JSD)] = paste("JSD",samples[combinations[1,i]],samples[combinations[2,i]], sep="_")
    }

    
    # for mutational function dNdS
    subTi = vector()
    for (i in 1:length(submutIndex)){
        subTi = union(subTi, submutIndex[[i]])
    }
    pubTi = allpubclone
        
    gNIndex = match("geneName",colnames(sampAB))
    gLIndex = match("geneLoc",colnames(sampAB))
    fCIndex = match("functionalClass",colnames(sampAB))
    subMutGeneFunc = sampAB[subTi, c(gNIndex, gLIndex, fCIndex)]
    pubMutGeneFunc = sampAB[pubTi, c(gNIndex, gLIndex, fCIndex)]

    
    # for FST
    FST = vector()
    for (i in 1:dim(combinations)[2]) {
      mutsSubPair = which()
      FST = c(FST, mean(fst.wc84(mutsSubPair, minAF=statsAF)))
      names(JSD)[length(JSD)] = paste("JSD",samples[combinations[1,i]],samples[combinations[2,i]], sep="_")
    }

    mutsSub = sampAB[which((sampAB[,mafaAi] >= statsAF | sampAB[,mafaBi] >= statsAF) & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 &
                     ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 & (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                     (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),]
    mutsSub = data.frame(maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi])
    FST = mean(fst.wc84(mutsSub, minAF=statsAF))


    # for other stats
    mutsA2 = sampAB[which( sampAB[,mafaAi] >= statsAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaAi]
    mutsAh2 = sampAB[which( sampAB[,mafaAi] >= highAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaAi]
    mutsASp2 = sampAB[which( sampAB[,mafaAi] >= statsAF & sampAB[,mafaBi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15),mafaAi]
    mutsASph2 = sampAB[which( sampAB[,mafaAi] >= highAF & sampAB[,mafaBi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15),mafaAi]
    mutsB2 = sampAB[which( sampAB[,mafaBi] >= statsAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaBi]
    mutsBh2 = sampAB[which( sampAB[,mafaBi] >= highAF & sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 & 
            ((sampAB[,ccfAi]+3.09*sampAB[,ccfsdAi]) < 1 | (sampAB[,ccfBi]+3.09*sampAB[,ccfsdBi]) < 1) &
                (sampAB[,mafaAi] < 0.25 | sampAB[,mafaBi] < 0.25)),mafaBi]
    mutsBSp2 = sampAB[which( sampAB[,mafaBi] >= statsAF & sampAB[,mafaAi] == 0 &
                               sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 ),mafaBi]
    mutsBSph2 = sampAB[which( sampAB[,mafaBi] >= highAF & sampAB[,mafaAi] == 0 &
                                 sampAB[,depthAi] >= 15 & sampAB[,depthBi] >= 15 ),mafaBi]
        
    # list for output
    muts = list(submutIndex=submutIndex, submutMafa=submutMafa, subMutGeneFunc=subMutGeneFunc, pubMutGeneFunc=pubMutGeneFunc)
    muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, fstInput=mutsSub,
        lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
        lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
        lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
        lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
        FST=FST, JSD=JSD, subMutGeneFunc = subMutGeneFunc, pubMutGeneFunc = pubMutGeneFunc)
    return(muts)
}


fst.wc84 <- function(af, minAF=0.1) {
       mafis = which(grepl("maf", colnames(af)))
       keep = as.vector(apply(af, 1, function(x, mafis) {
                                  maxmaf = max(as.numeric(x[mafis]))
                                  if (maxmaf >= minAF) {
                                      TRUE
                                  } else {
                                      FALSE
                                  }
                              }, mafis=mafis))     #filter data
       af = af[keep,]
       fst_array = c()
       for(k in 1:nrow(af)) {
            n1 = af$depth1[k]
            n2 = af$depth2[k]
            p1 = af$maf1[k]
            p2 = af$maf2[k]
            A = 2*n1*n2/(n1+n2)
            B = 1/(n1+n2-2)*(n1*p1*(1-p1)+n2*p2*(1-p2))
            C = n1*n2/(n1+n2)*(p1-p2)^2
            fstk =1 - A*B/(C+(A-1)*B)
            fst_array = c(fst_array,fstk)
       }
       return(fst_array)
}

fst.nei86 <- function(af) {
       fst_array = c()
       for(k in 1:nrow(af)){
            p1 = af$maf1[k]
            p2 = af$maf2[k]
            p_aver = (p1+p2)/2
            fstk = (p1-p2)^2/(2*p_aver*(1-p_aver))
            fst_array = c(fst_array,fstk)
       }
       return(fst_array)
}
fst.hudson <- function(af){
    fst_array = c()
     for(k in 1:nrow(af)){
            n1 = af$depth1[k]
            n2 = af$depth2[k]
            p1 = af$maf1[k]
            p2 = af$maf2[k]
            A = (p1-p2)^2-(p1*(1-p1))/(n1-1)-(p2*(1-p2))/(n2-1)
            B = p1*(1-p2)+p2*(1-p1)
            fstk = A/B
            fst_array = c(fst_array,fstk)
    }
    return(fst_array)
}


dNdS <- function(sampAB, g.dnds) {
    g.table1 = table(sampAB$geneName[which(sampAB$pubOrSub == "public" &
                                     (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])
    g.table1 = g.table1[which(g.table1 > 0)]
    g.match1 = names(g.table1) %in% g.dnds$gene
    g.table1 = g.table1[g.match1]
    subset.dnds1 = g.dnds$NS[match(names(g.table1),g.dnds$gene)]
    g.total1 = sum(g.table1)
    g.frac1 = g.table1/g.total1
    g.norm1 = sum(g.frac1*subset.dnds1)
    f.table1 = table(sampAB$functionalClass[which(sampAB$pubOrSub == "public" &
                                                 (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])
    nonsyn1 = as.numeric(f.table1["nonsynonymous SNV"] + f.table1["stopgain"] + f.table1["stoploss"])
    syn1 = as.numeric(f.table1["synonymous SNV"])
    res.dnds1 = nonsyn1/syn1
    res.dnds1 = res.dnds1/g.norm1

    cadd.g.table1 = table(sampAB$geneName[which(sampAB$pubOrSub == "public" & sampAB$CADD_phred != '.')])
    cadd.g.table1 = cadd.g.table1[which(cadd.g.table1 > 0)]
    cadd.g.table1 = cadd.g.table1[names(cadd.g.table1) %in% g.dnds$gene]
    subset.cadd1 = g.dnds$MFLF[match(names(cadd.g.table1),g.dnds$gene)]
    cadd.g.total1 = sum(cadd.g.table1)
    cadd.g.frac1 = cadd.g.table1/cadd.g.total1
    cadd.g.norm1 = sum(cadd.g.frac1*subset.cadd1)
    cadd.f.more1 = length(which(sampAB$pubOrSub == "public" & sampAB$CADD_phred != '.' & sampAB$CADD_phred >= 20))
    cadd.f.less1 = length(which(sampAB$pubOrSub == "public" & sampAB$CADD_phred != '.' & sampAB$CADD_phred < 20))
    res.cadd1 = (cadd.f.more1/cadd.f.less1)/cadd.g.norm1
    
    g.table2 = table(sampAB$geneName[which(sampAB$pubOrSub != "public" & sampAB$pubOrSub != "unknown" &
                                         (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])
    g.table2 = g.table2[which(g.table2 > 0)]
    g.match2 = names(g.table2) %in% g.dnds$gene
    g.table2 = g.table2[g.match2]
    subset.dnds2 = g.dnds$NS[match(names(g.table2),g.dnds$gene)]
    g.total2 = sum(g.table2)
    g.frac2 = g.table2/g.total2
    g.norm2 = sum(g.frac2*subset.dnds2)
    f.table2 = table(sampAB$functionalClass[which(sampAB$pubOrSub != "public" & sampAB$pubOrSub != "unknown" &
                                                      (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])

    nonsyn2 = as.numeric(f.table2["nonsynonymous SNV"] + f.table2["stopgain"] + f.table2["stoploss"])
    syn2 = as.numeric(f.table2["synonymous SNV"])
    res.dnds2 = nonsyn2/syn2
    res.dnds2 = res.dnds2/g.norm2

    cadd.g.table2 = table(sampAB$geneName[which(grepl("private",sampAB$pubOrSub) & sampAB$CADD_phred != '.')])
    cadd.g.table2 = cadd.g.table2[which(cadd.g.table2 > 0)]
    cadd.g.table2 = cadd.g.table2[names(cadd.g.table2) %in% g.dnds$gene]
    subset.cadd2 = g.dnds$MFLF[match(names(cadd.g.table2),g.dnds$gene)]
    cadd.g.total2 = sum(cadd.g.table2)
    cadd.g.frac2 = cadd.g.table2/cadd.g.total2
    cadd.g.norm2 = sum(cadd.g.frac2*subset.cadd2)
    cadd.f.more2 = length(which(grepl("private",sampAB$pubOrSub) & sampAB$CADD_phred != '.' & sampAB$CADD_phred >= 20))
    cadd.f.less2 = length(which(grepl("private",sampAB$pubOrSub) & sampAB$CADD_phred != '.' & sampAB$CADD_phred < 20))
    res.cadd2 = (cadd.f.more2/cadd.f.less2)/cadd.g.norm2
    
    res = c(res.dnds1, nonsyn1, syn1, g.norm1, res.dnds2, nonsyn2, syn2, g.norm2, res.cadd1, cadd.f.more1, cadd.f.less1, cadd.g.norm1, res.cadd2, cadd.f.more2, cadd.f.less2, cadd.g.norm2)
    names(res) = c("pub.dnds", "pub.nonsyn", "pub.syn", "pub.norm", "sub.dnds", "sub.nonsyn", "sub.syn", "sub.norm", "pub.dmfdlf", "pub.mf", "pub.lf", "pub.cadd.norm", "sub.dmfdlf", "sub.mf", "sub.lf", "sub.cadd.norm")
    return(res)
}


pycloneTable <- function(sampAB, samples, outdir, minAF = 0.05, minDepth = 30) {

    colnames = colnames(sampAB)
    keep = as.vector(apply(sampAB, 1, function(x, coln){
                         maxAF = max(as.numeric(x[match(paste(samples, "mafa", sep=""), coln)]))
                         minDP = min(as.numeric(x[match(paste(samples, "d", sep=""), coln)]))      
                         if (maxAF >= minAF & minDP >= minDepth) TRUE
                         else FALSE
                     },coln=colnames))
    sampAB = sampAB[keep,]
    colnames(sampAB) = colnames

    for (i in 1:length(samples)) {
        sn = samples[i]
        outTsv = paste(outdir, "/", sn, ".tsv", sep="")
        refci = match(paste(sn, "refc", sep=""), colnames(sampAB))
        altci = match(paste(sn, "altc", sep=""), colnames(sampAB))
        minor_cni = match(paste(sn, "nb", sep=""), colnames(sampAB))
        major_cni = match(paste(sn, "nt", sep=""), colnames(sampAB))
        
        mutation_id = paste(sampAB$chr, sampAB$pos, sep=":")
        refc = sampAB[,refci]
        altc = sampAB[,altci]
        minor_cn = round(sampAB[,minor_cni])
        major_cn = sampAB[,major_cni] - minor_cn
        otsv = data.frame(mutation_id=mutation_id, ref_counts=refc, var_counts=altc, normal_cn=2, minor_cn=minor_cn, major_cn=major_cn)
        write.table(otsv, file=outTsv, row.names=F, quote=F, sep="\t")
    }
}


scicloneTable <- function(sampAB, samples, outdir, minAF = 0.05, minDepth = 30, titanPath="./titan/") {

    colnames = colnames(sampAB)
    keep = as.vector(apply(sampAB, 1, function(x, coln) {
                         maxAF = max(as.numeric(x[match(paste(samples, "mafa", sep=""), coln)]))
                         minDP = min(as.numeric(x[match(paste(samples, "d", sep=""), coln)]))      
                         if (maxAF >= minAF & minDP >= minDepth) TRUE
                         else FALSE
                     },coln=colnames))
    sampAB = sampAB[keep,]
    colnames(sampAB) = colnames

    for (i in 1:length(samples)) {
        sn = samples[i]
        outTsv = paste(outdir, "/", sn, ".tsv", sep="")
        refci = match(paste(sn, "refc", sep=""), colnames(sampAB))
        altci = match(paste(sn, "altc", sep=""), colnames(sampAB))
        mafci = match(paste(sn, "mafc", sep=""), colnames(sampAB))
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        minor_cni = match(paste(sn, "nb", sep=""), colnames(sampAB))
        major_cni = match(paste(sn, "nt", sep=""), colnames(sampAB))
        
        mutation_id = paste(sampAB$chr, sampAB$pos, sep=":")
        refc = sampAB[,refci]
        altc = sampAB[,altci]
        mafc = sampAB[,mafci]*100
        mafa = sampAB[,mafai]*100
        minor_cn = round(sampAB[,minor_cni])
        major_cn = sampAB[,major_cni] - minor_cn

        #chr, pos, ref_reads, var_reads, vaf
        otsv = data.frame(chr=sampAB$chr, pos=sampAB$pos, ref_reads=refc, var_reads=altc, vaf=mafa)
        write.table(otsv, file=outTsv, row.names=F, quote=F, sep="\t")

        outSeg = paste(outdir, "/", sn, "_seg.tsv", sep="")
        cnv.inputA = paste(titanPath, sn, "_nclones1.TitanCNA.segments.txt", sep="")
        cnvA = read.delim(cnv.inputA)
        cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]
        #chr, start, stop, segment_mean
        otsg = data.frame(chr=cnvA$chrom, start=cnvA$loc.start, stop=cnvA$loc.end, segment_mean=cnvA$seg.mean)
        write.table(otsg, file=outSeg, row.names=F, quote=F, sep="\t")
    }
    
}


prepareLichee <- function(samples, nmaf, SampAB) {
    depthCols = paste(samples, "d", sep="")
    licheeCol = c("chr","pos","ref","alt","geneName",nmaf,paste(samples,"mafa",sep=""))
    licheeRow = !is.na(SampAB$functionalClass)
    for (i in 1:length(depthCols)) {  #Depth
        licheeRow = licheeRow & SampAB[,match(depthCols[i], colnames(SampAB))] >= 20
    }
    #for (i in 1:length(samples)) {  #LOH
    #    sn = samples[i]
    #   if (grepl("CRCTumor", sn)){
    #        next
    #    } else {
    #        nbCol = paste(samples[i], "nb", sep="")
    #        licheeRow = licheeRow & SampAB[,match(nbCol, colnames(SampAB))] != 0
    #    }
    #}
    licheeInput = SampAB[licheeRow,match(licheeCol,colnames(SampAB))]
    licheeInput = data.frame(licheeInput, name=paste(licheeInput$geneName,
                  paste(licheeInput$chr,licheeInput$pos,licheeInput$ref,sep=":"), licheeInput$alt, sep="_"))
    licheeInput = licheeInput[,match(c("chr","pos","name",nmaf,paste(samples,"mafa",sep="")),colnames(licheeInput))]
    licheeInput[,4] = 0
    #for (i in 5:dim(licheeInput)[2]) {
    #    licheeInput[,i] = as.numeric(licheeInput[,i])/2
    #}
    colnames(licheeInput) = gsub("mafa", "", colnames(licheeInput))
    colnames(licheeInput) = gsub("maf", "", colnames(licheeInput))
    return(licheeInput)
}


generateUpSetList <- function(sampAB, samples, minAF=0.08) {
    outlist = list()
    for(i in 1:length(samples)) {
        sn = samples[i]
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        needr = which(as.numeric(sampAB[,mafai]) >= minAF)
        outlist[[i]] = paste(sampAB$chr[needr],sampAB$pos[needr],sampAB$ref[needr],sampAB$alt[needr],sep=":")
    }
    names(outlist) = samples
    return(outlist)
}


spruceInput <- function(samp, samples, minMaf) {
    depthCols = paste(samples, "d", sep="")
    licheeRow = !is.na(samp$functionalClass)
    for (i in 1:length(depthCols)) {
        licheeRow = licheeRow & samp[,match(depthCols[i], colnames(samp))] >= 15
    }
    samp = samp[licheeRow,]
    result = data.frame()
    mutindex = 0
    for (r in 1:dim(samp)[1]) {    #each mutation
        maxMaf = 0
        for (i in 1:length(samples)) {
            sn = samples[i]
            mafi = match(paste(sn,"mafa",sep=""), colnames(samp))
            maf = as.numeric(samp[r,mafi])
            if (maf > maxMaf) {
                maxMaf = maf
            }
        }
        if (maxMaf < minMaf)  {
            next
        }
        mutn = as.character(samp[r, match("geneName", colnames(samp))])   #mutation name
        for (i in 1:length(samples)) {   #each tumor
            sn = samples[i]
            mafi = match(paste(sn,"mafc",sep=""), colnames(samp))
            depi = match(paste(sn,"d",sep=""), colnames(samp))
            nti = match(paste(sn,"nt",sep=""), colnames(samp))
            nbi = match(paste(sn,"nb",sep=""), colnames(samp))
            pui = match(paste(sn,"pu",sep=""), colnames(samp))
            segi = match(paste(sn,"seg",sep=""), colnames(samp))
            maf = as.numeric(samp[r,mafi])
            depth= as.numeric(samp[r,depi])
            nt = as.numeric(samp[r,nti])
            nb = round(as.numeric(samp[r,nbi]))
            na = nt-nb
            pu = as.numeric(samp[r,pui])
            nc = 1-pu
            seg = as.numeric(samp[r,segi])
            maflb = quantile(rbinom(1000,depth,maf),prob=seq(0,1,0.01))["1%"]/depth
            mafub = quantile(rbinom(1000,depth,maf),prob=seq(0,1,0.01))["99%"]/depth
            if (na == 1 & nb == 1){
                result = rbind(result, c(i-1, sn, mutindex, mutn, maflb, maf, mafub, 1, 1, 1, NA, NA, NA, seg),
                    stringsAsFactors=FALSE)
            } else {
                result = rbind(result, c(i-1, sn, mutindex, mutn, maflb, maf, mafub, 1, 1, nc, na, nb, pu, seg),
                    stringsAsFactors=FALSE)
            }
        } #each tumor
        mutindex = mutindex+1
    } #each mutation
    colnames(result) = c("sample_index","sample_label","character_label","character_index",
                "vaf_lb","vaf_mean","vaf_ub","x","y","mu","x","y","mu","seg")
    return(result)
}

#lod calculation
calLOD <- function(e,v,d) {
    f = round(v/d, 7)
    r = d-v
    Pr = f*(e/3) + (1-f)*(1-e)
    Pm = f*(1-e) + (1-f)*(e/3)
    Pr0 = 1-e
    Pm0 = e/3
    lod = log10((Pr^r*Pm^v)/(Pr0^r*Pm0^v))
    return(lod)
}


calNormalLOD <- function(e,v,d) {
    f = round(v/d, 7)
    fg = 0.5
    r = d-v
    Pr0 = 1-e
    Pm0 = e/3
    Prg = fg*(e/3) + (1-fg)*(1-e)
    Pmg = fg*(1-e) + (1-fg)*(e/3)
    lod = log10((Pr0^r*Pm0^v)/(Prg^r*Pmg^v))
    return(lod)
}
