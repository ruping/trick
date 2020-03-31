# plot multi sample mutation Table
script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir, "smkey.R", sep="/"))
library(caTools)
library(KernSmooth)
library(RColorBrewer)

getSampMutMulti <- function(samples, normal, d, cmedianTh, original, cmemeTh=6.5) {
    numsamples = length(samples)
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
    gerpindex = match("GERP_RS", colnames(d))
    siftindex = match("SIFT_score", colnames(d))
    polyphenindex = match("Polyphen2_HVAR_pred", colnames(d))
    somindex = match("somatic", colnames(d))
    repindex = match("rep",colnames(d))

    if (! is.na(caddindex)) {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,caddindex,gerpindex,siftindex,polyphenindex,dronindex,repindex,somindex,ncindex)]
    } else {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,dronindex,repindex,somindex,ncindex)]
    }
    
    resColnames = colnames(res)
    resAdd = t(data.frame(apply(res, 1, function(x, cindex, original, resColnames) {     #x is every row
              maxMaf = 0
              maxTlod = 0
              ssb = 0                    #combined strand bias check
              ssbc = 0
              ssbp = 0
              refdav = 0
              altdav = 0
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
                  sampleNameTime = paste(sampleName, "time", sep="")
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
                  strandBiases = strsplit(ss[[1]][4], ",")
                  strandBias = as.numeric(strandBiases[[1]][1])
                  strandBiasRef = 0
                  strandBiasFisherP = -1
                  if (length(strandBiases[[1]]) > 1){
                      strandBiasRef = as.numeric(strandBiases[[1]][2])
                      strandBiasFisherP = as.numeric(strandBiases[[1]][3])
                  }
                  mappingBias = as.numeric(ss[[1]][5])
                  cmedianav = as.numeric(cmeme[[1]][2])
                  cmemeSum = sum(as.numeric(cmeme[[1]]))
                  
                  #decide a b allele count
                  refnow = round(as.numeric(x[j+2])*(1-mafTmp))
                  altnow = round(as.numeric(x[j+2])*mafTmp)

                  if ( mafTmp > 0 ) {
                      ssb = ssb + strandBias*altnow
                      ssbp = ssbp + strandBiasFisherP*altnow
                      refdav = refdav + refnow
                      altdav = altdav + altnow
                      ssbc = ssbc + 1
                  }
                  
                  #decide mafNow
                  if ( mafTmp == 0 ) {
                      mafNow = mafTmp
                  } else if (endBias < 0.9 & ((strandBias != 0 & strandBias != 1) | (strandBiasFisherP > 0.7 & refnow >= 10 & altnow >= 5 & mafTmp >= 0.1)) & mappingBias < 0.8 & cmemeSum < cmemeTh & cmedianav < cmedianTh) {  
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
                  resVector = c(resVector, c(mafNow, 0, 0, 0, 0, refnow, altnow))
                  names(resVector)[(length(resVector)-6):length(resVector)] = c(sampleNameMaf,sampleNameMafa,sampleNameCcf,sampleNameCcfSD,sampleNameTime,sampleNameRef,sampleNameAlt)   #time added
              } #for each sample
              
              ssb = ssb/altdav
              ssbp = ssbp/altdav
              #message(paste(x[1],x[2], sep="  "))
              #message(paste(ssb, ssbp, sep="  "))
              altdav = altdav/ssbc
              refdav = refdav/ssbc
              #if (((ssb >= 0.95 | ssb <= 0.05) & ssbp < 0.1) | ssbp < 0.001) {                         #multiple sample strand bias
              if (((ssb > 0.88 | ssb < 0.12) & numsamples > 1) | ((ssb > 0.9 | ssb < 0.1) & numsamples == 1)) {                         #multiple sample strand bias
              #if ((ssb >= 0.95 | ssb <= 0.05)| ssbp < 0.001) {                                        #multiple sample strand bias 
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


getSampMutMultiIndel <- function(samples, normal, d) {
    numsamples = length(samples)
    rindexTF = vector()
    cindex = vector()
    for (si in 1:length(samples)) {
        if (length(rindexTF) > 0) {
            rindexTF = rindexTF | grepl(samples[si], d$traceSomatic)
        } else {
            rindexTF = grepl(samples[si], d$traceSomatic)
        }
        cindex = append(cindex, match(paste(samples[si], "maf", sep=""), colnames(d)))
    }
    rindex = which(rindexTF)
    cindex = as.vector(sapply(cindex, function(x){c(x-1,x,x+1)}))
    
    ncindex = match(paste(normal, "maf", sep=""), colnames(d))
    ncindex = as.vector(sapply(ncindex, function(x){c(x,x+1)}))  #normal samples maf and depth indexes
    
    #dronindex = match("dron", colnames(d))
    gnindex = match("geneName", colnames(d))
    glindex = match("geneLoc", colnames(d))
    gfindex = match("functionalClass", colnames(d))
    #caddindex = match("CADD_phred", colnames(d))
    #gerpindex = match("GERP_RS", colnames(d))
    #siftindex = match("SIFT_score", colnames(d))
    #polyphenindex = match("Polyphen2_HVAR_pred", colnames(d))
    somindex = match("traceSomatic", colnames(d))
    germindex = match("traceGermline", colnames(d))
    repindex = match("rep",colnames(d))
    scindex = match("sc",colnames(d))

    res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,repindex,scindex,somindex,germindex,ncindex)]
    
    resColnames = colnames(res)
    resAdd = t(data.frame(apply(res, 1, function(x, cindex, resColnames) {     #x is every row
              #ssb = 0                                                         #combined strand bias check
              ssbc = 0
              #ssbp = 0
              refdav = 0
              altdav = 0
              resVector = vector()
              for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample
                  sampleName = resColnames[j]
                  sampleNameMaf = paste(sampleName, "mafc", sep="")
                  sampleNameMafa = paste(sampleName, "mafa", sep="")
                  sampleNameCcf = paste(sampleName, "ccf", sep="")
                  sampleNameCcfSD = paste(sampleName, "ccfSD", sep="")
                  sampleNameTime = paste(sampleName, "time", sep="")
                  sampleNameAlt = paste(sampleName, "altc", sep="")
                  sampleNameRef = paste(sampleName, "refc", sep="")

                  altd = 0
                  refd = 0
                  tdep = 0
                  mafTmp = 0
                  cmeme = 0
                  endBias = 0
                  strandBias = 0
                  mappingBias = 0
                  if ( grepl("\\|", as.character(x[j])) ) {    #sample Name column
                      ori = strsplit(as.character(x[j]), "\\|")
                      oriDep = strsplit(ori[[1]][3], ",")
                      mafTmp = as.numeric(ori[[1]][1])
                      altd = as.numeric(oriDep[[1]][1])
                      tdep = as.numeric(oriDep[[1]][2])
                  } else if ( grepl("\\|", as.character(x[j+1])) ) { #samle Name maf column
                      ss = strsplit(as.character(x[j+1]), "\\|")
                      mafTmp = min(c(as.numeric(ss[[1]][1]),1))
                      cmeme = strsplit(ss[[1]][3], ",")
                      endBias = as.numeric(ss[[1]][2])
                      strandBias = as.numeric(ss[[1]][4])
                      mappingBias = as.numeric(ss[[1]][5])
                      tdep = as.numeric(x[j+2])
                      refd = max(c(round(tdep*(1-mafTmp)),0))
                      altd = round(tdep*mafTmp)
                  }
                  
                  #cmedianav = as.numeric(cmeme[[1]][2])
                  #cmemeSum = sum(as.numeric(cmeme[[1]]))

                  if ( mafTmp > 0 ) {
                      #ssb = ssb + strandBias*altnow
                      #ssbp = ssbp + strandBiasFisherP*altnow
                      refdav = refdav + refd
                      altdav = altdav + altd
                      ssbc = ssbc + 1
                  }
                  
                  #decide mafNow
                  mafNow = mafTmp
                  
                  resVector = c(resVector, c(mafNow, 0, 0, 0, 0, refd, altd))
                  names(resVector)[(length(resVector)-6):length(resVector)] = c(sampleNameMaf,sampleNameMafa,sampleNameCcf,sampleNameCcfSD,sampleNameTime,sampleNameRef,sampleNameAlt)   #time added
              } #for each sample
              
              #ssb = ssb/altdav
              #ssbp = ssbp/altdav
              #message(paste(x[1],x[2], sep="  "))
              #message(paste(ssb, ssbp, sep="  "))
              altdav = altdav/ssbc
              refdav = refdav/ssbc

              #if (((ssb >= 0.95 | ssb <= 0.05) & ssbp < 0.1) | ssbp < 0.001) {                         #multiple sample strand bias
              #if (((ssb > 0.88 | ssb < 0.12) & numsamples > 1) | ((ssb > 0.9 | ssb < 0.1) & numsamples == 1)) {                         #multiple sample strand bias
              #if ((ssb >= 0.95 | ssb <= 0.05)| ssbp < 0.001) {                                        #multiple sample strand bias 
              #    for (j in seq(6,(3+length(cindex)), by=3)) {
              #        sampleName = resColnames[j]
              #        sampleNameMaf = paste(sampleName, "mafc", sep="")
              #        resVector[sampleNameMaf] = 0
              #    }
              #}

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

          }, cindex=cindex, resColnames=resColnames)))
    

    res = cbind(res, resAdd)
    res = res[which(res$totalMaf != 0),]
    return(res)
}



#adjust CCF titan for multi samples
mergeCNA <- function(titanPath="../OS/titanresult/", sn, skipchunk = 19) {
    #prepare CNA
    cnv.inputA = paste(titanPath, sn, "_nclones1.TitanCNA.segments.txt", sep="")
    cnvA = read.delim(cnv.inputA)
    cnvA = cnvA[which(cnvA$num.mark > skipchunk),]                               #skip two few marks
    cp2 <- c(which(cnvA$logcopynumberratio[-1] != cnvA$logcopynumberratio[-nrow(cnvA)] |
                       cnvA$chrom[-1] != cnvA$chrom[-nrow(cnvA)] |
                           cnvA$LOHcall[-1] != cnvA$LOHcall[-nrow(cnvA)]),
             nrow(cnvA))
    cp1 <- c(1,cp2[-length(cp2)]+1)
    cnvA2 <- data.frame(chrom=cnvA$chrom[cp1],
                        loc.start=cnvA$loc.start[cp1],
                        loc.end=cnvA$loc.end[cp2],
                        num.mark=cnvA$num.mark[cp1],         #recal later
                        seg.mean=cnvA$seg.mean[cp1],         #recal later
                        copynumber=cnvA$copynumber[cp1],
                        minor_cn=cnvA$minor_cn[cp1],
                        major_cn=cnvA$major_cn[cp1],
                        allelicratio=cnvA$allelicratio[cp1], #recal later
                        LOHcall=cnvA$LOHcall[cp1],
                        cellularprevalence=cnvA$cellularprevalence[cp1],
                        ploidy=cnvA$ploidy[cp1],
                        normalproportion=cnvA$normalproportion[cp1],
                        logcopynumberratio=cnvA$logcopynumberratio[cp1])
    for (j in 1:length(cp1)) {
        cnvA2$num.mark[j] <- sum(cnvA$num.mark[cp1[j]:cp2[j]])
        cnvA2$seg.mean[j] <- mean(cnvA$seg.mean[cp1[j]:cp2[j]])
        cnvA2$allelicratio[j] <- mean(cnvA$allelicratio[cp1[j]:cp2[j]])
    }
    return(cnvA2)
}

cnaTiming <- function(sn="SPCG-OS052_11D", titanPath = "../OS/titanresult/", sampAB, public=FALSE, skipchunk = 19,
                      mmut=10, method="fullMLE") {

    #mergeCNA
    cnvA2 = mergeCNA(titanPath=titanPath, sn=sn, skipchunk = skipchunk)

    #overlap between CNA and SNV
    cnvSeqNames = cnvA2$chrom
    if ( !grepl("chr", cnvA2$chrom[1]) ) {
        cnvSeqNames = paste("chr", cnvA2$chrom, sep="")
    }
    snvSeqNames = sampAB$chr
    if ( !grepl("chr", sampAB$chr[1]) ) {
        snvSeqNames = paste("chr", sampAB$chr, sep="")
    }
    cnvRangeA = GRanges(seqnames = cnvSeqNames, ranges = IRanges(cnvA2$loc.start, end=cnvA2$loc.end), strand=rep('+',dim(cnvA2)[1]))
    snvRange = GRanges(seqnames = snvSeqNames, ranges = IRanges(sampAB$pos, end=sampAB$pos), strand=rep('+',dim(sampAB)[1]))
    foA = findOverlaps(cnvRangeA, snvRange)
    cnvHits = queryHits(foA)
    snvHits = subjectHits(foA)

    #produce results
    result = list()
    resultTable = vector()
    li = 1
    for (i in 1:dim(cnvA2)[1]) {
        cnchrom = cnvA2$chrom[i]
        cnstart = cnvA2$loc.start[i]
        cnend = cnvA2$loc.end[i]
        cnLOHcall = as.character(cnvA2$LOHcall[i])
        cnminor = cnvA2$minor_cn[i]
        cnmajor = cnvA2$major_cn[i]
        cntotal = cnvA2$copynumber[i]
        if (cntotal > 4 | cntotal < 2 | (cnminor == 1 & cnmajor == 1)){
            #message(paste("skip cnv", i, "due to cnv configuration", cnmajor, cnminor, sep=" "))
            next
        }
        if (length(which(cnvHits == i)) < mmut) {
            #message(paste("skip cnv", i, sep=" "))
            next     #skip too few mutations 1
        }
        
        chopped = sampAB[snvHits[which(cnvHits == i)],]
        if (public == TRUE) {
            chopped = chopped[which(chopped$pubOrSub == "public"),]
        } else {
            chopped = chopped[which(chopped[,paste(sn,"ccf",sep="")] + 2.58*chopped[,paste(sn,"ccfSD",sep="")] >= 1),]
        }
        if (dim(chopped)[1] < mmut) {
            #message(paste("skip cnv", i, "second round", dim(chopped)[1], sep=" "))
            next    #skip too few mutations 2
        }
        normCont = 1-chopped[,paste(sn,"pu",sep="")][1]
        
        onlyMuts = data.frame(chromosome=chopped$chr, position=chopped$pos, refbase=chopped$ref, mutbase=chopped$alt,
            rsID=chopped$id, t_ref_count=chopped[,paste(sn,"refc",sep="")], t_alt_count=chopped[,paste(sn,"altc",sep="")],
            allelefreq=chopped[,paste(sn,"mafc",sep="")],
            t_depth=chopped[,paste(sn,"refc",sep="")] + chopped[,paste(sn,"altc",sep="")])
        
        if (cnmajor == 2 & cnminor == 0) {
            hmatrix = makeEventHistory(type="LOH",copies=c(cnmajor,cnminor),onlyIdentifiable=F)[[1]]
            cntype = "CNLOH"
        } else if (cnmajor == 2 & cnminor == 2) {   #convert to CNLOH type
            hmatrix = matrix(c(0,2,1,0),ncol=2,nrow=2,byrow=TRUE)
            cntype = "CNLOH"
            cntotal = 2
            onlyMuts$allelefreq = sapply(onlyMuts$allelefreq, function(x){ min(1, 2*x)})
            onlyMuts$t_alt_count = round(onlyMuts$allelefreq*onlyMuts$t_depth)
            onlyMuts$t_ref_count = onlyMuts$t_depth - onlyMuts$t_alt_count
        } else if (cnmajor == 3 & cnminor == 0) {
            hmatrix = matrix(c(1,2,3,0,1,0,1,0,0), byrow=T, ncol=3)
            cntype = "ALOH"
        } else if (cnmajor == 4 & cnminor == 0) {
            hmatrix = matrix(c(1,2,3,4,0,0,1,0,0,1,0,0,1,0,0,0), byrow=T, ncol=4)
            cntype = "ALOH"
        } else {
            hmatrix = makeEventHistory(type="gain",copies=c(cnmajor,cnminor),onlyIdentifiable=F)[[1]]
            cntype = "gain"
        }
        
        x<-eventTiming(x=onlyMuts$t_alt_count, m=onlyMuts$t_depth, history=hmatrix, totalCopy=cntotal, minMutations=mmut,
                       type=cntype,normCont=normCont, bootstrapCI="parametric",CILevel=0.9, method=method)
        x = c(x, cnid=i, cnchrom=cnchrom, cnstart=cnstart, cnend=cnend, cnmajor=cnmajor, cnminor=cnminor, cnLOHcall=cnLOHcall)

        result[[li]] = x
        names(result)[[li]] = i
        currentline = c(cnchrom, cnstart, cnend, cnmajor, cnminor, x$summaryTable[1],
            as.numeric(x$pi[1]), x$piCI[1,1], x$piCI[1,2], as.numeric(x$pi[2]), x$piCI[2,1], x$piCI[2,2])
        names(currentline) = c("chrom","loc.start","loc.end","major_cn","minor_cn","nmut",
                 "p0","p0l","p0h","p1","p1l","p1h")
        message(paste(currentline, collapse=" "))
        if (length(resultTable) == 0) {
            resultTable = currentline
        } else {
            resultTable = rbind(resultTable, currentline)
        }
        li=li+1
    }
    resultTable = data.frame(resultTable)
    return(list(result, resultTable, cnvA2))    
}



adjust.ccf.titan.multi <- function(sampAB, samples, t, titanPath="./titan/", correctColname=FALSE, overadj=1.6, sigTh=0.9, skipchunk=19) {
    
    numsamples = length(samples)
    purities = vector()
    
    for (i in 1:length(samples)) {
        sn = samples[i]
        message(sn)
        cnv.inputA = paste(titanPath, sn, "_nclones1.TitanCNA.segments.txt", sep="")
        cnvA = read.delim(cnv.inputA)
        cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]                 #skip NA
        cnvA = cnvA[which(cnvA$num.mark > 9),]                               #skip two few marks
        cnvA$nt = cnvA$copynumber
        if ("minor_cn" %in% colnames(cnvA)) {
            cnvA$nb = cnvA$minor_cn
        } else { 
            cnvA$nb = partialRound(cnvA$copynumber*(
                1-((cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion) - (1-cnvA$cellularprevalence)*0.5)
                /cnvA$cellularprevalence))
        }
        
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
        
        #pa,pu,nt,nb,seg
        message("info table building")
        queHits = queryHits(foA)
        subHits = subjectHits(foA)
        infos = t(sapply(1:dim(sampAB)[1], function(x, queHits, subHits) {                           
                             if (x %in% queHits){
                                 queMIndex = match(x, queHits)
                                 subMIndex = subHits[queMIndex]
                                 c(cnvA$cellularprevalence[subMIndex],     #pa
                                   cnvA$nt[subMIndex],                     #nt
                                   cnvA$nb[subMIndex],                     #nb
                                   subMIndex)                              #seg
                             } else {
                                 c(0,2,1,0)
                             }}, queHits = queHits, subHits = subHits))
        infos = data.frame(infos)
        message("info table built")

        pa1 = infos[,1]
        nt1 = infos[,2]
        nb1 = infos[,3]
        seg1 = infos[,4]
        con1 = as.numeric(cnvA$normalproportion[1])          #pu1 = 1-cnvA$normalproportion[1]
        pu1 = 1 - (con1 + (1-max(as.numeric(pa1)))*(1-con1))
        sAGP = pa1*(1-con1)
        message(pu1)
        purities = append(purities, pu1)
        names(purities)[length(purities)] = sn
        
        sampAB = data.frame(sampAB, pu=pu1, pa=pa1, sAGP=sAGP, nt=nt1, nb=nb1, seg=seg1)
        colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]] = paste(sn,colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]], sep="")
        if ( correctColname == TRUE ) {
            colnames(sampAB) = gsub("\\.","-",colnames(sampAB))
        }
    }
    
    
    for( i in 1:dim(sampAB)[1]) {  # rescale the maf and calculate CCF

        if (i %% 1000 == 0) {
            message(i)
        }
        
        foundSites = 0             # count how many sites found
        depthTotal = 0
        mafaTotal = 0
        ccfTotal = 0
        ccfsdTotal = 0
        for (j in 1:length(samples)) {
            sn = samples[j]
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            if (maf1 > t) {
                foundSites = foundSites+1
            }
        }
        
        for (j in 1:length(samples)) {
            sn = samples[j]

            pa1 = as.numeric(sampAB[i, match(paste(sn, "pa", sep=""), colnames(sampAB))])
            nt1 = as.numeric(sampAB[i, match(paste(sn, "nt", sep=""), colnames(sampAB))])
            nb1 = as.numeric(sampAB[i, match(paste(sn, "nb", sep=""), colnames(sampAB))])
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            refc1 = as.numeric(sampAB[i, match(paste(sn, "refc", sep=""), colnames(sampAB))])
            altc1 = as.numeric(sampAB[i, match(paste(sn, "altc", sep=""), colnames(sampAB))])
            pu1 = as.numeric(sampAB[i, match(paste(sn, "pu", sep=""), colnames(sampAB))])       #cell purity
            #if (nt1 > 0) pu1 = nt1*pu1/(nt1*pu1+2*(1-pu1))                                     #effective purity
            sAGP = as.numeric(sampAB[i, match(paste(sn, "sAGP", sep=""), colnames(sampAB))])    #segmental aneu- ploid genome proportion
            if (nt1 > 0 & nt1 != 2) pu1 = nt1*sAGP/(nt1*sAGP+2*(1-sAGP))                                   #effective purity
            
            if (maf1 > 0) {
                if ((maf1 > t & foundSites >= 2) | numsamples == 1) {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"unknown",overadj=overadj,sigTh=sigTh)
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                    sampAB[i, match(paste(sn, "time", sep=""), colnames(sampAB))] = CCF1[2]   #need to change
                } else {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"late",overadj=overadj,sigTh=sigTh)
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                    sampAB[i, match(paste(sn, "time", sep=""), colnames(sampAB))] = CCF1[2]
                }
            }

            #for merged MAF
            dep1 = as.numeric(sampAB[i, match(paste(sn, "d", sep=""), colnames(sampAB))])
            mafa1 = as.numeric(sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))])
            ccf1 = as.numeric(sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))])
            ccfsd1 = as.numeric(sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))])
            depthTotal = depthTotal + dep1
            mafaTotal = mafaTotal + dep1*mafa1
            ccfTotal = ccfTotal + dep1*ccf1
            ccfsdTotal = ccfsdTotal + dep1*ccfsd1
            #for merged MAF
            
        }     #for each sample
        sampAB[i, match("mergeMAFA", colnames(sampAB))] = round(mafaTotal/depthTotal, 5)
        sampAB[i, match("mergeCCF", colnames(sampAB))] = round(ccfTotal/depthTotal, 5)
        sampAB[i, match("mergeCCFsd", colnames(sampAB))] = round(ccfsdTotal/depthTotal, 5)
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

plotRes.multi.matrix.pdf <- function(sampAB, sroot, samples, pdfsize = 16, plotType = "AF", snr="", sns=vector(),ssAF=0, pobs=rep("pubOrSub",dim(combn(length(samples),2))[2]), png=FALSE) {
    combinations = combn(length(samples),2)
    nplots = dim(combinations)[2]
    ndims = length(samples)-1

    resStats = list()
    outfile = paste(sroot,"multi_hist.pdf",sep="")
    if (plotType == "Scatter") {
        outfile = paste(sroot,"multi_scatter.pdf",sep="")
    }
    if (plotType == "Density") {
        outfile = paste(sroot,"multi_density.pdf",sep="")
    }

    if (png){
        png(file = outfile, width=pdfsize*100, height=pdfsize*100)
    } else {
        pdf(file = outfile, width=pdfsize, height=pdfsize)
    }
    layout(matrix(seq(1,ndims^2), ndims, ndims, byrow = FALSE))
    
    pindex = 0
    for (ci in 1:ndims) {                           #for each column
        noplotRow = vector()
        if (ci > 1) {
            noplotRow = (1:ndims)[1:(ci-1)]         #which row do not plot?
        }
        for (ri in 1:ndims) {                       #for each row
            if (ri %in% noplotRow) {                #no plot
                plot(NULL,NULL,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
            } else {                                #plot
                pindex = pindex + 1
                pair = combinations[,pindex]
                sn1 = samples[pair[1]]              #samplename1
                sn1mafa = paste(sn1, "mafa",sep="")
                sn1d = paste(sn1, "d",sep="")
                sn1short = gsub("CRCTumor","",sn1)
                sn1short = gsub("Recurrence","Rec",sn1short)
                sn1short = gsub("Primary","Pri",sn1short)
                sn1short = gsub(paste(sroot,"_",sep=""),"",sn1short)
                sn1short = gsub(paste(sroot,"-",sep=""),"",sn1short)
                sn2 = samples[pair[2]]              #samplename2
                sn2mafa = paste(sn2, "mafa",sep="")
                sn2d = paste(sn2, "d",sep="")
                sn2short = gsub("CRCTumor","",sn2)
                sn2short = gsub("Recurrence","Rec",sn2short)
                sn2short = gsub("Primary","Pri",sn2short)
                sn2short = gsub(paste(sroot,"_",sep=""),"",sn2short)
                sn2short = gsub(paste(sroot,"-",sep=""),"",sn2short)
                main.title = paste(snr, sns[pair[1]], "vs", sns[pair[2]], sep=" ")
                statName = paste(sroot, sn1short, sn2short, "stats", sep="_")
                if (plotType == "AF") {
                    resStats[[statName]] = plotRes.multi.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], sn1mafa, sn2mafa, 0.05,
                                plotDepth=F,plotDensity=F,plotScatter=F, pdf=F,ssAF=ssAF,pob=pobs[pindex])
                } else if (plotType == "Scatter") {
                    resStats[[statName]] = plotRes.multi.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], sn1mafa, sn2mafa, 0.05,
                                plotDepth=F,plotDensity=F,plotAF=F,pdf=F,ssAF=ssAF,pob=pobs[pindex])
                } else if (plotType == "Density") {
                    resStats[[statName]] = plotRes.multi.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], sn1mafa, sn2mafa, 0.05,
                                plotDepth=F,plotScatter=F,plotAF=F,pdf=F,ssAF=ssAF,pob=pobs[pindex])
                }
            }
        }
    }
    dev.off()
    return(resStats)
}



plotRes.multi.pdf <- function(sampAB, sampName, main=sampName, sn1n="", sn2n="", sn1, sn2, minAF, ratio=1, plotAF=TRUE, plotDensity=TRUE, plotDepth=TRUE, plotScatter=TRUE, pdf=TRUE,
                              alpha=1, binw=0, widthadj=0, heightadj=0, nohistlegend = FALSE, ssAF=0, pob="pubOrSub") {

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

  #For recording AUC
  AUCstuff = subclonalMut(sampAB, sn1s, sn2s, minAF, ssAF=ssAF, pob=pob)
  
  #check depth power to reject a presence of a mutation
  depthPowerKeep <- as.vector(apply(sampAB, 1, function(x,mafa1i,mafa2i,dp1i,dp2i) {
                                        if(as.numeric(x[mafa1i]) == 0){vaf = as.numeric(x[mafa2i])
                                                    if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.numeric(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else if (as.numeric(x[mafa2i]) == 0){vaf = as.numeric(x[mafa1i])
                                                    if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.numeric(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else {TRUE}
                                    }, mafa1i=mafa1Index,mafa2i=mafa2Index,dp1i=dp1i,dp2i=dp2i))
  sampAB = sampAB[depthPowerKeep,]

  subMuts = subclonalMut(sampAB, sn1s, sn2s, minAF, ssAF=ssAF, pob=pob)       #subclonal mutations
  subMuts$rAUC = AUCstuff$rAUC
  subMuts$weightAF = AUCstuff$weightAF
  subMuts$subArow = AUCstuff$subArow
  subMuts$subBrow = AUCstuff$subBrow
  
  #allA_Rows = which(sampAB[,maf1Index] > minAF & sampAB[,mafa1Index] > minAF & sampAB[,maf1Index] <= 1)
  allA_Rows = intersect(union(subMuts$subAi, subMuts$pubTi), which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))
  #message(paste(subMuts$subArow[match(allA_Rows[which(sampAB[allA_Rows, maf1Index] <= 0.05)], subMuts$subAi)], collapse="  "))
  subA_Rows = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))
  ssA_Rows  = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1 & sampAB[,maf2Index] <= ssAF))
  BinWidthA = round(dpih(sampAB[allA_Rows, maf1Index]/ratio),2)
  if ((BinWidthA < 0.02 & length(allA_Rows) < 3000) | BinWidthA == 0) { BinWidthA = 0.02 }
  
  #allB_Rows = which(sampAB[,maf2Index] > minAF & sampAB[,mafa2Index] > minAF & sampAB[,maf2Index] <= 1)
  allB_Rows = intersect(union(subMuts$subBi, subMuts$pubTi), which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))
  allAB_Rows = union(allA_Rows,allB_Rows)
  subB_Rows = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))
  ssB_Rows  = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1 & sampAB[,maf1Index] <= ssAF))
  
  BinWidthB = round(dpih(sampAB[allB_Rows, maf2Index]/ratio),2)
  if ((BinWidthB < 0.02 & length(allB_Rows) < 3000) | BinWidthB == 0) { BinWidthB = 0.02 }
  BinWidth = min(c(BinWidthA, BinWidthB, 0.1))
  if ( binw != 0 ) {
      BinWidth = binw
  }
  message(paste("bin width: ", BinWidthA, BinWidthB, BinWidth, sep=" "))
  nbreaksA = round((max(sampAB[allA_Rows, maf1Index]/ratio)-min(sampAB[allA_Rows, maf1Index]/ratio)+0.01)/BinWidth)
  nbreaksB = round((max(sampAB[allB_Rows, maf2Index]/ratio)-min(sampAB[allA_Rows, maf2Index]/ratio)+0.01)/BinWidth)
  nbreaks = ceiling(diff(range(minAF,1))/BinWidth)
  message(paste("nbreaks: ", nbreaksA, nbreaksB, nbreaks, sep=" "))
  breaksA = seq(minAF,1,length.out=nbreaks)
  breaksB = seq(minAF,1,length.out=nbreaks)
  
  
  sampAh = hist(sampAB[allA_Rows, maf1Index]/ratio, breaks=breaksA,  plot=F)
  sampAhsub = hist(sampAB[subA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  sampAhss = hist(sampAB[ssA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  ylimup = max(sampAh$count)
  
  sampBh = hist(sampAB[allB_Rows, maf2Index]/ratio, breaks=breaksB, plot=F)
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
  
  if (plotAF == TRUE) {

      if (pdf == TRUE) {
          pdf(file = paste(sampName, "hist.pdf", sep="_"), width = 8+widthadj, height = 8+heightadj, useDingbats=FALSE)
      }
      par(mar=c(4.5,5,4.5,0))

      plot( sampAh, col=rgb(0,0,0,1/4), xlim=c(0, 1), ylim=c(ylimdown,ylimup), border=F, ylab="# of Mutations", xlab="VAF", axes = F, main = main, cex.lab = 2.3, cex.main = 2.3)  # first histogram
      plot( sampAhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green set border
      plot( sampAhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
      
      plot( sampBh,col=rgb(0,0,0,1/4), border=F, add=T )  # second histogram
      plot( sampBhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green
      plot( sampBhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
      sampAhss$counts = 0                                  #make black line
      sampBhss$counts = 0                                  #make black line
      plot( sampAhss, col="black", add=T, border=F)        #make black line
      plot( sampBhss, col="black", add=T, border=F)        #make black line

      sn1s = sn1n
      if (sn1s == ""){
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
      }
      sn2s = sn2n
      if (sn2s == ""){
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
      }
      message(sn1s)
      message(sn2s)
      
      axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=1.7)
      axis(side=2,at=decideTickAt(ylimdown, ylimup),labels=allAbs(decideTickAt(ylimdown, ylimup)),cex.axis=1.7)
      text(x=0.5,y=ylimdown,labels=paste(sn2s,",",sum(subMuts$pubTn,subMuts$sharedBn,subMuts$ssBn), "SSNVs", sep = " "), cex=2.2)
      text(x=0.5,y=ylimup,labels=paste(sn1s,",",sum(subMuts$pubTn,subMuts$sharedAn,subMuts$ssAn), "SSNVs", sep = " "), cex=2.2)

      #stats starting from here!
      fHsub = round(mean(c(subMuts$ratioHighSubA,subMuts$ratioHighSubB)),3)
      fst = round(subMuts$FST,3)
      ksd = round(subMuts$KSD,3)
      text(x=0.8,y=3*ylimup/4, labels=bquote(paste("fH"["sub"], " = ", .(fHsub))),cex=2.2)
      text(x=0.8,y=(3/4-0.167)*ylimup, labels=bquote(paste("FST", " = ", .(fst))),cex=2.2)
      text(x=0.8,y=(3/4-0.334)*ylimup, labels=bquote(paste("KSD", " = ", .(ksd))),cex=2.2)

      npub = subMuts$pubTn
      legendText = c(paste("Public ","(",npub,")",sep=""),"Pvt-Shared","Pvt-Rgn Specific")
      legendXpos = 0.57
      legendYpos = ylimdown/3
      lengedCol = c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1))
      if (nohistlegend == TRUE){
          legendText = c(paste("",npub,sep=""))
          legendXpos = 0.75
          legendYpos = ylimdown/2
          lengedCol = c(rgb(0,0,0,1/4))
      }
      legend(legendXpos,ylimdown/3, legend=legendText, col=lengedCol, pch=15, bty="n", cex=1.7)
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotDensity == TRUE) {
      #two way density plot
      if (pdf == TRUE) {
          pdf(file = paste(sampName, "density.pdf", sep="_"), width = 8+widthadj, height = 8+heightadj, useDingbats=FALSE)
      }
      par(mar=c(4.5,5,4.5,2))
      if (pdf == FALSE) {
          par(mar=c(5,5,5,4))
      }
      smkey(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("VAF",sn1s,sep=" "), ylab=paste("VAF",sn2s,sep=" "),
            main = main, xlim=c(0,1), ylim=c(0,1),cex.lab = 2.3, cex.main = 2.3,cex.axis=1.7)
      if (length(which(sampAB$dron != 0)) > 0) {
          dronindex = which(sampAB$dron != 0)
          pointLabel(sampAB[dronindex,maf1Index],sampAB[dronindex,maf2Index],labels=as.character(sampAB$geneName[dronindex]),
                     col=decideDriverColor(sampAB[dronindex,maf1Index],sampAB[dronindex,maf2Index]),cex=1.7,font=2)
      }
      segments(0,0,0.8,0.8)
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotScatter == TRUE) {
      #scatter plot with density color
      if (pdf == TRUE) {
          pdf(file = paste(sampName, "scatter.pdf", sep="_"), width = 5 + widthadj/2, height = 5 + heightadj/2, useDingbats=FALSE)
      }
      drx = vector()
      dry = vector()
      drlabels = vector()
      if (length(which(sampAB$dron == 1)) > 0) {
          dronindex = which(sampAB$dron != 0)
          knownDronIndex = which(sampAB$dron == 2)
          knownDrivers = as.character(sampAB$geneName[knownDronIndex])          #only for known drivers
          dronindex = intersect(allAB_Rows, dronindex)
          drx = sampAB[dronindex,maf1Index]
          dry= sampAB[dronindex,maf2Index]
          drlabels = as.character(sampAB$geneName[dronindex])
          drlabels = as.vector(sapply(drlabels, function(x, knownDrivers){if (x %in% knownDrivers){x} else {""}}, knownDrivers=knownDrivers))
      }
      pub_Rows = setdiff(allAB_Rows, union(subA_Rows, subB_Rows))
      pub_Rows = match(pub_Rows, allAB_Rows)
      shared_Rows = setdiff(union(subA_Rows, subB_Rows), union(ssA_Rows, ssB_Rows))
      shared_Rows = match(shared_Rows, allAB_Rows)
      ssAB_Rows = union(ssA_Rows, ssB_Rows)
      ssAB_Rows = match(ssAB_Rows, allAB_Rows)
      allSub_Rows = match(union(subA_Rows, subB_Rows), allAB_Rows)
      if (pdf == TRUE) {
          scatterDensityPlot(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("VAF",sn1s,sep=" "), ylab=paste("VAF",sn2s,sep=" "),
                             main = main, cex=1.2, cex.lab = 2.3, cex.main = 2.3, cex.axis=1.7, drx=drx, dry=dry, drlabels = drlabels,
                             groups=list(a=pub_Rows, b=shared_Rows, c=ssAB_Rows),alpha=alpha,
                             groupColors=list(a=brewer.pal(9, "Greys")[3:9], b=brewer.pal(9, "Greens")[2:5], c=brewer.pal(9, "Blues")[3:9]))
      } else {
          scatterDensityPlot(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("VAF",sn1s,sep=" "), ylab=paste("VAF",sn2s,sep=" "),
                             main = main, cex=1.2, cex.lab = 2.3, cex.main = 2.3, cex.axis=1.7, drx=drx, dry=dry, drlabels = drlabels, layout = FALSE,
                             groups=list(a=pub_Rows, b=shared_Rows, c=ssAB_Rows),alpha=alpha,
                             groupColors=list(a=brewer.pal(9, "Greys")[3:9], b=brewer.pal(9, "Greens")[2:5], c=brewer.pal(9, "Blues")[3:9]))
      }
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotDepth == TRUE) {
      maint = gsub("\\s.+","",main)
      pdf(file = paste(sampName, "depth.pdf", sep="_"), width = 10, height = 6.6, useDingbats=FALSE)
      par(mfrow=c(2,3))
      ymaxdepth = max(quantile(sampAB[which(sampAB[,mafa1Index] > minAF),dp1i], prob=seq(0,1,0.05))["95%"],
          quantile(sampAB[which(sampAB[,mafa1Index] > minAF),dp2i], prob=seq(0,1,0.05))["95%"])
      calledA = which(sampAB[,mafa1Index] > minAF & sampAB[,mafa1Index] <= 1)
      calledB = which(sampAB[,mafa2Index] > minAF & sampAB[,mafa2Index] <= 1)
      plot(sampAB[calledA,mafa1Index],sampAB[calledA,dp1i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("VAF",sn1s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste(maint,",","sSNVs in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      legend("topright",legend=c("Pvt-Site Specific","Pvt-Shared","Public"),col=c(rgb(31/255,120/255,180/255,1),rgb(178/255,223/255,138/255,1),rgb(0,0,0,1/4)),pch=19, cex=1.3)
      plot(sampAB[calledA,mafa1Index],sampAB[calledA,dp2i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("VAF",sn1s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste("sSNVs in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledA,dp1i],sampAB[calledA,dp2i],xlim=c(0,ymaxdepth),ylim=c(0,ymaxdepth),col = subSSColor(calledA, subA_Rows, ssA_Rows),
           xlab=paste("depth",sn1s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste("sSNVs in",sn1s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,mafa2Index],sampAB[calledB,dp2i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("VAF",sn2s,sep=" "), ylab=paste("depth",sn2s,sep=" "), main=paste(maint,",","sSNVs in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,mafa2Index],sampAB[calledB,dp1i],xlim=c(0,1),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("VAF",sn2s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste("sSNVs in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      plot(sampAB[calledB,dp2i],sampAB[calledB,dp1i],xlim=c(0,ymaxdepth),ylim=c(0,ymaxdepth),col = subSSColor(calledB, subB_Rows, ssB_Rows),
           xlab=paste("depth",sn2s,sep=" "), ylab=paste("depth",sn1s,sep=" "), main=paste("sSNVs in",sn2s,sep=" "),cex.main=1.5,cex.lab=1.5)
      dev.off()
  }

  return(subMuts)
  #return(ssfs)
}


subclonalMut <- function(sampAB, snA, snB, minAF=0.08, statsAF=0.08, highAF=0.2, ratio=1, crv=2.58, ssAF=0, pob="pubOrSub")  {                   #determinine subclonal mutations
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
    pobi = match(pob, colnames(sampAB))

    gNIndex = match("geneName",colnames(sampAB))
    gLIndex = match("geneLoc",colnames(sampAB))
    fCIndex = match("functionalClass",colnames(sampAB))
    
    # subclonal mutations
    subAi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > minAF & ((sampAB[,mafaBi] == 0 & (sampAB[,nbBi] != 0 | sampAB[,nbAi] == 0)) | sampAB[,mafaBi] != 0) )
    subArow = rownames(sampAB)[subAi]
    mutsA = sampAB[subAi,mafaAi]/ratio
    ssAi  = intersect(subAi, which( sampAB[,mafaAi] > minAF & sampAB[,mafaBi] <= ssAF ))
    
    subBi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaBi] > minAF & ((sampAB[,mafaAi] == 0 & (sampAB[,nbAi] != 0 | sampAB[,nbBi] == 0)) | sampAB[,mafaAi] != 0) )
    subBrow = rownames(sampAB)[subBi]
    mutsB = sampAB[subBi,mafaBi]/ratio
    ssBi  = intersect(subBi, which( sampAB[,mafaBi] > minAF & sampAB[,mafaAi] <= ssAF ))
    
    KSD = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)
    
    #for rAUC
    allSubRows = union(subAi,subBi)
    mafs = paste(c(snA,snB), "mafa", sep="")
    depths = paste(c(snA,snB), "d", sep="")
    rAUCout = rAUC(sampAB[allSubRows,], 0.04, mafs, depths)     #need to add for lower bound
    rAUC = round(rAUCout$rAUC,8)
    weightAF = rAUCout$weightAF

    # for PUB
    pubTi = which( sampAB[,pobi] == "public" &
                      !((sampAB[,mafaAi] < 0.15 & sampAB[,nbAi] == 0) | (sampAB[,mafaBi] < 0.15 & sampAB[,nbBi] == 0)) )
    pubTrow = rownames(sampAB)[pubTi]

    # for FST
    mutsSub = sampAB[allSubRows,]
    mutsSub = data.frame( maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi] )
    FST = fst.hudson(mutsSub, minAF=statsAF)

    
    # for other stats

    mutsA2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF )), mafaAi]
    mutsAh2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF )), mafaAi]
    mutsASr2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] > 0.02)), mafaAi]   #shared
    
    mutsASp2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] == 0)), mafaAi]
    mutsASph2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF & sampAB[,mafaBi] == 0)), mafaAi]
    
    mutsB2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF )), mafaBi]
    mutsBh2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF )), mafaBi]
    mutsBSr2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] > 0.02)), mafaBi]  #shared
    
    mutsBSp2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] == 0)), mafaBi]
    mutsBSph2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF & sampAB[,mafaAi] == 0)), mafaBi]


    # mutational counts
    pubTn = length(pubTi)
    ssAn = length(ssAi)
    ssBn = length(ssBi)
    sharedAn = length(mutsASr2)
    sharedBn = length(mutsBSr2)
    #sharedn = length(union(subAi, subBi))-length(union(ssAi,ssBi))
        
    # list for output    
    muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, ssAi=ssAi, ssBi=ssBi, pubTi=pubTi, subArow=subArow, subBrow=subBrow, pubTrow = pubTrow,
        lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
        lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
        lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
        lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
        lenSharedA=length(mutsASr2), lenSharedB=length(mutsBSr2), ratioSharedA=length(mutsASr2)/length(mutsA2), ratioSharedB=length(mutsBSr2)/length(mutsB2),
        FST=FST, KSD=KSD, rAUC=rAUC, weightAF=weightAF, pubTn=pubTn, sharedAn=sharedAn, sharedBn=sharedBn, ssAn=ssAn, ssBn=ssBn)
    return(muts)
}


subSSColor <- function(allindex, subindex, ssindex) {
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
    for (i in 1:length(x)) {
        x[i] = abs(x[i])
    }
    return(x)
}


computeCCF <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90) {

    ccf = 0
    ccf2 = 0
    sd = 0
    cc <- seq(0.02, 1, by = 0.01)
    evoType = "A1/A2/B/C"
    N = A + S
    nc = nt * pa + 2 * (1 - pa)
    nc2 = nt * sAGP + 2 * (1 - sAGP)

    #message(nc)
    #message(paste(c(f,A,S,pu,pa,sAGP,nt,nb,nc),collapse=" "))
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
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        if (Ptot > 0) {
            cp.A <- pEarly.a/Ptot
            cp.CD <- 1 - cp.A
            cp.C <- pLate/Ptot
            cp.D <- pEuploid/Ptot
            cp.AC <- 1 - cp.D
            cp.AD <- 1 - cp.C
                                        #message(paste(c(Ptot), collapse=" "))
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
        }
    } else if (nb == 0 & nt == 0) {        #homozygous deletion
        evoType <- "C"
        ccf = (f*nc2)/sAGP
        ff.C = cc[cc<=(1-sAGP)]/nc2
        Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=(1-sAGP)])
        if (sAGP > 0.98) {
            ff.C = min(1, cc[cc<=0.02]/nc2)
            Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=0.02])
            message(paste(ff.C, collapse="   "))
        }
        ccf2 <- (Ms.C$M1)/pu                    #dbinom
        
        sd <- Ms.C$SD
        if (is.nan(sd)) {
            sd = 0.01
        }
        #message(paste(c(ccf,ccf2,sd,f),collapse=" "))
    } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cpEarly.a <- pEarly.a/Ptot
        cpLate.eup <- 1 - cpEarly.a
        cpLate <- pLate/Ptot
        cpEup <- pEuploid/Ptot
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
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
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
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2     #present in major
        fl.ea <- (sAGP * (nt - nb))/nc2                
        fh.eb <- (nb * sAGP + 1 - sAGP)/nc2            #present in minor
        fl.eb <- nb * sAGP/nc2

        pEarly.ib = vector()
        if ((nt-nb-1) > nb) {                 #present in between gains
            for (ib in (nb+1):(nt-nb-1)) {    #for each state, low to high CN
                fh.ebt = (ib * sAGP + 1 - sAGP)/nc2
                fl.ebt = (ib * sAGP)/nc2
                pEarly.ib = append(pEarly.ib, pbeta(fh.ebt, S+1, A+1) - pbeta(fl.ebt, S+1, A+1))
            }
        }
        
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pEarly.b + pLate + pEuploid

        if (length(pEarly.ib) > 0) {         #some states in between
            Ptot = Ptot + sum(pEarly.ib)
        }
        
        cp.A <- pEarly.a/Ptot
        cp.B <- pEarly.b/Ptot
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        cp.I <- 0
        if (length(pEarly.ib) > 0) {
           cp.I <- sum(pEarly.ib)/Ptot
        }
        #message(paste(c(cp.A, cp.B, cp.C, cp.D), collapse=" "))
        cp.AB <- 1 - cp.C - cp.D - cp.I
        cp.AC <- 1 - cp.B - cp.D - cp.I
        cp.AD <- 1 - cp.B - cp.D - cp.I
        cp.BC <- 1 - cp.A - cp.D - cp.I
        cp.BD <- 1 - cp.A - cp.C - cp.I
        cp.CD <- 1 - cp.A - cp.B - cp.I
        cp.ABC <- 1 - cp.D - cp.I
        cp.ABD <- 1 - cp.C - cp.I
        cp.ACD <- 1 - cp.B - cp.I
        cp.BCD <- 1 - cp.A - cp.I
        if (Ptot > 0) {
            if (cp.A >= sigTh) {                   # earl A
                evoType = "A1"
            } else if (cp.B >= sigTh){
                evoType <- "A2"
            } else if (cp.C >= sigTh){
                evoType <- "B"
            } else if (cp.D >= sigTh){
                evoType <- "C"
            } else if (cp.I >= sigTh){
                evoType <- "I"
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
        }
        allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid, sum(pEarly.ib))
        names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid", "pEarly.ib")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
            ccf = (f/pu)*nc - (nt - nb - 1)*pa          #early A1
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            #ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc2
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            #ccf2 = ccf2/pu
            sd <- Ms.A$SD
        #} else if (maxType == "pEarly.b" & prior != "late") {
        } else if (evoType == "I") {                   # early intermediate
            ibm = ((nb+1):(nt-nb-1))[which.max(pEarly.ib)]
            ccf = (f/pu)*nc - (ibm - 1)*pa             #early IB
            ff.I <- pu*(cc - pa +  ibm * pa)/nc        #dbinom
            Ms.I <- computeSD(N, S, ff.I)              #dbinom
            ccf2 <- Ms.I$M1                            #dbinom
            sd <- Ms.I$SD
        } else if (evoType == "A2" & prior != "late") {    
            ccf = (f/pu)*nc - (nb - 1)* pa         #early A2
            ff.B <- pu*(cc - pa + nb * pa)/nc      #dbinom
            Ms.B <- computeSD(N, S, ff.B)          #dbinom
            ccf2 <- Ms.B$M1                        #dbinom
            sd <- Ms.B$SD
        } else {
            ccf = (f/pu)*nc                        #other
            ff.C <- pu*cc/nc                       #dbinom
            Ms.C <- computeSD(N, S, ff.C)          #dbinom
            ccf2 <- Ms.C$M1                        #dbinom
            sd <- Ms.C$SD
        }
    }
    if ( f > 0.1 & ccf >= overadj ) {    #correct for over-adjustment
        if (evoType != "A1") {
            if ( (nt-nb) >= 3 ) {
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
    if ( f > 0.7 & ccf2 < 0.05 & is.nan(sd) ) {  #probably un-identified LOH
        ccf = f*2
        ccf2 = 1
        sd = 0.01
    }
    if (is.nan(sd)) {
        sd = 0.01
    }
    return(c(ccf, evoType, ccf2, sd))
}



computeCCF2 <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90) {  #after

    ccf = 0
    ccf2 = 0
    sd = 0
    cc <- seq(0.02, 1, by = 0.01)
    evoType = "A1/A2/B/C"
    N = A + S
    #nc = nt * pa + 2 * (1 - pa)
    nc2 = nt * sAGP + 2 * (1 - sAGP)                 #now adjust after getting mixed ccf
    nc = nc2

    message(paste("nc =", nc, sep=" "))
    #message(paste(c(f,A,S,pu,pa,sAGP,nt,nb),collapse=" "))
    if (nb == 1 & nt == 2) {     #normal diploid
        ccf = 2*f
        ff = cc/2
        Ms = computeSD(N, S, ff)
        ccf2 <- Ms$M1
        sd <- Ms$SD
    }
    else if (nt == 1) {          #het deletion
        ccf = f/nc
        ff.C <- cc/nc
        Ms.C <- computeSD(N, S, ff.C)          #dbinom
        ccf2 <- Ms.C$M1                        #dbinom
        sd <- Ms.C$SD
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
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
    } else if (nb == 0 & nt == 0) {        #homozygous deletion
        evoType <- "C"
        ccf = f*nc
        ff.C = cc[cc<=(1-sAGP)]/nc
        Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=(1-sAGP)])
        ccf2 <- (Ms.C$M1)/pu                    #dbinom
        sd <- Ms.C$SD
        #message(paste(c(ccf,ccf2),collapse=" "))
    } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cpEarly.a <- pEarly.a/Ptot
        cpLate.eup <- 1 - cpEarly.a
        cpLate <- pLate/Ptot
        cpEup <- pEuploid/Ptot
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
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
            ccf = f*nc - (nt - nb - 1)*sAGP
            #ff.early = (cc[cc>=p] - p + (nt - nb) * p)/nc
            ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc    #dbinom
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            sd <- Ms.A$SD
        } else {
            ccf = f*nc
            #ff.late = cc/nc
            ff.C <- cc/nc                        #dbinom
            Ms.C <- computeSD(N, S, ff.C)           #dbinom
            ccf2 <- Ms.C$M1                         #dbinom
            sd <- Ms.C$SD
        }
    } else if (nb >= 1 & nt > 2) {
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.eb <- (nb * sAGP + 1 - sAGP)/nc2
        fl.eb <- nb * sAGP/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
        cp.A <- pEarly.a/Ptot
        cp.B <- pEarly.b/Ptot
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        #message(paste(c(cp.A, cp.B, cp.C, cp.D), collapse=" "))
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
        if (Ptot > 0) {
            if (cp.A >= sigTh) {                   # earl A
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
        }
        allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid)
        names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
            ccf = f*nc - (nt - nb - 1)*sAGP         #early A1
            #ff.A <- (cc[cc>=p] - p + (nt - nb) * p)/nc
            ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc    #dbinom
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            #ccf2 = ccf2/pu
            sd <- Ms.A$SD
        #} else if (maxType == "pEarly.b" & prior != "late") {
        } else if (evoType == "A2" & prior != "late") {    
            ccf = f*nc - (nb - 1)* sAGP         #early A2
            #ff.B <- (cc[cc>=p] - p + nb * p)/nc
            ff.B <- (cc - sAGP + nb * sAGP)/nc     #dbinom
            Ms.B <- computeSD(N, S, ff.B)          #dbinom
            ccf2 <- Ms.B$M1                        #dbinom
            sd <- Ms.B$SD
        } else {
            ccf = f*nc                        #other
            ff.C <- cc/nc                       #dbinom
            Ms.C <- computeSD(N, S, ff.C)          #dbinom
            ccf2 <- Ms.C$M1                        #dbinom
            sd <- Ms.C$SD
        }
    }
    ccf = ccf/sAGP
    ccf2 = ccf/sAGP
    if ( f > 0.1 & ccf >= overadj ) {    #correct for over-adjustment
        if (evoType != "A1") {
            if ( (nt-nb) >= 3 ) {
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
    if (is.nan(sd)) {
        sd = 0.02
    }
    return(c(ccf, evoType, ccf2, sd))
}


computeSD <- function(N, S, f, cc=seq(0.02, 1, by = 0.01)) {
    M1list <- c()
    M2list <- c()
    MLElist <- c()
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
    bin = seq(minAF,1,0.02)                     # the bins I used
    mhist = hist(muts, breaks = bin, plot=F)
    mfreq = mhist$counts/sum(mhist$counts)

    for(k in 1:length(mfreq)) {
        if(mfreq[k]==0) {
            mfreq[k] = 0.00000001
        }
    }
    return(mfreq)
}


JS.divergence <- function(sub1,sub2, minAF=0.04) {
    sub1 = sub1[which(sub1 <= 1 & sub1 >= minAF)]
    sub2 = sub2[which(sub2 <= 1 & sub2 >= minAF)]
    freqs1 = pmf(sub1, minAF)
    freqs2 = pmf(sub2, minAF)
    m<-0.5*(freqs1 + freqs2)
    JS<-0.5*(sum(freqs1*log2(freqs1/m)) + sum(freqs2*log2(freqs2/m)))
    return(JS)
}



pubOrSub <- function(sampAB, samples, minAF=0.05, minDepTotal=5*length(samples), groupName = "", pAF=0.25, setpub=FALSE) {

    originalColNames = colnames(sampAB)

    maxPa = vector()
    maxsAGP = vector()
    for (i in 1:length(samples)) {           #each sample
        sn = samples[i]
        pai = match(paste(sn, "pa", sep=""), originalColNames)
        sAGPi = match(paste(sn, "sAGP", sep=""), originalColNames)
        
        maxPa = c(maxPa, max(sampAB[,pai]))
        maxsAGP = c(maxsAGP, max(sampAB[,sAGPi]))
    }
    message(paste(maxPa,collapse="\t"))
    message(paste(maxsAGP,collapse="\t"))

    pstype = as.vector(apply(sampAB, 1, function(x, coln, samples, maxPa, maxsAGP, minAF, minDepTotal) {
                                 mutVector = as.vector(x)
                                 cppres = pubOrSub.Calc(x, coln, samples, maxPa, maxsAGP, pAF=pAF, setpub=setpub)

                                 cpstype = "unknown"
                                 totalDepth = sum(as.numeric(mutVector[match(paste(samples,"d",sep=""),coln)]))
                                 if (is.nan(cppres[1]) | is.nan(cppres[3])) {
                                     message(paste(x,collapse="\t"))
                                     message(paste(cppres, collapse="\t"))
                                 }
                                 if (totalDepth < minDepTotal) {
                                     cpstype = "unknown"
                                 } else if (cppres[1] >= 0.05 | cppres[3] == 1) {   #accept public
                                     cpstype = "public"
                                 } else if (cppres[2] >= 0.05) {                    #accept absent
                                     cpstype = "absent"
                                 } else {                                           #subclonal
                                     foundSamples = ""
                                     for (i in 1:length(samples)) {                 #get subclonal info for each sample
                                         sn = samples[i]
                                         mafai = match(paste(sn, "mafa", sep=""), coln)
                                         nbi = match(paste(sn, "nb", sep=""), coln)
                                         depthi = match(paste(sn, "d", sep=""), coln)
                                         if (as.numeric(mutVector[mafai]) > minAF) {
                                             foundSamples = paste(foundSamples, sn, ",", sep="")
                                         }
                                     }
                                     if (foundSamples != ""){
                                         foundSamples = gsub(",$","", foundSamples)
                                         cpstype = paste("private", foundSamples, sep="=")
                                     }
                                 }
                                 cpstype
                             }, coln=originalColNames, samples=samples, maxPa = maxPa, maxsAGP = maxsAGP, minAF=minAF, minDepTotal=minDepTotal))

    if (groupName == "") {
        if ("pubOrSub" %in% originalColNames) {
            sampAB$pubOrSub = pstype
        } else {
            sampAB = data.frame(sampAB, pubOrSub=pstype)
            colnames(sampAB) = c(originalColNames,"pubOrSub")
        }
    } else {
        if (groupName %in% originalColNames) {
            sampAB[,groupName] = pstype
        } else {
            sampAB = data.frame(sampAB, pubOrSubTmp=pstype)
            colnames(sampAB) = c(originalColNames,groupName)
        }
    }
    return(sampAB)

}


pubOrSub.Calc <- function(mutVector, originalNames, samples, maxPa = vector(), maxsAGP = vector(), pAF=0.25, setpub=FALSE) {

    CCFaboveOne = all((as.numeric(mutVector[match(paste(samples, "ccf", sep=""), originalNames)]) +
                           2.58*as.numeric(mutVector[match(paste(samples, "ccfSD", sep=""), originalNames)])) >= 1)
    aboveContri = length(which((as.numeric(mutVector[match(paste(samples, "ccf", sep=""), originalNames)]) +
                                    2.58*as.numeric(mutVector[match(paste(samples, "ccfSD", sep=""), originalNames)])) >= 1))
    VAFaboveQua = all(as.numeric(mutVector[match(paste(samples, "mafa", sep=""), originalNames)]) >= pAF)

    cpop = sapply(samples, function(x, samples, originalNames, mutVector, aboveContri, maxPa, maxsAGP) {
                      sn = x
                      si = match(x, samples)
                      ccfi = match(paste(sn, "ccf", sep=""), originalNames)
                      ccfsdi = match(paste(sn, "ccfSD", sep=""), originalNames)
                      mafai = match(paste(sn, "mafa", sep=""), originalNames)
                      refci = match(paste(sn, "refc", sep=""), originalNames)
                      altci = match(paste(sn, "altc", sep=""), originalNames)
                      pui = match(paste(sn, "pu", sep=""), originalNames)
                      pai = match(paste(sn, "pa", sep=""), originalNames)
                      sAGPi = match(paste(sn, "sAGP", sep=""), originalNames)
                      nti = match(paste(sn, "nt", sep=""), originalNames)
                      nbi = match(paste(sn, "nb", sep=""), originalNames)
                      depthi = match(paste(sn, "d", sep=""), originalNames)

                      cPa = as.numeric(mutVector[pai])
                      prior = "unknown"
                      if (cPa == maxPa[si] | cPa == 0) {
                          cPa = 1
                      }
                      if ((aboveContri >= 2 | (aboveContri >= 1 & length(samples) == 2)) & setpub == TRUE) {                    #looks like public, change prior
                          prior = "pub"
                      }
                      cpop.probs = pubOrSub.prob.binom(as.integer(mutVector[refci]), as.integer(mutVector[altci]), as.numeric(maxsAGP[si]), as.numeric(mutVector[sAGPi]), as.integer(mutVector[nti]), as.integer(mutVector[nbi]), pa=cPa, prior=prior)
                      cpop.probs
                  }, originalNames=originalNames, samples=samples, mutVector=mutVector, aboveContri=aboveContri, maxPa=maxPa, maxsAGP=maxsAGP)

    cpp = prod(cpop[1,])
    cpa = prod(cpop[2,])
    if (is.nan(cpp)){
        message(paste(cpop[1,],collapse="\t"))
    }
    
    return(c(cpp, cpa, as.numeric(CCFaboveOne | VAFaboveQua)))

}

pubOrSub.simu <- function(sampAB, samples, minAF=0.05, minDepTotal=5*length(samples), groupName = "", pAF=0.25) {    #for simulated data

    originalColNames = colnames(sampAB)

    pstype = as.vector(apply(sampAB, 1, function(x, coln, samples, pAF, minAF, minDepTotal) {
                                 mutVector = x
                                 names(mutVector) = coln
                                 cppres = pubOrSub.Calc.sim(mutVector, samples, rep(1,length(samples)), rep(1,length(samples)), pAF=pAF)
                                 cpstype = "unknown"
                                 totalDepth = sum(as.numeric(mutVector[paste("depth",samples,sep="")]))
                                 if (totalDepth < minDepTotal) {
                                     cpstype = "unknown"
                                 } else if (cppres[1] >= 0.05 | cppres[3] == 1) {   #accept public
                                     cpstype = "public"
                                 } else if (cppres[2] >= 0.05) {  #accept absent
                                     cpstype = "absent"
                                 } else {     #subclonal
                                     foundSamples = ""
                                     for (i in 1:length(samples)) {                                     #get subclonal ones 2nd round: Refine
                                         sn = samples[i]
                                         mafai = match(paste("maf", sn, sep=""), coln)
                                         depthi = match(paste("depth", sn, sep=""), coln)
                                         if (as.numeric(mutVector[mafai]) > minAF) {
                                             foundSamples = paste(foundSamples, sn, ",", sep="")
                                         }
                                     }
                                     if (foundSamples != "") {
                                         foundSamples = gsub(",$","", foundSamples)
                                         cpstype = paste("private", foundSamples, sep="=")
                                     }
                                 }
                                 cpstype
                             }, coln=originalColNames, samples=samples, pAF=pAF, minAF=minAF, minDepTotal=minDepTotal))
    
    if (groupName == "") {
        if ("pubOrSub" %in% originalColNames) {
            sampAB$pubOrSub = pstype
        } else {
            sampAB = data.frame(sampAB, pubOrSub=pstype)
            colnames(sampAB) = c(originalColNames,"pubOrSub")
        }
    } else {
        if (groupName %in% originalColNames) {
            sampAB[,groupName] = pstype
        } else {
            sampAB = data.frame(sampAB, pubOrSubTmp=pstype)
            colnames(sampAB) = c(originalColNames,groupName)
        }
    }
    return(sampAB)
}


pubOrSub.Calc.sim <- function(mutVector, samples, maxPa = vector(), maxsAGP = vector(), pAF=0.25) {

    originalNames = names(mutVector)

    aboveContri = 0
    VAFaboveQua = all(as.numeric(mutVector[match(paste("maf", samples, sep=""), originalNames)]) >= pAF)

    cpop = sapply(samples, function(x, originalNames, mutVector) {
                      sn = x
                      mafai = match(paste("maf", sn, sep=""), originalNames)
                      depthi = match(paste("depth", sn, sep=""), originalNames)
                      depth = as.numeric(mutVector[depthi])
                      maf = as.numeric(mutVector[mafai])
                      refc = round(depth * (1-maf))
                      altc = round(depth * maf)

                      cPa = 1
                      prior = "unknown"
                      
                      cpop.probs = pubOrSub.prob.binom(refc, altc, 1, 1, 2, 1, pa=1, prior=prior)
                      cpop.probs
                  }, originalNames=originalNames, mutVector=mutVector)
    
    cpp = prod(cpop[1,])
    cpa = prod(cpop[2,])

    return(c(cpp, cpa, as.numeric(VAFaboveQua)))

}


pubOrSub.prob <- function(A, S, pu, sAGP, nt, nb, pa=1, prior="unknown") {   #based on beta, bayesian model
    
    N = A + S
    
    nc = nt * sAGP + 2 * (1 - sAGP)
    #fh.pub <- (sAGP * (nt - nb) + 1 - sAGP)/nc
    fh.pub <- 1
    fl.sub <- 0.02
    fl.pub <- (pu * (nt - nb))/nc
    fh.sub <- pu/nc
    
    #set for special cases
    if (nb == 0 & nt == 0) {        #homo deletion, wierd case
        fh.pub <- pu
        fl.pub <- 0.5*pu
        fh.sub <- 0.5*pu
    } else if (nb == 0 | nt == 2 * nb) {
        if (prior == "pub"){
            fl.pub <- (pu - sAGP)/nc
        }
    } else if (nb >= 1 & nt > 2) {
        fl.pub <- nb * pu/nc
    }

    #pp = pbinom(S,N,prob)
    p.pub = pbeta(fh.pub,S+1,A+1) - pbeta(fl.pub,S+1,A+1)
    p.sub = pbeta(fh.sub,S+1,A+1) - pbeta(fl.sub,S+1,A+1)
    p.abs = pbeta(fl.sub,S+1,A+1)
    return(c(p.pub, p.sub, p.abs))
    
}


pubOrSub.prob.binom <- function(A, S, pu, sAGP, nt, nb, pa=1, prior="unknown") {      #using binomial test
    
    N = A + S
    
    nc = nt * sAGP + 2 * (1 - sAGP)
    f.abs <- 0.02
    f.pub <- (pu * (nt - nb))/nc

    f.pub = ifelse(nb == 0 & nt == 0, 0.5*pu,
            ifelse(nb == 0 & prior == "pub", max((pu - sAGP)/nc, 0),
            ifelse(nt >= 2, pu/nc, (pu * (nt - nb))/nc)))

    f.pub = min(1,f.pub)
    #p.pub = pbeta(fh.pub,S+1,A+1) - pbeta(fl.pub,S+1,A+1)
    p.pub = pbinom(S, N, f.pub)       #p for pub (under pub vaf, the prob reaching S)
    #p.abs = 1 - pbinom(S, N, f.abs)   #p for abs (under abs, the prob reaching S)
    p.abs = pbeta(f.abs, S+1, A+1)
    return(c(p.pub, p.abs))
    
}


earlyMut.prob.binom <- function(A, S, pu, sAGP, nt, nb, pa=1, prior="unknown") {      #using binomial test
    
    N = A + S
    
    nc = nt * sAGP + 2 * (1 - sAGP)
    f.abs <- 0.02
    f.pub <- (pu * (nt - nb))/nc

    f.pub = ifelse(nb == 0 & nt == 0, 0.5*pu,
            ifelse(nb == 0 & prior == "pub", max((pu - sAGP)/nc, 0),
            ifelse(nt >= 2, (pu * (nt - nb))/nc, (pu * (nt - nb))/nc)))     #as long as more than 2 copies

    f.pub = min(1,f.pub)
    p.pub = pbinom(S, N, f.pub)       #p for pub (under pub vaf, the prob reaching S)
    #p.abs = 1 - pbinom(S, N, f.abs)   #p for abs (under abs, the prob reaching S)
    p.abs = pbeta(f.abs, S+1, A+1)
    return(c(p.pub, p.abs))
    
}



fst.wc84 <- function(af, minAF=0.08) {   #wc84 is sentitive to sample size
       mafis = which(grepl("maf", colnames(af)))
       keep = as.vector(apply(af, 1, function(x, mafis) {
                                  maxmaf = max(as.numeric(x[mafis]))
                                  if (maxmaf > minAF) {
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


fst.hudson <- function(af, minAF=0.08) {  #combine for multi SSNV sites using ratio of averages
    mafis = which(grepl("maf", colnames(af)))
    keep = as.vector(apply(af, 1, function(x, mafis) {
                               maxmaf = max(as.numeric(x[mafis]))
                               if (maxmaf > minAF) {
                                   TRUE
                               } else {
                                   FALSE
                               }
                           }, mafis=mafis))     #filter data
    af = af[keep,]
    
    Ns = c()
    Ds = c()
    for(k in 1:nrow(af)) {
        n1 = af$depth1[k]
        n2 = af$depth2[k]
        p1 = af$maf1[k]
        p2 = af$maf2[k]
        N = (p1-p2)^2-(p1*(1-p1))/(n1-1)-(p2*(1-p2))/(n2-1)     # covariance
        D = p1*(1-p2)+p2*(1-p1)                                 # standard deviations
        Ns = c(Ns, N)
        Ds = c(Ds, D)
    }
    Fst.h = mean(Ns)/mean(Ds)
    return(Fst.h)
}


dNdS <- function(sampAB, g.dnds, ndriver = 220, talt = 4, pscn="pubOrSub") {
    
    # pub dNdS
    cateTF = sampAB[,pscn] == "public" & sampAB$totalAlt >= talt
    res.pub = dNdS.cate(sampAB, cateTF, g.dnds, ndriver)
    
    # sub dNdS
    cateTF = grepl("private",sampAB[,pscn]) & sampAB$totalAlt >= talt
    res.sub = dNdS.cate(sampAB, cateTF, g.dnds, ndriver)
    
    res = c(res.pub, res.sub)
    names(res) = c("pub.dnds", "pub.nonsyn", "pub.syn", "pub.norm", "pub.dmfdlf.cadd", "pub.mf.cadd", "pub.lf.cadd", "pub.cadd.norm","pub.dmfdlf.gerp", "pub.mf.gerp", "pub.lf.gerp", "pub.gerp.norm",
             "pub.dmfdlf.sift", "pub.mf.sift", "pub.lf.sift", "pub.sift.norm","pub.dmfdlf.polyphen", "pub.mf.polyphen", "pub.lf.polyphen", "pub.polyphen.norm","npub", "npubDriver",  "pub.driver.enrich.p", "pub.driver.fc",
             "sub.dnds", "sub.nonsyn", "sub.syn", "sub.norm", "sub.dmfdlf.cadd", "sub.mf.cadd", "sub.lf.cadd", "sub.cadd.norm","sub.dmfdlf.gerp", "sub.mf.gerp", "sub.lf.gerp", "sub.gerp.norm",
             "sub.dmfdlf.sift", "sub.mf.sift", "sub.lf.sift", "sub.sift.norm","sub.dmfdlf.polyphen", "sub.mf.polyphen", "sub.lf.polyphen", "sub.polyphen.norm", "nsub", "nsubDriver",  "sub.driver.enrich.p","sub.driver.fc")
    return(res)
}

dNdS2 <- function(sampAB, g.dnds, ndriver = 220, talt = 4, pscn=c("pubOrSub_P","pubOrSub_BM")) {

    #############################################################################
    cateTF = sampAB[,pscn[2]] == "public" & sampAB[,pscn[1]] == "public" & sampAB$totalAlt >= talt
    # pub dNdS NS
    res.pub = dNdS.cate(sampAB, cateTF, g.dnds, ndriver)

    #############################################################################
    cateTF = sampAB[,pscn[2]] == "public" & grepl("private",sampAB[,pscn[1]]) & sampAB$totalAlt >= talt
    # sub dNdS
    res.sub = dNdS.cate(sampAB, cateTF, g.dnds, ndriver)

    #############################################################################
    cateTF = sampAB[,pscn[2]] == "public" & sampAB[,pscn[1]] == "absent" & sampAB$totalAlt >= talt
    # specific dNdS
    res.spec = dNdS.cate(sampAB, cateTF, g.dnds, ndriver)
    
    res = c(res.pub, res.sub, res.spec)
    names(res) = c("pub.dnds", "pub.nonsyn", "pub.syn", "pub.norm", "pub.dmfdlf.cadd", "pub.mf.cadd", "pub.lf.cadd", "pub.cadd.norm","pub.dmfdlf.gerp", "pub.mf.gerp", "pub.lf.gerp", "pub.gerp.norm",
             "pub.dmfdlf.sift", "pub.mf.sift", "pub.lf.sift", "pub.sift.norm","pub.dmfdlf.polyphen", "pub.mf.polyphen", "pub.lf.polyphen", "pub.polyphen.norm","npub", "npubDriver",  "pub.driver.enrich.p", "pub.driver.fc",
             "sub.dnds", "sub.nonsyn", "sub.syn", "sub.norm", "sub.dmfdlf.cadd", "sub.mf.cadd", "sub.lf.cadd", "sub.cadd.norm","sub.dmfdlf.gerp", "sub.mf.gerp", "sub.lf.gerp", "sub.gerp.norm",
             "sub.dmfdlf.sift", "sub.mf.sift", "sub.lf.sift", "sub.sift.norm","sub.dmfdlf.polyphen", "sub.mf.polyphen", "sub.lf.polyphen", "sub.polyphen.norm", "nsub", "nsubDriver",  "sub.driver.enrich.p","sub.driver.fc", 
             "spec.dnds", "spec.nonsyn", "spec.syn", "spec.norm", "spec.dmfdlf.cadd", "spec.mf.cadd", "spec.lf.cadd", "spec.cadd.norm","spec.dmfdlf.gerp", "spec.mf.gerp", "spec.lf.gerp", "spec.gerp.norm",
             "spec.dmfdlf.sift", "spec.mf.sift", "spec.lf.sift", "spec.sift.norm","spec.dmfdlf.polyphen", "spec.mf.polyphen", "spec.lf.polyphen", "spec.polyphen.norm","nspec", "nspecDriver","spec.driver.enrich.p", "spec.driver.fc")
    return(res)
}

dNdS.cate <- function(sampAB, cateTF, g.dnds, ndriver) {

    g.dnds2 = g.dnds[which(g.dnds$MFLFsift < 100),]
    
    # dNdS NS
    g.table1 = table(sampAB$geneName[which(cateTF & (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])
    g.norm1 = dNdS.norm(g.table1, "NS", g.dnds)
    f.table1 = table(sampAB$functionalClass[which(cateTF & (grepl("synonymous",sampAB$functionalClass) | grepl("stop",sampAB$functionalClass)))])
    nonsyn1 = as.numeric(sum(as.numeric(f.table1[which(names(f.table1) != "synonymous SNV" & names(f.table1) != "unknown")])))
    syn1 = as.numeric(f.table1["synonymous SNV"])
    res.dnds1 = nonsyn1/(syn1+0.1)
    if (g.norm1 != 0){
        res.dnds1 = res.dnds1/g.norm1
    }
    
    # cadd MFLFcadd
    cadd.g.table1 = table(sampAB$geneName[which(cateTF & !is.na(sampAB$CADD_phred))])
    cadd.g.norm1 = dNdS.norm(cadd.g.table1, "MFLFcadd", g.dnds)
    cadd.f.more1 = length(which(cateTF & !is.na(sampAB$CADD_phred) & sampAB$CADD_phred >= 20))
    cadd.f.less1 = length(which(cateTF & !is.na(sampAB$CADD_phred) & sampAB$CADD_phred < 20))
    res.cadd1 = cadd.f.more1/(cadd.f.less1+0.1)
    if (cadd.g.norm1 != 0){
        res.cadd1 = res.cadd1/cadd.g.norm1
    }

    # GERP MFLFgerp
    gerp.g.table1 = table(sampAB$geneName[which(cateTF & !is.na(sampAB$GERP_RS))])
    gerp.g.norm1 = dNdS.norm(gerp.g.table1, "MFLFgerp", g.dnds)
    gerp.f.more1 = length(which(cateTF & !is.na(sampAB$GERP_RS) & sampAB$GERP_RS >= 5.2))
    gerp.f.less1 = length(which(cateTF & !is.na(sampAB$GERP_RS) & sampAB$GERP_RS < 5.2))
    res.gerp1 = gerp.f.more1/(gerp.f.less1+0.1)
    if (gerp.g.norm1 != 0){
        res.gerp1 = res.gerp1/gerp.g.norm1
    }

    # SIFT MFLFsift
    sift.g.table1 = table(sampAB$geneName[which(cateTF & !is.na(sampAB$SIFT_score))])
    sift.g.norm1 = dNdS.norm(sift.g.table1, "MFLFsift", g.dnds2)
    sift.f.more1 = length(which(cateTF & !is.na(sampAB$SIFT_score) & sampAB$SIFT_score <= 0.05))
    sift.f.less1 = length(which(cateTF & !is.na(sampAB$SIFT_score) & sampAB$SIFT_score > 0.05))
    res.sift1 = sift.f.more1/(sift.f.less1+0.1)
    if (sift.g.norm1 != 0) {
        res.sift1 = res.sift1/sift.g.norm1
    }

    # POLYPHEN MFLFpolyphen
    polyphen.g.table1 = table(sampAB$geneName[which(cateTF & !is.na(sampAB$Polyphen2_HVAR_pred))])
    polyphen.g.norm1 = dNdS.norm(polyphen.g.table1, "MFLFpolyphen", g.dnds)
    polyphen.f.more1 = length(which(cateTF & !is.na(sampAB$Polyphen2_HVAR_pred) & sampAB$Polyphen2_HVAR_pred != "B"))
    polyphen.f.less1 = length(which(cateTF & !is.na(sampAB$Polyphen2_HVAR_pred) & sampAB$Polyphen2_HVAR_pred == "B"))
    res.polyphen1 = polyphen.f.more1/(polyphen.f.less1+0.1)
    if (polyphen.g.norm1 != 0){
        res.polyphen1 = res.polyphen1/polyphen.g.norm1
    }

    # driver
    nmut = length(which(cateTF & sampAB$geneLoc != "intergenic" & sampAB$functionalClass != "synonymous SNV"))
    nmutDriver = length(which(cateTF & sampAB$geneLoc != "intergenic" & sampAB$functionalClass != "synonymous SNV" & sampAB$dron != 0))
    mut.driver.enrich.p = round(phyper(nmutDriver, ndriver, 20000-ndriver, nmut, lower.tail=F),8)
    mut.driver.fc = round((nmutDriver/(nmut+0.1))/(ndriver/20000),8)

    return(c(res.dnds1, nonsyn1, syn1, g.norm1, res.cadd1, cadd.f.more1, cadd.f.less1, cadd.g.norm1, res.gerp1, gerp.f.more1, gerp.f.less1, gerp.g.norm1,
             res.sift1, sift.f.more1, sift.f.less1, sift.g.norm1, res.polyphen1, polyphen.f.more1, polyphen.f.less1, polyphen.g.norm1,
             nmut, nmutDriver, mut.driver.enrich.p, mut.driver.fc))

}


dNdS.norm <- function(g.table, coln, g.dnds) {
    g.table = g.table[which(g.table > 0)]
    g.match = names(g.table) %in% g.dnds$gene
    g.table = g.table[g.match]
    subset.dnds = g.dnds[match(names(g.table),g.dnds$gene), match(coln,colnames(g.dnds))]
    g.total = sum(g.table)
    g.frac = g.table/g.total
    g.norm = sum(g.frac*subset.dnds)
    return(g.norm)
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



scicloneTable <- function(sampAB, samples, outdir, minAF = 0.05, maxAF = 0.8, minDepth = 30, titanPath="./titan/", log=FALSE) {

    colnames = colnames(sampAB)
    keep = as.vector(apply(sampAB, 1, function(x, coln, minAF, maxAF, minDepth) {
                               maxSampAF = max(as.numeric(x[match(paste(samples, "mafa", sep=""), coln)]))
                               minSampAF = min(as.numeric(x[match(paste(samples, "mafa", sep=""), coln)]))
                               minSampDP = min(as.numeric(x[match(paste(samples, "d", sep=""), coln)]))      
                               if (maxSampAF > minAF & minSampAF < maxAF & minSampDP >= minDepth) TRUE
                               else FALSE
                           }, coln=colnames, minAF=minAF, maxAF=maxAF, minDepth=minDepth))
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
        if (log){
            mafa = log2(mafa + 1) * 10
        }
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


mergeVAF <- function(sampAB, samples) {
    rindexTF = vector()
    cimafa = match(paste(samples,"mafa",sep=""), colnames(sampAB))
    cid = match(paste(samples,"d",sep=""), colnames(sampAB))
    for (si in 1:length(samples)) {
        if (length(rindexTF) > 0) {
            rindexTF = rindexTF | grepl(paste(samples[si],"\\[",sep=""), sampAB$somatic)
        } else {
            rindexTF = grepl(paste(samples[si],"\\[",sep=""), sampAB$somatic)
        }
    }
    rindex = which(rindexTF)
    mergeRes = data.frame(chr = sampAB[rindex,"chr"], pos = sampAB[rindex,"pos"], ref_reads = 0, var_reads=0, vaf = 0)
    for (ri in 1:length(rindex)) {
        cd = sampAB[rindex[ri],cid]
        cmafa = sampAB[rindex[ri],cimafa]
        mergedMAFA = sum(cmafa * cd)/sum(cd)
        mergedDEP = sum(cd)
        mergedALT = round(mergedDEP * mergedMAFA)
        mergedREF = mergedDEP - mergedALT
        mergeRes$ref_reads[ri] = mergedREF
        mergeRes$var_reads[ri] = mergedALT
        mergeRes$vaf[ri] = round(mergedMAFA * 100,2)
    }
    return(mergeRes)
}



prepareLichee <- function(samples, nmaf, SampAB, minDepth=20, clonal=FALSE, correctColname=FALSE, toPhylip=FALSE) {
    depthCols = paste(samples, "d", sep="")
    ccfCols = paste(samples, "ccf", sep="")
    ccfSDcols = paste(samples, "ccfSD", sep="")
    licheeCol = c("chr","pos","ref","alt","geneName",nmaf,paste(samples,"mafc",sep=""))
    #licheeRow = !is.na(SampAB$functionalClass)
    licheeRow = rep(TRUE, dim(SampAB)[1])
    for (i in 1:length(depthCols)) {  #Depth
        licheeRow = licheeRow & SampAB[,match(depthCols[i], colnames(SampAB))] >= minDepth
    }
    #message(length(which(licheeRow)))
    #message(ccfCols)
    
    if (clonal) {
        licheeRow2 = rep(FALSE, dim(SampAB)[1])
        for (i in 1:length(ccfCols)) {  #clonality
            licheeRow2 = licheeRow2 | (SampAB[,match(ccfCols[i], colnames(SampAB))] + 1.96*SampAB[,match(ccfSDcols[i], colnames(SampAB))]) >= 1
        }
        licheeRow = licheeRow & licheeRow2
    }
    licheeInput = SampAB[licheeRow,match(licheeCol,colnames(SampAB))]
    
    licheeInput = data.frame(licheeInput, name=paste(licheeInput$geneName,
                                              paste(licheeInput$chr,licheeInput$pos,licheeInput$ref,sep=":"), licheeInput$alt, sep="_"))
    if ( correctColname == TRUE ) {
        colnames(licheeInput) = gsub("\\.","-", colnames(licheeInput))
    }
    #message(paste(colnames(licheeInput), collapse="\t"))
    #message(paste(dim(licheeInput), collapse="\t"))
    licheeInput = licheeInput[,match(c("chr","pos","name",nmaf,paste(samples,"mafc",sep="")),colnames(licheeInput))]
    licheeInput[,4] = 0
    
    colnames(licheeInput) = gsub("mafc", "", colnames(licheeInput))
    colnames(licheeInput) = gsub("maf", "", colnames(licheeInput))
    return(licheeInput)

}

prepareTreeomics <- function(samples, nmaf, nd, SampAB, minDepth=20, clonal=FALSE, correctColname=FALSE) {
    
    depthCols = paste(samples, "d", sep="")
    ccfCols = paste(samples, "ccf", sep="")
    ccfSDcols = paste(samples, "ccfSD", sep="")
    licheeCol = c("chr","pos","ref","alt","geneName","geneLoc","functionalClass",nmaf,paste(samples,"altc",sep=""),depthCols,nd)

    licheeRow = rep(TRUE, dim(SampAB)[1])
    for (i in 1:length(depthCols)) {  #Depth
        licheeRow = licheeRow & SampAB[,match(depthCols[i], colnames(SampAB))] >= minDepth
    }
    
    if (clonal) {
        licheeRow2 = rep(FALSE, dim(SampAB)[1])
        for (i in 1:length(ccfCols)) {  #clonality
            licheeRow2 = licheeRow2 | (SampAB[,match(ccfCols[i], colnames(SampAB))] + 1.96*SampAB[,match(ccfSDcols[i], colnames(SampAB))]) >= 1
        }
        licheeRow = licheeRow & licheeRow2
    }
    licheeInput = SampAB[licheeRow,match(licheeCol,colnames(SampAB))]

    idxgeneName = match("geneName",colnames(licheeInput))
    idxgeneLoc = match("geneLoc",colnames(licheeInput))
    idxfunClass = match("functionalClass",colnames(licheeInput))
    GeneNames = as.vector(apply(licheeInput, 1, function(x, ign, igl, ifc){
                                    if (x[igl] == "intergenic" | x[igl] == "ncRNA_exonic" | x[igl] == "ncRNA_splicing" | x[igl] == "splicing") {
                                        x[ign]
                                    } else if (x[igl] == "exonic" & x[ifc] != "synonymous SNV") {
                                        x[ign]
                                    } else if (x[igl] == "exonic" & x[ifc] == "synonymous SNV") {
                                        paste(x[ign],"[synonymous]",sep="")
                                    } else {
                                        paste(x[ign],"[",x[igl],"]",sep="")
                                    }
                                                            }, ign=idxgeneName, igl=idxgeneLoc, ifc=idxfunClass))
    
    licheeInput = data.frame(licheeInput, Change=paste(licheeInput$ref, licheeInput$alt, sep=">"), Gene=GeneNames)
    if ( correctColname == TRUE ) {
        colnames(licheeInput) = gsub("\\.","-", colnames(licheeInput))
    }

    licheeInput1 = licheeInput[,match(c("chr","pos","Change","Gene",nmaf,paste(samples,"altc",sep="")),colnames(licheeInput))]
    licheeInput1[,5] = 0
    colnames(licheeInput1)[1:4] = c("Chromosome","Position","Change","Gene")
    colnames(licheeInput1) = gsub("altc", "", colnames(licheeInput1))
    colnames(licheeInput1) = gsub("maf", "", colnames(licheeInput1))
    colnames(licheeInput1) = gsub("SPCG-", "", colnames(licheeInput1))

    licheeInput2 = licheeInput[,match(c("chr","pos","Change","Gene",nd,paste(samples,"d",sep="")),colnames(licheeInput))]
    colnames(licheeInput2)[1:4] = c("Chromosome","Position","Change","Gene")
    colnames(licheeInput2) = gsub("d$", "", colnames(licheeInput2))
    colnames(licheeInput2) = gsub("SPCG-", "", colnames(licheeInput2))
    return(list(licheeInput1, licheeInput2))

}



prepareLichee2 <- function(samples, sgSamples, nmaf, sampAB) {
    
    licheeCol = c("chr","pos","ref","alt","geneName",nmaf,paste(c(samples, sgSamples),"mafa",sep=""))
    sn1 = samples[1]
    sn2 = samples[2]
    submuts = subclonalMut(sampAB, sn1, sn2, minAF=0.01)
    needr = vector()
    sgHighi = vector()

    depthSample = vector()
    for (i in 1:length(samples)) {  #Depth of original samples
        sn = samples[i]
        depi = match(paste(sn, "d", sep=""), colnames(sampAB))
        if (length(depthSample) > 0) {
            depthSample = (depthSample & sampAB[,depi] >= 20)
        } else {
            depthSample = (sampAB[,depi] >= 20)
        }
    }
    
    for (i in 1:length(sgSamples)) {
        sn = sgSamples[i]
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        depi = match(paste(sn, "d", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        nbSN1i = match(paste(sn1, "nb", sep=""), colnames(sampAB))
        nbSN2i = match(paste(sn2, "nb", sep=""), colnames(sampAB))
        if (length(needr) == 0) {
            sgHighi = which(as.numeric(sampAB[,mafai]) > 0.15)
            needr = which(as.numeric(sampAB[,depi]) >= 20 & (as.numeric(sampAB[,nbi]) != 0 | as.numeric(sampAB[,nbSN1i]) == 0 | as.numeric(sampAB[,nbSN2i]) == 0))
        } else {
            sgHighi = union(sgHighi, which(as.numeric(sampAB[,mafai]) > 0.15))
            needr = intersect(needr, which(as.numeric(sampAB[,depi]) >= 20 & (as.numeric(sampAB[,nbi]) != 0 | as.numeric(sampAB[,nbSN1i]) == 0 | as.numeric(sampAB[,nbSN2i]) == 0)) )
        }
    }
    
    pubTi = intersect(needr, submuts$pubTi)
    subAi = intersect(needr, submuts$subAi)
    subBi = intersect(needr, submuts$subBi)
    sgHighi = intersect(needr, sgHighi)
    licheeRow = union(sgHighi, union(pubTi, union(subAi, subBi)))
    licheeRow = intersect(licheeRow, which(depthSample))
    licheeRow = intersect(licheeRow, which(!grepl("chrM", sampAB$chr)))

    
    licheeInput = sampAB[licheeRow,match(licheeCol,colnames(sampAB))]
    licheeInput = data.frame(licheeInput, name=paste(licheeInput$geneName,
                  paste(licheeInput$chr,licheeInput$pos,licheeInput$ref,sep=":"), licheeInput$alt, sep="_"))
    licheeInput = licheeInput[,match(c("chr","pos","name",nmaf,paste(c(samples,sgSamples),"mafa",sep="")),colnames(licheeInput))]
    licheeInput[,4] = 0

    colnames(licheeInput) = gsub("mafa", "", colnames(licheeInput))
    colnames(licheeInput) = gsub("maf", "", colnames(licheeInput))
    colnames(licheeInput) = gsub("ccf", "", colnames(licheeInput))

    profile = as.vector(apply(licheeInput, 1, function(x, samps, sgSamples, colnames) {
                        profile = rep(0, length(c(samples,sgSamples)))
                        profile[which(as.numeric(x[match(samps, colnames)]) > 0.015)] = 1
                        profile[length(samps) + which(as.numeric(x[match(sgSamples, colnames)]) > 0.15)] = 1
                        profile = paste(profile, sep="", collapse="")
                        profile
                    }, colnames = colnames(licheeInput), samps = samples, sgSamples = sgSamples))
    licheeInput = data.frame(licheeInput, profile=paste("0",profile,sep=""))
    licheeInput = licheeInput[,match(c("chr","pos","name","profile",colnames(licheeInput)[4],c(samples,sgSamples)),colnames(licheeInput))]

    needProfile = names(table(licheeInput$profile)[which(table(licheeInput$profile) >= 4)])
    needProfile = needProfile[which(grepl("1", needProfile))]
    licheeInput = licheeInput[licheeInput$profile %in% needProfile,]

    clusters = vector()
    for (i in 1:length(needProfile)) {
        cpro = needProfile[i]
        rowsi = which(licheeInput$profile == cpro)
        rows = paste(rowsi, sep=",", collapse=",")
        cdata = licheeInput[rowsi, match(c(colnames(licheeInput)[5], samples,sgSamples), colnames(licheeInput))]
        colmeans = as.vector(apply(cdata, 2, median))
        clusters = rbind(clusters, c(cpro, colmeans, rows))
    }
    return(list(li=licheeInput, cl=clusters))
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


generateUpSetList2 <- function(sampAB, sn1, sn2, samples) {
    
    outlist = list()
    
    submuts = subclonalMut(sampAB, sn1, sn2, minAF=0.015)
    depthBulkSample = vector()
    for (i in 1:length(c(sn1,sn2))) {  #Depth of original samples
        sn = c(sn1,sn2)[i]
        depi = match(paste(sn, "d", sep=""), colnames(sampAB))
        if (length(depthBulkSample) > 0) {
            depthBulkSample = (depthBulkSample & sampAB[,depi] >= 20)
        } else {
            depthBulkSample = (sampAB[,depi] >= 20)
        }
    }
    depthBulkSampleGoodi = which(depthBulkSample)

    needr = vector()
    sgHighi = vector()
    for (i in 1:length(samples)) {
        sn = samples[i]
        depi = match(paste(sn, "d", sep=""), colnames(sampAB))
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        nbi = match(paste(sn, "nb", sep=""), colnames(sampAB))
        nbSN1i = match(paste(sn1, "nb", sep=""), colnames(sampAB))
        nbSN2i = match(paste(sn2, "nb", sep=""), colnames(sampAB))
        if (length(needr) == 0) {
            sgHighi = which(as.numeric(sampAB[,mafai]) > 0.15)
            needr = which(as.numeric(sampAB[,depi]) >= 20 & (as.numeric(sampAB[,nbi]) != 0 | as.numeric(sampAB[,nbSN1i]) == 0 | as.numeric(sampAB[,nbSN2i]) == 0))
        } else {
            sgHighi = union(sgHighi, which(as.numeric(sampAB[,mafai]) > 0.15))
            needr = intersect(needr, which(as.numeric(sampAB[,depi]) >= 20 & (as.numeric(sampAB[,nbi]) != 0 | as.numeric(sampAB[,nbSN1i]) == 0 | as.numeric(sampAB[,nbSN2i]) == 0)) )
        }
    }

    sgHighi = intersect(depthBulkSampleGoodi, intersect(needr, sgHighi))
    pubTi = intersect(intersect(needr,depthBulkSampleGoodi), submuts$pubTi)
    ssAi = intersect(intersect(needr,depthBulkSampleGoodi), submuts$ssAi)
    ssBi = intersect(intersect(needr,depthBulkSampleGoodi), submuts$ssBi)
    allSubi = intersect(intersect(needr,depthBulkSampleGoodi), union(submuts$subAi,submuts$subBi))
    sharedi = setdiff(allSubi,union(ssAi,ssBi))
    outlist[[1]] = paste(sampAB$chr[pubTi],sampAB$pos[pubTi],sampAB$ref[pubTi],sampAB$alt[pubTi],sep=":")
    outlist[[2]] = paste(sampAB$chr[sharedi],sampAB$pos[sharedi],sampAB$ref[sharedi],sampAB$alt[sharedi],sep=":")
    outlist[[3]] = paste(sampAB$chr[ssAi],sampAB$pos[ssAi],sampAB$ref[ssAi],sampAB$alt[ssAi],sep=":")
    outlist[[4]] = paste(sampAB$chr[ssBi],sampAB$pos[ssBi],sampAB$ref[ssBi],sampAB$alt[ssBi],sep=":")
    needallr = union(sharedi, union(pubTi,union(ssAi,ssBi)))
    needallr = union(sgHighi,needallr)
    
    for (i in 1:length(samples)) {
        sn = samples[i]
        mafai = match(paste(sn, "mafa", sep=""), colnames(sampAB))
        foundr = which(as.numeric(sampAB[,mafai]) > 0.15)
        sampler = intersect(foundr, needallr)
        outlist[[i+4]] = paste(sampAB$chr[sampler],sampAB$pos[sampler],sampAB$ref[sampler],sampAB$alt[sampler],sep=":")
    }
    
    names(outlist) = c(paste(sn1,sn2,"public",sep="_"),"Pvt_Shared",paste(sn1,"specific",sep="_"),paste(sn2,"specific",sep="_"),samples)
    return(outlist)

}


spruceInput <- function(samp, samples, minMaf, clonal=FALSE) {
    depthCols = paste(samples, "d", sep="")
    licheeRow = !is.na(samp$functionalClass)
    for (i in 1:length(depthCols)) {
        licheeRow = licheeRow & samp[,match(depthCols[i], colnames(samp))] >= 15
    }
    samp = samp[licheeRow,]
    result = data.frame()
    mutindex = 0
    #maxnanbs = 2   #max nanb status
    for (r in 1:dim(samp)[1]) {    #each mutation
        maxMaf = 0
        isClonal = FALSE
        #ntnbs = vector()   #all nt nb stats
        for (i in 1:length(samples)) {
            sn = samples[i]
            mafi = match(paste(sn,"mafa",sep=""), colnames(samp))
            maf = as.numeric(samp[r,mafi])
            ccfi = match(paste(sn,"ccf",sep=""), colnames(samp))
            ccf = as.numeric(samp[r,ccfi])
            ccfSDi = match(paste(sn,"ccfSD",sep=""), colnames(samp))
            ccfSD = as.numeric(samp[r,ccfSDi])
            #nti = match(paste(sn,"nt",sep=""), colnames(samp))
            #nt = as.character(samp[r,nti])
            #nbi = match(paste(sn,"nb",sep=""), colnames(samp))
            #nb = as.character(round(as.numeric(samp[r,nbi])))
            #na = nt-nb
            #nanb = paste(na,nb,sep="_")
            #nanbs = append(nanbs, nanb)
            if (maf > maxMaf) {
                maxMaf = maf
            }
            if ((ccf + 1.96*ccfSD) >= 1) {
                isClonal = TRUE
            }
        }
        if (maxMaf < minMaf)  {
            next
        }
        if (clonal ==TRUE & isClonal == FALSE){
            next
        }

        #nanbst = table(nanbs)
        #nstatus = length(nanbst)
        #if (nstatus > maxnanbs) {
        #    maxnanbs = nstatus     #update max ntnb status
        #}
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
            if (na == 1 & nb == 1){  #WT type diploid
                result = rbind(result, c(i-1, sn, mutindex, mutn, maflb, maf, mafub, 1, 1, 1, NA, NA, NA, seg),
                    stringsAsFactors=FALSE)
            } else {
                result = rbind(result, c(i-1, sn, mutindex, mutn, maflb, maf, mafub, 1, 1, nc, na, nb, pu, seg),
                    stringsAsFactors=FALSE)
            }
        } #each tumor
        mutindex = mutindex+1
    } #each mutation
    colnames(result) = c("sample_index","sample_label","character_index","character_label",
                "vaf_lb","vaf_mean","vaf_ub","x","y","mu","x","y","mu","seg")
    return(result)
}


rAUC <- function(data, lower, mafs, depths) {

    if (length(mafs) <= 8) {
        lower = round((0.08/length(mafs)),2)
    } else {
        lower = max(0.005, round((0.08/length(mafs)),3))
    }
    #message(paste("lower: ", lower,sep=""))
    weightAF = weightAFs(data, mafs, depths)
    #message(paste(weightAF, collapse=" "))
    weightAF = weightAF[which(weightAF >= lower & weightAF <= 0.25)]
    #message(paste("weightedMuts: ", length(weightAF), sep=""))
    vafs = seq(lower,0.25,by=0.01)
    nstep = length(vafs)
    counts = vector()
    ncounts = ((1/vafs)-(1/0.25))/((1/lower)-(1/0.25))
    
    for ( i in 1:length(vafs) ) {
        counts = append(counts, length(which(weightAF > vafs[i]))-length(which(weightAF > vafs[nstep])))
    }
    counts = counts/counts[1]
    bei = bezierCurve(vafs, counts, length(vafs))
    AUC = trapz(bei$x,bei$y)
    #AUC = trapz(mafs,counts)
    nAUC = trapz(vafs,ncounts)
    rAUC = AUC/nAUC
    rAUC = round(rAUC,8)
    return(list(rAUC=rAUC,weightAF=weightAF))

}


weightAFs <- function(af.data, mafs, depths, minAF=0) {
    
    keep = vector()
    for (i in 1:length(mafs)) {
        cmaf = mafs[i]
        cmafi = match( cmaf, colnames(af.data) )
        if (length(keep) == 0) {
            keep = af.data[,cmafi] > minAF
        } else {
            keep = keep | (af.data[,cmafi] > minAF)
        }
    }
    af.data = af.data[keep,]
    #message(paste("keepRows4AUC: ", dim(af.data)[1], sep=""))

    tdepth = rep(0,dim(af.data)[1])
    talt = rep(0,dim(af.data)[1])
    for (i in 1:length(mafs)) {
        cmaf = mafs[i]
        cmafi = match(cmaf,colnames(af.data))
        cdep = depths[i]
        cdepi = match(cdep,colnames(af.data))
        tdepth = tdepth + af.data[,cdepi]
        #message(paste(tdepth, collapse=" "))
        talt = talt + af.data[,cdepi]*af.data[,cmafi]
    }
    tmaf = talt/tdepth
    return(tmaf)
}

plotAFS <- function(muts, sn="sample", power=1, add=FALSE, color="black", lwd=2, lty = 1, up=0.25, down, last=FALSE, countsup=vector(), colup="red") {

    if ( last == TRUE ) {
        curve((1/x-1/up)/(1/down-1/up), from=down, to=up,add=T, col=color, lwd=lwd, lty=lty)
        return(1)
    }

    mafs = seq(down,up,by=0.01)
    nstep = length(mafs)
    #neutral = 1/mafs^power-1/mafs[nstep]^power

    #cumulative fraction
    counts = vector()
    for ( i in 1:length(mafs) ) {
        counts = append(counts, length(which(muts > mafs[i]))-length(which(muts > mafs[nstep])))
    }
    counts = counts/counts[1]

    if ( length(countsup) > 0 ) {
        lines(bezierCurve(mafs, countsup, length(mafs)), col=colup, lwd=lwd, lty=lty)
        return(2)
    }
    
    if (add == FALSE) {
        plot(bezierCurve(mafs,counts,length(mafs)), type="l", lwd=lwd,lty=lty, axes = F, xlab="",xlim=c(down,up),ylim=c(0,1),
             ylab="Fraction SNVs in [f, fmax]", col=color, main=sn, cex.main=1.7,cex.lab=1.6)
        title(xlab="f (merged VAF)", line=4, cex.lab=1.6)
        axis(side=1,at=c(seq(down,0.22,by=0.02),0.25),labels=c(seq(down,0.22,by=0.02),0.25),las=2, cex.axis=1.3)
        axis(side=2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2),cex.axis=1.3)
        #axis(side=1,at=c(1/c(0.25,0.12,0.08,0.06,0.05)-1/0.25),labels=paste("f=",c(0.25,0.12,0.08,0.06,0.05),sep=""))
    }
    else {
        lines(bezierCurve(mafs, counts, length(mafs)), col=color, lwd=lwd, lty=lty)
        #abline(0,-1/(1/down-1/up),lwd=3)
    }
    #ruption.lm = lm(counts ~ neutral)
    #rs = summary(ruption.lm)$r.squared
    #text(1/0.22^power, max(counts), labels = bquote(R^2 == .(rs), list(rs=round(rs,5))))
    return(counts)

}


bezierCurve <- function(x, y, n=10) {
    outx <- NULL
    outy <- NULL

    i <- 1
    for (t in seq(0, 1, length.out=n))
        {
            b <- bez(x, y, t)
            outx[i] <- b$x
            outy[i] <- b$y

            i <- i+1
        }

    return (list(x=outx, y=outy))
}


bez <- function(x, y, t){
    outx <- 0
    outy <- 0
    n <- length(x)-1
    for (i in 0:n)
        {
        outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
        outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
        }

    return (list(x=outx, y=outy))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# plotting simVAF data


#for plotting simulated vaf data

plotRes.simVAF.matrix.pdf <- function(sampAB, samples, depths, pdfsize = 16, plotType = "AF", snr="", sns=vector()) {
    combinations = combn(length(samples),2)
    nplots = dim(combinations)[2]
    ndims = length(samples)-1

    resStats = list()
    outfile = paste(snr,"multi_hist.pdf",sep="")
    if (plotType == "Scatter") {
        outfile = paste(snr,"multi_scatter.pdf",sep="")
    }
    if (plotType == "Density") {
        outfile = paste(snr,"multi_density.pdf",sep="")
    }
    pdf(file = outfile, width=pdfsize, height=pdfsize)
    layout(matrix(seq(1,ndims^2), ndims, ndims, byrow = FALSE))
    #par(mfcol=c(ndims,ndims))
    pindex = 0
    for (ci in 1:ndims) {                           #for each column
        noplotRow = vector()
        if (ci > 1) {
            noplotRow = (1:ndims)[1:(ci-1)]         #which row do not plot?
        }
        for (ri in 1:ndims) {                       #for each row
            if (ri %in% noplotRow) {                #no plot
                plot(NULL,NULL,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="")
            } else {                                #plot
                pindex = pindex + 1
                pair = combinations[,pindex]
                sn1 = samples[pair[1]]              #samplename1
                sn1mafa = sn1
                sn1d = depths[pair[1]]
                sn2 = samples[pair[2]]              #samplename2
                sn2mafa = sn2
                sn2d = depths[pair[2]]
                main.title = paste(snr, sns[pair[1]], "vs", sns[pair[2]], sep=" ")
                statName = paste(snr, sns[pair[1]], sns[pair[2]], "stats", sep="_")
                if (plotType == "AF") {
                    resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                sn1=sn1mafa, sn2=sn2mafa, plotDensity=F, plotScatter=F, pdf=F)
                } else if (plotType == "Scatter") {
                    resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                sn1=sn1mafa, sn2=sn2mafa, plotDensity=F, plotAF=F, pdf=F)
                } else if (plotType == "Density") {
                    resStats[[statName]] = plotRes.simVAF.pdf(sampAB[which(sampAB[,sn1d] >= 15 & sampAB[,sn2d] >= 15),], main.title, main=main.title, sn1n=sns[pair[1]], sn2n=sns[pair[2]], dp1=depths[pair[1]], dp2=depths[pair[2]],
                                sn1=sn1mafa, sn2=sn2mafa, plotScatter=F, plotAF=F, pdf=F)
                }
            }
        }
    }
    dev.off()
    return(1)
}



plotRes.simVAF.pdf <- function(sampAB, sampName, main=sampName, sn1n="A", sn2n="B", sn1="maf1", sn2="maf2", dp1="depth1", dp2="depth2",minAF=0.05, ratio=1, plotAF=TRUE, plotDensity=TRUE, plotScatter=TRUE, pdf=TRUE, alpha=1, binw=0, reportSSR=FALSE) {

  sn1s = sn1n     #name of sample 1
  sn2s = sn2n     #name of sample 2
  
  dp1i = match(dp1,colnames(sampAB))
  dp2i = match(dp2,colnames(sampAB))
  
  maf1Index = match(sn1, colnames(sampAB))
  maf2Index = match(sn2, colnames(sampAB))

  #check depth power to reject a presence of a mutation
  depthPowerKeep <- as.vector(apply(sampAB, 1, function(x,maf1i,maf2i,dp1i,dp2i) {
                                        if(as.numeric(x[maf1i]) == 0){vaf = as.numeric(x[maf2i])
                                                   if (pbinom(0,as.numeric(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else if (as.numeric(x[maf2i]) == 0){vaf = as.numeric(x[maf1i])
                                                   if (pbinom(0,as.numeric(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
                                        else {TRUE}
                                    }, maf1i=maf1Index,maf2i=maf2Index,dp1i=dp1i,dp2i=dp2i))
  sampAB = sampAB[depthPowerKeep,]

  subMuts = subclonalMutSim(sampAB, sn1, sn2, dp1, dp2, minAF=minAF)
  
  allA_Rows = which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1)                                                       #all sample A rows
  subA_Rows = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1))                             #sample A sub rows
  ssA_Rows  = intersect(subMuts$subAi, which(sampAB[,maf1Index] > minAF & sampAB[,maf1Index] <= 1 & sampAB[,maf2Index] == 0))   #sample A specific rows
  
  allB_Rows = which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1)                                                       #all sample B rows
  subB_Rows = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1))                             #sample B sub rows
  ssB_Rows  = intersect(subMuts$subBi, which(sampAB[,maf2Index] > minAF & sampAB[,maf2Index] <= 1 & sampAB[,maf1Index] == 0))   #sample B specific rows
  allAB_Rows = union(allA_Rows, allB_Rows)

  ssR = length(union(ssA_Rows,ssB_Rows))/length(union(subA_Rows,subB_Rows))
  if (reportSSR == TRUE) {
      return(ssR)
  }
  
  BinWidthA = round(dpih(sampAB[allA_Rows, maf1Index]/ratio),2)
  if (BinWidthA == 0) { BinWidthA = 0.02 }
  BinWidthB = round(dpih(sampAB[allB_Rows, maf2Index]/ratio),2)
  if (BinWidthB == 0) { BinWidthB = 0.02 }
  BinWidth = min(c(BinWidthA, BinWidthB, 0.1))
  if ( binw != 0 ) {
      BinWidth = binw
  }
  #message(paste("bin width: ", BinWidthA, BinWidthB, BinWidth, sep=" "))
  nbreaksA = round((max(sampAB[allA_Rows, maf1Index]/ratio)-min(sampAB[allA_Rows, maf1Index]/ratio)+0.01)/BinWidth)
  nbreaksB = round((max(sampAB[allB_Rows, maf2Index]/ratio)-min(sampAB[allA_Rows, maf2Index]/ratio)+0.01)/BinWidth)
  nbreaks = ceiling(diff(range(minAF,1))/BinWidth)
  #message(paste("nbreaks: ", nbreaksA, nbreaksB, nbreaks, sep=" "))
  breaksA = seq(minAF,1,length.out=nbreaks)
  breaksB = seq(minAF,1,length.out=nbreaks)

  
  sampAh = hist(sampAB[allA_Rows, maf1Index]/ratio, breaks=breaksA,  plot=F)
  sampAhsub = hist(sampAB[subA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  sampAhss = hist(sampAB[ssA_Rows, maf1Index]/ratio, breaks=sampAh$breaks, plot=F)
  ylimup = max(sampAh$count)
  

  sampBh = hist(sampAB[allB_Rows, maf2Index]/ratio, breaks=breaksB, plot=F)
  sampBh$counts = sampBh$counts*(-1)
  sampBhsub = hist(sampAB[subB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhsub$counts = sampBhsub$count*(-1)
  sampBhss = hist(sampAB[ssB_Rows, maf2Index]/ratio, breaks=sampBh$breaks, plot=F)
  sampBhss$counts = sampBhss$counts*(-1)
  ylimdown = min(sampBh$count)

  if (abs(ylimup) >= abs(ylimdown)) {
      ylimdown = (-1)*ylimup
  } else {
      ylimup = (-1)*ylimdown
  }
  
  if (plotAF == TRUE) {

      if (pdf == TRUE) {
          pdf(file = paste(sampName, "hist.pdf", sep="_"), width = 8, height = 8, useDingbats=FALSE)
      }
      par(mar=c(4.5,5,4.5,0))

      plot( sampAh, col=rgb(0,0,0,1/4), xlim=c(0, 1), ylim=c(ylimdown,ylimup), border=F, ylab="# of Mutations", xlab="VAF", axes = F, main = main,cex.lab = 2.9, cex.main = 2.9)    # first histogram
      plot( sampAhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green no border
      plot( sampAhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue

      
      plot( sampBh,col=rgb(0,0,0,1/4), border=F, add=T )  # second histogram
      plot( sampBhsub, col=rgb(178/255,223/255,138/255,1), add=T, border=F )    #subclonal green
      plot( sampBhss, col=rgb(31/255,120/255,180/255,1), add=T, border=F )      #site specific blue
      sampAhss$counts = 0                                  #make black line
      sampBhss$counts = 0                                  #make black line
      plot( sampAhss, col="black", add=T, border=F)        #make black line
      plot( sampBhss, col="black", add=T, border=F)        #make black line
      
      axis(side=1,at=seq(0,1,by=0.1),labels=seq(0,1,by=0.1),cex.axis=2.4)       #x-axis
      axis(side=2,at=decideTickAt(ylimdown, ylimup),labels=allAbs(decideTickAt(ylimdown, ylimup)),cex.axis=2.4)
      text(x=0.5,y=ylimdown,labels=paste(sn2s,",",length(which(sampAB[,maf2Index] > minAF)), "SSNVs", sep = " "), cex=2.8)
      text(x=0.5,y=ylimup,labels=paste(sn1s,",",length(which(sampAB[,maf1Index] > minAF)), "SSNVs", sep = " "), cex=2.8)

      #stats starting from here!
      fHsub = round(mean(c(subMuts$ratioHighSubA,subMuts$ratioHighSubB)),3)
      fst = round(subMuts$FST,3)
      ksd = round(subMuts$KSD,3)
      text(x=0.8,y=3*ylimup/4, labels=bquote(paste("fH"["sub"], " = ", .(fHsub))),cex=2.6)
      text(x=0.8,y=(3/4-0.167)*ylimup, labels=bquote(paste("FST", " = ", .(fst))),cex=2.6)
      text(x=0.8,y=(3/4-0.334)*ylimup, labels=bquote(paste("KSD", " = ", .(ksd))),cex=2.6)

      npub = subMuts$pubTn
      legend(0.5,ylimdown/4, legend=c(paste("Public ","(",npub,")",sep=""),"Pvt-Shared","Rgn Specific"),
             col=c(rgb(0,0,0,1/4),rgb(178/255,223/255,138/255,1),rgb(31/255,120/255,180/255,1)), pch=15, bty="n", cex=2.4)
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotDensity == TRUE) {
      #two way density plot
      if (pdf == TRUE) {
          pdf(file = paste(sampName, "density.pdf", sep="_"), width = 8, height = 8, useDingbats=FALSE)
      }
      par(mar=c(4.5,5,4.5,2))
      if (pdf == FALSE) {
          par(mar=c(5,5,5,4))
      }
      smkey(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("VAF",sn1s,sep=" "), ylab=paste("VAF",sn2s,sep=" "),
            main = main, xlim=c(0,1), ylim=c(0,1),cex.lab = 2.3, cex.main = 2.3,cex.axis=1.7)
      segments(0,0,0.8,0.8)
      if (pdf == TRUE) {
          dev.off()
      }
  }

  if (plotScatter == TRUE) {
      #scatter plot with density color
      if (pdf == TRUE) {
          pdf(file = paste(sampName, "scatter.pdf", sep="_"), width = 5, height = 5, useDingbats=FALSE)
      }
      pub_Rows = setdiff(allAB_Rows, union(subA_Rows, subB_Rows))
      pub_Rows = match(pub_Rows, allAB_Rows)
      shared_Rows = setdiff(union(subA_Rows, subB_Rows), union(ssA_Rows, ssB_Rows))
      shared_Rows = match(shared_Rows, allAB_Rows)
      ssAB_Rows = union(ssA_Rows, ssB_Rows)
      ssAB_Rows = match(ssAB_Rows, allAB_Rows)
      allSub_Rows = match(union(subA_Rows, subB_Rows), allAB_Rows)
      if (pdf == TRUE) {
          scatterDensityPlot(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("MAF",sn1s,sep=" "), ylab=paste("MAF",sn2s,sep=" "),
                             main = main, cex=1.2, cex.lab = 2.3, cex.main = 2.3, cex.axis=1.7,
                             groups=list(a=pub_Rows, b=shared_Rows, c=ssAB_Rows),alpha=alpha,
                             groupColors=list(a=brewer.pal(9, "Greys")[3:9], b=brewer.pal(9, "Greens")[2:5], c=brewer.pal(9, "Blues")[3:9]))
      } else {
          scatterDensityPlot(sampAB[allAB_Rows,maf1Index],sampAB[allAB_Rows,maf2Index],xlab=paste("MAF",sn1s,sep=" "), ylab=paste("MAF",sn2s,sep=" "),
                             main = main, cex=1.2, cex.lab = 2.3, cex.main = 2.3, cex.axis=1.7, layout = FALSE,
                             groups=list(a=pub_Rows, b=shared_Rows, c=ssAB_Rows),alpha=alpha,
                             groupColors=list(a=brewer.pal(9, "Greys")[3:9], b=brewer.pal(9, "Greens")[2:5], c=brewer.pal(9, "Blues")[3:9]))
      }
      if (pdf == TRUE) {
          dev.off()
      }
  }

  return(subMuts)
}



subclonalMutSim <- function(sampAB, snA, snB, dpA, dpB, minAF=0.05, statsAF=0.08, highAF=0.2, ratio=1, pob="pubOrSub")   {                  #determinine subclonal mutations
    mafaAi = match(snA, colnames(sampAB))
    mafaBi = match(snB, colnames(sampAB))
    depthAi = match(dpA, colnames(sampAB))
    depthBi = match(dpB, colnames(sampAB))
    pobi = match(pob, colnames(sampAB))
    
    # for KSD
    subAi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > minAF )    
    ssAi  = intersect(subAi, which( sampAB[,mafaAi] > minAF & sampAB[,mafaBi] == 0 ))
    mutsA = sampAB[subAi,mafaAi]/ratio
    
    subBi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaBi] > minAF )
    ssBi  = intersect(subBi, which( sampAB[,mafaBi] > minAF & sampAB[,mafaAi] == 0 ))
    mutsB = sampAB[subBi,mafaBi]/ratio
    
    KSD = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)

    # pub
    pubTi = which( sampAB[,pobi] == "public" )

    # FST
    allSubRows = union(subAi,subBi)
    mutsSub = sampAB[allSubRows,]
    mutsSub = data.frame( maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi] )
    FST = fst.hudson(mutsSub, minAF=statsAF)

    # for other stats
    mutsA2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF )), mafaAi]
    mutsAh2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF )), mafaAi]
    mutsASr2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] > 0.02)), mafaAi]   #shared
    
    mutsASp2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > statsAF & sampAB[,mafaBi] == 0)), mafaAi]
    mutsASph2 = sampAB[intersect(subAi, which( sampAB[,mafaAi] > highAF & sampAB[,mafaBi] == 0)), mafaAi]
    
    mutsB2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF )), mafaBi]
    mutsBh2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF )), mafaBi]
    mutsBSr2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] > 0.02)), mafaBi]  #shared
    
    mutsBSp2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > statsAF & sampAB[,mafaAi] == 0)), mafaBi]
    mutsBSph2 = sampAB[intersect(subBi, which( sampAB[,mafaBi] > highAF & sampAB[,mafaAi] == 0)), mafaBi]

    # mutational counts
    pubTn = length(pubTi)
    ssAn = length(ssAi)
    ssBn = length(ssBi)
    sharedAn = length(mutsASr2)
    sharedBn = length(mutsBSr2)
        
    # list for output    
    muts = list(A=mutsA,B=mutsB,subAi=subAi,subBi=subBi, ssAi=ssAi, ssBi=ssBi, pubTi=pubTi, #subArow=subArow, subBrow=subBrow, pubTrow = pubTrow,
        lenSubA=length(mutsA2),lenSubAh=length(mutsAh2),ratioHighSubA=length(mutsAh2)/length(mutsA2),
        lenSubB=length(mutsB2),lenSubBh=length(mutsBh2),ratioHighSubB=length(mutsBh2)/length(mutsB2),
        lenSsA=length(mutsASp2),lenHighSsA=length(mutsASph2),ratioHighSsA=length(mutsASph2)/length(mutsASp2),pSsA=length(mutsASp2)/length(mutsA2),
        lenSsB=length(mutsBSp2),lenHighSsB=length(mutsBSph2),ratioHighSsB=length(mutsBSph2)/length(mutsBSp2),pSsB=length(mutsBSp2)/length(mutsB2),
        lenSharedA=length(mutsASr2), lenSharedB=length(mutsBSr2), ratioSharedA=length(mutsASr2)/length(mutsA2), ratioSharedB=length(mutsBSr2)/length(mutsB2),
        FST=FST, KSD=KSD, pubTn=pubTn, sharedAn=sharedAn, sharedBn=sharedBn, ssAn=ssAn, ssBn=ssBn)

    return(muts)
}


QT95 <- function(gd){
    qd = quantile(gd, prob=seq(0,1,0.025))
    qdl = round(as.numeric(qd["2.5%"]),2)
    qdh = round(as.numeric(qd["97.5%"]),2)
    return(c(qdl, qdh))
}

QT90 <- function(gd){
    qd = quantile(gd, prob=seq(0,1,0.025))
    qdl = round(as.numeric(qd["5%"]),2)
    qdh = round(as.numeric(qd["95%"]),2)
    return(c(qdl, qdh))
}

samplemean <- function(x, d) {
  return(mean(x[d]))
}


outMutTable <- function(data, samples) {
    if (length(samples) == 1) {
        keep = data[,paste(samples,"ccf",sep="")] > 0
    } else {
        keep = rowSums(data[,paste(samples,"ccf",sep="")]) > 0
    }
    data2 = data[keep,]
    #data2 = data[!(data$pubOrSub == "unknown"),]
    #data2 = data2[,c("chr","pos","ref","alt","geneName","geneLoc","functionalClass",
    #as.vector(t(outer(samples,c("mafc","mafa","ccf","ccfSD","refc","altc","pu","pa","nt","nb"), paste, sep=""))),
    #    colnames(data)[grepl("pubOrSub", colnames(data2))])]
    data2 = data2[,c("chr","pos","ref","alt","geneName","geneLoc","functionalClass",
    as.vector(t(outer(samples,c("mafc","ccf","ccfSD","refc","altc","pu"), paste, sep=""))))]
    return(data2)
}

outTable.maf <- function(data, samples) {
    outmaf = data.frame()
    for(r in 1:dim(data)[1]) {                                          #each row
        Hugo_Symbol = gsub("\\(dist=\\d+\\)", "", data[r,"geneName"])
        Chromosome = data[r,"chr"]
        Start_Position = data[r,"pos"]
        End_Position = data[r,"pos"]
        Variant_Classification = data[r,"functionalClass"]
        if (is.na(Variant_Classification)) {
            Variant_Classification = data[r,"geneLoc"]
        }
        Variant_Type = "DEL"
        Reference_Allele = data[r,"ref"]
        if (Reference_Allele == '-') {
            Variant_Type = "INS"
        }
        Tumor_Seq_Allele2 = data[r,"alt"]
        
        for (i in 1:length(samples)) {
            tmp = data.frame()
            sn = samples[i]
            Tumor_Sample_Barcode = sn
            ccfindex = match(paste(sn, "ccf", sep=""), colnames(data))
            ccfsdindex = match(paste(sn, "ccfSD", sep=""), colnames(data))
            somindex = match("traceSomatic", colnames(data))
            if ((data[r,ccfindex] + 1.96*data[r,ccfsdindex]) >= 1) {
                tmp = rbind(tmp, c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Variant_Type,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode))
                colnames(tmp) =  c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
                outmaf = rbind(outmaf, tmp)
            }
        }
    }
    return(outmaf)
}

