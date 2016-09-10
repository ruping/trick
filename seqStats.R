library(maptools)
library(RColorBrewer)


###########################################################
# cumulative target coverage                            ###
###########################################################

plotLorenz <- function(need, outFile, legend, otype="pdf",plotCov = TRUE, plotLorenz = TRUE, cex.main=1.7, cex.lab=1.5, cex.axis=1.3, targetName="Targeted Region") {

    colorGroup = decideColorGroup(legend, FALSE)
    if (length(legend) > 22) {
        colorGroup = decideColorGroup(legend)
        #stop("length of legend is greater than 22")
    }

    if (otype == "pdf"){
        pdf(file=outFile, width=15, height=7)
    } else if (otype == "png") {
        png(file=outFile, width=1200, height=580)
    }
    if (plotCov & plotLorenz){
        layout(matrix(c(1,2),ncol=2))
    }
    if (plotCov){
        plot(NULL, xlim=c(0, 5),ylim=c(0,1), type="l", axes=F, xlab=expression("Read Depth">=""), ylab=paste("Cumulative Fraction of",targetName,sep=" "), main=paste("Depth Coverage in",targetName,sep=" "), cex.lab=cex.lab, cex.main=cex.main)
        axis(side=1, at=c(0,1,log10(30),2,3,4,5), labels=c(1,10,30,100,expression(10^3),expression(10^4),expression(10^5)), tck=0.02, cex.axis=cex.axis)
        axis(side=2, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), tck=0.02, cex.axis=cex.axis)
        abline(v=log10(30),lty=2)

        tlabels = data.frame(xat=3.7, yat=rep(1,length(need)))
        rownames(tlabels) = gsub(".lorenzNoDup","",basename(need))
        rownames(tlabels) = gsub("CRCTumor","",rownames(tlabels))
        #rownames(tlabels) = gsub("LOVO_","",rownames(tlabels))
        #rownames(tlabels) = gsub("HCT116_","",rownames(tlabels))
        for (i in 1:length(need)) {
            d = read.delim(need[i])
            colnames(d) = c("dep","depc","depr","cumc","cumr","cumc2","cumr2")
            lines(log10(d$dep), d$cumc, col=makecolor(basename(need[i]), colorGroup))
            tlabels$yat[i] = d$cumc[30]     #a line at depth of 30
        }
        tlabels = tlabels[order(tlabels$yat),]
        tlabels$yatNew = seq(0.2,0.9,by=(0.9-0.2)/(length(need)-1))
        llabels = gsub("HCT116_|LOVO_|WES-|TCGA-\\d+\\-|Patient","",rownames(tlabels))
        llabels = gsub("Recurrence","Rec",llabels)
        pointLabel(tlabels$xat, tlabels$yatNew, label= llabels, cex=cex.axis, col=makecolor(rownames(tlabels), colorGroup))
    }
    if (plotLorenz){
        plot(NULL, xlim=c(0,1),ylim=c(0,1), type="l", axes=F, xlab="Cumulative Fraction of Covered Region", ylab="Cumulative Fraction of Reads",main="Lorenz Curve of Read Consumption", cex.lab=cex.lab, cex.main=cex.main)
        axis(side=1, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), tck=0.02, cex.axis=cex.axis)
        axis(side=2, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), tck=0.02, cex.axis=cex.axis)
        abline(0,1,lty=2)

        sgini = 0
        for (i in 1:length(need)) {
            d = read.delim(need[i])
            colnames(d) = c("dep","depc","depr","cumc","cumr","cumc2","cumr2")
            lines(d$cumc2/max(d$cumc2),d$cumr2,col=makecolor(basename(need[i]), colorGroup))
            sgini = sgini + seqGini(d$cumc2/max(d$cumc2), d$cumr2)
        }
        sgini = round(sgini/length(need), 4)
        legend("topleft", legend=gsub("WES-|TCGA-\\d+\\-|Patient","",legend), col=makecolor(legend, colorGroup), bty="n",lwd=1,cex=cex.axis)
        text(0.5,0.5, labels = paste("Mean Gini Index = ", sgini, sep=""), cex=cex.axis)
    }
    if (otype != "none"){
        dev.off()
    }

}


###########################################################
# median coverage and duplicates                        ###
###########################################################

medianCoverage <- function (needMapStats, needLorenz, ofile, legend, otype = "pdf", cex.main=1.7, cex.lab=1.5, cex.axis=1.3) {

    colorGroup = decideColorGroup(legend, FALSE)
    if (length(legend) > 22){
        colorGroup = decideColorGroup(legend)
        #stop("length of legend is greater than 22")
    }
    
    if (otype == "pdf") {
        pdf(file=ofile,width=12,height=6)
    } else if (otype == "png") {
        png(file=ofile,width=1000,height=500)
    }
    par(mar=c(4.5,4.5,1.5,1.5))
    
    coverage = vector()
    duprate = vector()
    for (i in 1:length(needLorenz)) {
        samp = gsub(".lorenzNoDup","",basename(needLorenz[i]))
        samp = gsub("CRCTumor","",samp)
        #samp = gsub("LOVO_","",samp)
        #samp = gsub("HCT116_","",samp)
        d = read.delim(needLorenz[i])
        colnames(d) = c("dep","depc","depr","cumc","cumr","cumc2","cumr2")
        coverage = append(coverage, d[nearestIndex(d$cumc, 0.5),1])
        names(coverage)[i] = samp
    }
    for (i in 1:length(needMapStats)) {
        samp = gsub(".mapping.stats", "", basename(needMapStats[i]))
        samp = gsub("CRCTumor","",samp)
        #samp = gsub("LOVO_","",samp)
        #samp = gsub("HCT116_","",samp)
        d = read.table(needMapStats[i], header=F)
        mm = d[6,2]/d[2,2]
        duprate = append(duprate, mm)
        names(duprate)[i] = samp
    }
    covdup = data.frame(coverage=coverage, duprate=duprate[match(names(coverage),names(duprate))])

    plot(covdup$duprate, covdup$coverage, col=makecolor(rownames(covdup), colorGroup), pch=19, cex=1.2, xlim=c(-.15,1), ylim=c(0,nearestDecimal(max(coverage))+20),
         xlab="Read Duplication Rate", ylab="Median Coverage", main="Median Coverage vs Duplication Rate", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
    llabels = gsub("HCT116_|LOVO_|WES-|TCGA-\\d+\\-|Patient","",rownames(covdup))
    llabels = gsub("Recurrence","Rec",llabels)
    pointLabel(covdup$duprate, covdup$coverage,label= llabels, col=makecolor(rownames(covdup), colorGroup), cex=1.2)

    legend("topright", legend=gsub("WES-|TCGA-\\d+\\-|Patient","",legend), pch=19, col=makecolor(legend, colorGroup), bty="n",cex=1.5)
    if (otype != "none") {
        dev.off()
    }

}


###########################################################
#  Insert Size box plot                                 ###
###########################################################

plotInsBox <- function(needBox, ofile, legend, otype = "pdf", cex.axis=1.3, cex.main=1.7, cex.lab=1.5, xlab="Samples") {

    colorGroup = decideColorGroup(legend, FALSE)
    if (length(legend) > 22){
        colorGroup = decideColorGroup(legend)
        #stop("length of legend is greater than 22")
    }

    box = vector()
    for (i in 1:length(needBox)) {
        load(needBox[i])
        if (i == 1) {
            box = bb$stats
        } else {
            box = cbind(box, bb$stats)
        }
    }
    
    colnames(box) = gsub(".ins.rda", "", basename(needBox))
    colnames(box) = gsub("CRCTumor","",colnames(box))
    #colnames(box) = gsub("LOVO_","",colnames(box))
    #colnames(box) = gsub("HCT116_","",colnames(box))
    boxor = box[, order(box[3,], decreasing=T)]
    
    if (otype == "pdf") {
        pdf(file = ofile, width=18, height=7)
    } else if (otype == "png") {
        png(file = ofile, width=1440, height=560)
    }
    options(mar=c(0,5,1,2))
    boxplot(boxor, axes=F, ylim = c(-220, 500), main="Insert Size", ylab="bp", xlab=xlab, col=makecolor(colnames(boxor), colorGroup), cex.main=cex.main, cex.lab=cex.lab)
    axis(side=2, at=seq(0,500,by=50), labels=seq(0,500,by=50), cex.axis=1.3)
    llabels = gsub("HCT116_|LOVO_|WES-|TCGA-\\d+\\-|Patient","",colnames(boxor))
    llabels = gsub("Recurrence","Rec",llabels)
    text(1:length(needBox), -120, label= llabels, col = makecolor(colnames(boxor), colorGroup), srt=90, cex=cex.axis)
    #legend("topright", legend=legend, pch=15, col=makecolor(legend, colorGroup), bty="n", cex=cex.axis)
    if (otype != "none"){
        dev.off()
    }

}

###########################################################
#  Xeno mouse proportion                                ###
###########################################################

propMouseBar <- function(needpropm, ofile, legend, otype = "pdf", cex.main=1.7, cex.lab=1.5, cex.axis=1.3) {
    
    colorGroup = decideColorGroup(legend, FALSE)
    if (length(legend) > 22){
        colorGroup = decideColorGroup(legend)
        #stop("length of legend is greater than 22")
    }
    
    propm = vector()
    for (i in 1:length(needpropm)) {
        samp = gsub(".merged", "", basename(needpropm[i]))
        samp = gsub(".xenoStats", "", samp)
        #samp = gsub("LOVO_","",samp)
        #samp = gsub("HCT116_","",samp)
        d = read.table(needpropm[i], header=F)
        mi = which(as.character(d$V1) == "mouseReads:")
        hi = which(as.character(d$V1) == "humanReads:")
        propm = append(propm, d[mi,2]/(d[mi,2]+d[hi,2]))
        names(propm)[i] = samp
    }

    if (otype == "pdf") {
        pdf(file = ofile, width=10, height=6)
    } else if (otype == "png") {
        png(file = ofile, width=800, height=500)
    }
    options(mar=c(0,4,1,2))
    pmb = barplot(propm, col=makecolor(names(propm), colorGroup), axes=F, axisnames=F, ylim=c(-0.12, 0.2), main="Proportion Mouse Reads", cex.main=cex.main)
    axis(side=2, at=seq(0,0.2,by=0.05), labels=paste(seq(0,20,by=5),"%",sep=""),cex.axis=cex.axis)
    llabels = gsub("HCT116_|LOVO_|WES-|TCGA-\\d+\\-|Patient","",names(propm))
    llabels = gsub("Recurrence","Rec",llabels)
    text(pmb[,1], -0.06, label=llabels, col=makecolor(names(propm),colorGroup),srt=90, cex=cex.axis)
    legend("topright", legend=gsub("WES-|TCGA-\\d+\\-|Patient","",legend), pch=15, col=makecolor(legend, colorGroup), bty="n", cex=cex.axis)
    if (otype != "none"){
        dev.off()
    }
    
}


###########################################################
#  Purity and ploidy plot                               ###
###########################################################

purityPloidy <- function(needpp, ofile, legend, otype = "pdf", excludingSamples = c("none"), cex.main=1.7, cex.lab=1.5, cex.axis=1.3) {

    colorGroup = decideColorGroup(legend, FALSE)
    if (length(legend) > 22) {
        #stop("length of legend is greater than 22")
        colorGroup = decideColorGroup(legend, TRUE)
    }
    
    pp = data.frame(purity=rep(1,length(needpp)), ploidy=rep(1,length(needpp)))
    if (excludingSamples[1] != "none") {
        pp = data.frame(purity=rep(1,length(needpp)-length(excludingSamples)), ploidy=rep(1,length(needpp)-length(excludingSamples)))
    }
    
    #clear samples
    needppNew = vector()
    for (i in 1:length(needpp)) {
        samp = gsub("_nclones1.TitanCNA.segments.txt", "", basename(needpp[i]))
        if (samp %in% excludingSamples) {
            next
        } else {
            needppNew = append(needppNew, needpp[i])
        }
    }
    message(length(needppNew))
    
    for (i in 1:length(needppNew)) {
        samp = gsub("_nclones1.TitanCNA.segments.txt", "", basename(needppNew[i]))
        samp = gsub("CRCTumor","",samp)
        #samp = gsub("LOVO_","",samp)
        #samp = gsub("HCT116_","",samp)
        d = read.delim(needppNew[i], header=T)
        sploidy = d$ploidy[1]
        pa = max(as.numeric(na.omit(d$cellularprevalence)))
        con = d$normalproportion[1]
        spurity = 1 - (con + (1-pa)*(1-con))
        pp$purity[i] = spurity
        pp$ploidy[i] = sploidy
        rownames(pp)[i] = samp
    }
    if (otype == "pdf") {
        pdf(file = ofile, width=8, height=8)
    } else if (otype == "png") {
        png(file = ofile, width=700, height=700)
    }

    plot(pp$ploidy, pp$purity, col=makecolor(rownames(pp), colorGroup), ylab="Purity", xlab="Ploidy",main="Purity and Ploidy Estimates", pch=19,cex=1.3,xlim=c(1.5,4),ylim=c(0.1,1), cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
    llabels = gsub("HCT116_|LOVO_|WES-|TCGA-\\d+\\-|Patient","",rownames(pp))
    llabels = gsub("Recurrence","Rec",llabels)
    pointLabel(pp$ploidy, pp$purity, labels=llabels, col=makecolor(rownames(pp), colorGroup), cex=1.2)
    abline(h=c(0.5,0.6),lty=2)
    abline(v=c(2),lty=3)
    if (otype != "none"){
        dev.off()
    }
}



###########################################################
#  utils                                                ###
###########################################################

decideColorGroup <- function(legend, random=TRUE) {

    colorGroup = data.frame(legend=legend, col=rep(rgb(0,0,0,1/2), length(legend)))
    if (random == TRUE) {
        allColors = sampleColors(length(legend))
    } else {
        allColors = c(brewer.pal(9, "Set1")[c(1:5,7:9)], brewer.pal(8, "Dark2")[c(1:4,7)], brewer.pal(11, "Spectral")[1],
            brewer.pal(11, "BrBG")[c(1,11)], brewer.pal(9, "YlGnBu")[8], brewer.pal(9, "YlGnBu")[8], brewer.pal(8, "Set2")[1:4], "black")
        allColors = allColors[1:length(legend)]
    }
    colorGroup$col = apply(col2rgb(allColors),2,function(x){rgb(x[1]/255,x[2]/255,x[3]/255,3/4)})
    return(colorGroup)

}


makecolor <- function(x, colorGroup) {

    color = rep(rgb(0,0,0,1/3), length(x))
    for (i in 1:length(x)) {
        if (class(x) == "integer"){
            current = x[i]
            cc = as.vector(apply( colorGroup, 1, function(x, name=current){if(x[1] == current){x[2]} else {NA}} ))
        } else {
            current = gsub(".lorenzNoDup","",x[i])
            cc = as.vector(apply( colorGroup, 1, function(x, name=current){if(grepl(x[1], current)){x[2]} else {NA}} ))
        }
        cc = cc[!is.na(cc)]
        if (! is.na(cc[1])) {
            color[i] = cc[1]
        }
    }
    return(color)

}


sampleColors <- function(n) {

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist( mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals) ))

    if (n > length(col_vector)) {
        col_vector = rep(col_vector, 2)
    }
    return(sample(col_vector, n))

}


seqGini <- function (p, L) {
    prob = seq(0,1,0.01)      #equal interval
    A = 0
    B = 0
    for (i in 1:length(prob)) {
        ind = nearestIndex(p, prob[i])
        A = A + (p[ind]-L[ind])
        B = B + p[ind]
    }
    return(round(A/B, 4))
}

nearestIndex <- function (p, v) {
    diff = abs(v-p)
    index = which.min(diff)
    return(index)
}


nearestDecimal <- function(x){
    r = x %% 10
    nearD = 0
    if (r > 5){
        nearD = x + (10-r)
    } else {
        nearD = x - r
    }
    return(nearD)
}