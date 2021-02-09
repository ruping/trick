computeCCF <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90) {

    ccf = 0
    ccf2 = 0
    sd = 0
    cc <- seq(0.02, 1, by = 0.01)
    evoType = "A1/A2/B/C"
    N = A + S
    nc = nt * pa + 2 * (1 - pa)
    nc2 = nt * sAGP + 2 * (1 - sAGP)
    
    #message(paste(c(f,A,S,pu,pa,sAGP,nt,nb),collapse=" "))
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
        ccf = (f*nc2)/sAGP
        ff.C = cc[cc<=(1-sAGP)]/nc2
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
            ccf = (f/pu)*nc - (nt - nb - 1)*pa          #early A1
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            #ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc2
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            #ccf2 = ccf2/pu
            sd <- Ms.A$SD
        #} else if (maxType == "pEarly.b" & prior != "late") {
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

effectivePurity <- function(p, Nt) {
    pu = Nt*p/(Nt*p + 2*(1-p))
    return(pu)
}
