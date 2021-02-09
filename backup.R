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
