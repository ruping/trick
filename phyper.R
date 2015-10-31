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

setwd(path)
d = read.delim(file, header=F)
p = apply(d, 1, function(x){ phyper(x[3],x[1],(22000-x[1]),x[2],lower.tail=F) })
p.adj = p.adjust(p, method="BH")
e = data.frame(p = p, p.adj = p.adj)
write.table(e, file="pvalues", quote = F, sep = "\t", row.names = F, col.names = F)
