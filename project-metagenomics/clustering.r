
datafile <- "/groups/bioc6243/data/csv/MetaHIT_SangerSamples.genus.txt"
data=read.table(datafile, header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

clsfile <- "/groups/bioc6243/data/csv/MetaHIT_SangerSamples.genus.cls"
datacls=read.table(clsfile, header=T, row.names=1, dec=".", sep=",")


KLD <- function(x,y) sum(x * log(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) { 
            resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))
        }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
}


data.dist = dist.JSD(data)
pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
    return(cluster)
}

data.cluster=pam.clustering(data.dist, k=3)

#obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
#cat(obs.silhouette) #0.1899451

#data=noise.removal(data, percent=0.01)
#cat(data.cluster)

## plot 1
require(ade4)
k = 3

png('bet.png')
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
#dev.new()
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")
#s.class(obs.bet$ls, fac=as.factor(datacls$cls), grid=F,sub="Between-class analysis")


#plot 2
png('pca.png')
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
#dev.new()
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")
#s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))


