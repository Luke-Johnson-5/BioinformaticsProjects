
datafile <- "/groups/bioc6243/data/csv/MetaHIT_SangerSamples.genus.txt"
data=read.table(datafile, header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

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


require(clusterSim)
maxclusters = 20

nclusters = NULL
for (k in 1:maxclusters) {
    if (k==1) {
        nclusters[k]=NA 
    } else {
        data.cluster_temp=pam.clustering(data.dist, k)
        nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist, centrotypes = "medoids")
    }
}


png('sizes.png')
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
#dev.off()



