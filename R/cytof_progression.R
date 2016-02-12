#' Progression estimation of cytof expression data 
#' 
#' Apply isomap to estimate the order of cell progression
#' 
#' @param data Expression data matrix
#' @param method A vector of cluster results for the \code{data}
#' @param method Method for estimation of cell progression, \code{isomap} (by default), \code(tsne) or \code(pca)
#' @param uniformClusterSize the down sampled size of each cluster
#' @param seed the seed for random down sample of the clusters
#' @return a list of the progression data including \code{sampleData}, \code{sampleCluster} and \code{progressionData}.
#' @author Chen Jinmiao
#' @export
#' @examples
#' data(iris)
#' in_data <- iris[, 1:4]
#' out_data <- cytof_progression(in_data, cluster = iris[,5], uniformClusterSize = 50)
cytof_progression <- function(data, cluster, method="isomap", uniformClusterSize = 500, seed = 500){
    
    data <- as.matrix(data)
    
    if(!is.null(uniformClusterSize)){
        clusterStat <- as.data.frame(table(cluster))
        uniformClusterSize <- min(uniformClusterSize, min(clusterStat$Freq))
        
        set.seed(seed)
        clusterList <- split(1:length(cluster), factor(cluster))
        eachCluster <- llply(clusterList, function(x) sample(x, uniformClusterSize, replace = FALSE))
        sampleCellID <- do.call(base::c, eachCluster)
        
        nCluster <- rep(1:length(unique(cluster)), each = uniformClusterSize)
        sampleData <- data[sampleCellID, ,drop=FALSE]
    }else{
        sampleData <- data
        nCluster <- cluster
    }
    
    if(method == "isomap"){
        progressionData <- cytof_dimReduction(sampleData, method = "isomap", out_dim = 2)
    }else if(method == "tsne"){
        progressionData <- cytof_dimReduction(sampleData, method = "tsne", out_dim = 2)
    }else if(method == "pca"){
        progressionData <- cytof_dimReduction(sampleData, method = "pca", out_dim = 2)
    }else{
        return(NULL)
    }
    
    progressRes <- list(sampleData = sampleData, 
                        sampleCluster = nCluster, 
                        progressionData = progressionData)
    return(progressRes)
}