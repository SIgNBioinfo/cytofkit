#' Fast clustering by automaticly search and find of density peaks 
#' 
#' This package implement the clustering algorithm described by Alex Rodriguez
#' and Alessandro Laio (2014) with improvements of automatic peak detection and 
#' parallel implementation
#' 
#' @param data a data matrix for clustering
#' @param dimReduction dimenionality reduciton method
#' @param outDim number of dimensions will be used
#' @param distMethod distance calculaiton method, by default is euclidean
#' @param dc distance cutoff value
#' @param gaussian if apply gaussian to esitmate the density
#' @param alpha signance level for peak detection
#' @param detectHalos if detect the halos
#' @param parallel if run the algorithm in parallel
#' @param number of cores umployed for parallel compution
#' 
#' @return a object of \code{densityClusterX} class
#' 
#' @import doParallel
#' @import pdist pdist
#' @importFrom plyr llply
#' @export
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' densityClustRes <- densityClustX(data)
ClusterX <- function(data, dimReduction = NULL, outDim=2, distMethod = "euclidean", 
                          dc, gaussian=TRUE, alpha = 0.001, detectHalos = FALSE, 
                          parallel = FALSE, nCore = 4) {
    if(!is.null(dimReduction))
        data <- dimReduction(data, method = dimReduction, outDim = outDim)
    if(missing(dc))
        dc <- estimateDc(data, sampleSize = 10000, distMethod)
    if(parallel){
        require(doParallel)
        cl <- makeCluster(nCore)  
        registerDoParallel(cl)
    }
    rho <- localDensity(data, distMethod, dc, gaussian=gaussian, ifParallel = parallel)
    deltaWid <- minDistToHigher(data, rho, distMethod, ifParallel = parallel)
    delta <- deltaWid[[1]]
    higherID <- deltaWid[[2]]
    peakID <- peakDetect(rho, delta, alpha)
    cluster <- clusterAssign(peakID, higherID, rho)
    clusTable <- as.vector(table(cluster))
    if(sum(clusTable < length(rho)*0.0005) > 0){
        cat("Noise cluster removing, ")
        peakID <- peakID[clusTable >= length(rho)*alpha]
        cluster <- clusterAssign(peakID, higherID, rho)
    }
    
    if(detectHalos){
        halo <- haloDetect(data, rho, distMethod, cluster, peakID, dc)
    }else{
        halo <- NULL
    }
    
    if(parallel){
        stopCluster(cl)
    }
    
    if(ncol(data) < 3){
        plotData <- data
    }else{
        plotData <- data[ ,c(1,2)]
    }
        
    res <- list(cluster = cluster, dc = dc, rho = rho, delta = delta, peakID = peakID, 
                higherID = higherID, halo = halo, plotData = plotData)
    
    class(res) <- 'densityClusterX'
    return(res)
}


saveRes <- function(x, fileName, ...){
    UseMethod("saveRes", x)
}

saveRes.densityClusterX <- function(x, fileName = NULL, ...){
    if(is.null(x$plotData)){
        res <- data.frame()
    }else{
        res <- data.frame(x$plotData)
    }
   
    res$cluster <- if(is.null(x$cluster)) NA else x$cluster
    res$rho <- if(is.null(x$rho)) NA else x$rho
    res$delta <- if(is.null(x$delta)) NA else x$delta
    res$dc <- if(is.null(x$dc)) NA else x$dc
    res$halo <- if(is.null(x$halo)) NA else x$halo
    res$higherID <- if(is.null(x$higherID)) NA else x$higherID
    res$peakCheck <- seq_len(length(x$cluster)) %in% x$peakID 
     
    if(is.null(fileName))
      fileName <- "densityClusterX results.txt"
    
    write.table(res, file = fileName, col.names = TRUE, row.names = TRUE)
}

readDCXRes <- function(x){
    data <- read.table(x, header = TRUE, row.names = 1)
  
    res <- list()
    res$cluster <- data$cluster
    res$rho <- data$rho
    res$delta <- data$delta
    res$dc <- unique(data$dc)
    res$halo <- if(all(is.na(data$halo))) NULL else data$halo
    res$higherID <- data$higherID
    res$peakID <- seq_len(length(data$cluster))[data$peakCheck]
    res$plotData <- data[ ,setdiff(colnames(data), c("cluster", "rho", "delta", "dc", "halo", "peakCheck", "higherID"))]
  
    class(res) <- 'densityClusterX'
    return(res)
}



clusterPlot <- function (x, deleHalos=FALSE, addClusterLabel=TRUE, targetCenter=FALSE, point_size, ...) {
    UseMethod("clusterPlot", x)
}

clusterPlot.densityClusterX <- function(x, deleHalos = FALSE, addClusterLabel=TRUE, targetCenter=FALSE, point_size = 0.8, ...) {
    if(is.null(x$plotData))
        stop("Data can's be visualized due to high dimensionality, 
             Please try heatmapPlot to visualize the results! \n")
    if(ncol(x$plotData) == 2){
        df <- as.data.frame(x$plotData)
    }else if(ncol(x$plotData) > 2){
        warning("Plot data has more then 2 dimensions, only the first two columns are selected for plot!")
        df <- as.data.frame(x$plotData)[ ,c(1,2)]
    }else{
        stop("Wrong plot data!")
    }
    
    xvar <- colnames(df)[1]
    yvar <- colnames(df)[2]
    peakDf <- df[x$peakID, ]
    clusterText <- x$cluster[x$peakID]
    
    if(!(is.null(x$halo)) && deleHalos == FALSE){
        coreDf <- df[!(x$halo), ]
        haloDf <- df[x$halo, ]
        coreDf$colour <- as.factor(x$cluster[!(x$halo)])
        haloDf$colour <- as.factor(x$cluster)[x$halo]
        p <- ggplot(coreDf, aes_string(x=xvar, y=yvar)) + 
            geom_point(aes(colour = colour), size = point_size) + 
            geom_point(data=haloDf, aes(colour=colour), alpha = 0.5, size = point_size) 
    }else if(!(is.null(x$halo)) && deleHalos == TRUE){
        coreDf <- df[!(x$halo), ]
        coreDf$colour <- as.factor(x$cluster[!(x$halo)])
        p <- ggplot(coreDf, aes_string(x=xvar, y=yvar)) + 
            geom_point(aes(colour = colour), size = point_size)
    }else{
        df$colour <- as.factor(x$cluster)
        p <- ggplot(df, aes_string(x=xvar, y=yvar)) + 
            geom_point(aes(colour = colour), size = point_size)
    }
    
    P <- guides(colour = guide_legend(title = "Cluster", override.aes = list(size = 4))) + theme_bw()
    if(addClusterLabel == TRUE)
        p <- p + annotate("text", label = clusterText, x=peakDf[,1], y = peakDf[,2], size = 8, colour = "black")
    if(targetCenter == TRUE)
        p <- p + geom_point(data = peakDf, shape = 13, size = point_size * 2)
    
    p
}



densityPlot <- function (x, point_size = 0.8, addClusterLabel=FALSE, targetCenter=TRUE, ...) {
    UseMethod("densityPlot", x)
}

densityPlot.densityClusterX <- function(x, point_size = 0.8, addClusterLabel=FALSE, targetCenter=TRUE, ...) {
    if(is.null(x$plotData))
        stop("Data can's be visualized due to high dimensionality, 
             Please try heatmapPlot to visualize the cluster results! \n")
    if(ncol(x$plotData) == 2){
        df <- as.data.frame(x$plotData)
    }else if(ncol(x$plotData) > 2){
        warning("Plot data has more then 2 dimensions, only the first two columns are selected for plot!")
        df <- as.data.frame(x$plotData)[ ,c(1,2)]
    }else{
        stop("Wrong plot data!")
    }
    
    df$colour <- x$rho
    xvar <- colnames(df)[1]
    yvar <- colnames(df)[2]
    peakDf <- df[x$peakID, ]
    clusterText <- x$cluster[x$peakID]
    p <- ggplot(df, aes_string(x=xvar, y=yvar)) + 
        geom_point(aes(colour = colour), size = point_size) + theme_bw() + 
        scale_color_gradient2(low = "blue", high = "red", midpoint = median(x$rho))
    
    if(addClusterLabel == TRUE)
        p <- p + annotate("text", label = clusterText, x=peakDf[,1], y = peakDf[,2], size = 8, colour = "black")
    if(targetCenter == TRUE)
        p <- p + geom_point(data = peakDf, shape = 10, size = point_size * 4)
    
    p
}



peakPlot <- function (x, addClusterLabel=FALSE, ...) {
    UseMethod("peakPlot", x)
}
peakPlot.densityClusterX <- function(x, addClusterLabel=FALSE, ...) {
    df <- data.frame(rho = x$rho, delta = x$delta)
    peakDf <- df[x$peakID, ]
    clusterText <- x$cluster[x$peakID]
    peakDf$colour <- factor(1:length(x$peakID))
    p <- ggplot(df, aes(x=rho, y=delta)) + geom_point(size = 2) +
        geom_point(data = peakDf, aes(colour=colour), size = 4) + 
        guides(colour = guide_legend(title = "Cluster", override.aes = list(size = 4))) +
        theme_bw()
    
    if(addClusterLabel == TRUE)
        p <- p + annotate("text", label = clusterText, x=peakDf[,1], y = peakDf[,2], size = 6, colour = "black")
    
    p
}



#' Dimensionality reduction for high-dimensional data
#' 
#' When the dimensionality of the data matrix is high, a dimensionality 
#' reduction is suggested. A linear method PCA and a non-linear method 
#' t-SNE are provided(suggested). 
#' 
#' @param data Numeric matrix of data or data frame.
#' @param method Dimensionality reduction method, Options are: 'tsne' | 'pca'. default is 'tsne'.
#' @param outDim the output number of dimensionalty, default value is 2
#' 
#' @return a matrix with columns of outDim, row names matching the row names of data
dimReduction <- function(data, method = "tsne", outDim = 2) {
    
    cat(paste0("Dimension Transformation using ", method, " starts...  \n"))
    mapped <- NULL
    if (method == "pca") {
        pca <- prcomp(data, scale = TRUE)
        mapped <- pca$x[,1:outDim]
        colnames(mapped) <- paste("PCA_dim", c(1:ncol(mapped)), sep = "")
        rownames(mapped) <- row.names(data)
    }else if (method == "tsne") {
        tsne_out <- Rtsne(as.matrix(data), initial_dims = dim(as.matrix(data))[2],
                          dims = outDim, perplexity = 30, theta = 0.1,
                          check_duplicates = FALSE, pca = TRUE)
        
        mapped <- tsne_out$Y
        colnames(mapped) <- paste("tsne", c(1:ncol(mapped)), sep="_")
        rownames(mapped) <- row.names(data)
    }
    cat(paste0("Successfully transform the data to ", outDim, " dimensions using ", method,"!\n" ))
    return(mapped)
}



#' Estimate the distance cutoff (density neighbourhood) from down-sampled data
#' 
#' This function estimate a distance cutoff value from the down-samples data,
#' wchich meet the criteria that the average neighbor rate (number of points 
#' within the distance cutoff value) fall between the provided range. 
#' 
#' @param data Numeric matrix of data or data frame.
#' @param sampleSize The size of the down-sampled data.
#' @param neighborRateLow The lower bound of the neighbor rate (default 0.01).
#' @param neighborRateHigh The upper bound of the neighbor rate (default 0.15).
#' 
#' @return A numeric value giving the estimated distance cutoff value
#' 
#' @examples
#' data <- iris[,1:4]
#' estimateDc(data)
#' 
#' @export
#' 
estimateDc <- function(data, sampleSize = 10000, distMethod = "euclidean", 
                       neighborRateLow=0.01, neighborRateHigh=0.02) {
    data <- as.matrix(data)
    dataLens <- nrow(data)
    if(dataLens > sampleSize){
        sample <- data[sample(1:dataLens, sampleSize, replace = FALSE), ]
    }else{ sample <- data }
    
    if(distMethod == "euclidean"){
        comb <- as.matrix(dist(sample, method = "euclidean"))
    }else{
        suppressWarnings(comb <- ppdist(sample, sample))
    }
    size <- nrow(comb)
    dc <- min(comb)
    dcMod <- median(comb)*0.05
    
    while(TRUE) {
        neighborRate <- mean((apply(comb < dc, 1, sum)-1)/size)
        #neighborRate <- mean(apply((exp(-(comb/dc)^2)), 1, sum) - 1) / size
        if(neighborRate > neighborRateLow && neighborRate < neighborRateHigh) break  
        if(neighborRate >= neighborRateHigh) {
            dc <- dc - dcMod
            dcMod <- dcMod/2
        }else{
            dc <- dc + dcMod
        }
    }
    cat('Distance cutoff calculated to', dc, '\n')
    dc
}


#' Computes the local density of points in a data matrix
#' 
#' This function calculate the local density for each point in the matrix. 
#' With a rowise implementation of the pairwise distance calculation, makes 
#' the local density estimation faster and memory efficient. A big benifit 
#' is the aviliability for big data. Parallel computing is supported for 
#' fast calculation. The computation can either be done using a simple summation 
#' of the points with the distance cutoff for each observation, or by applying 
#' a gaussian kernel scaled by the distance cutoff (more robust for low-density data)
#' 
#' @param data Numeric matrix of data or data frame.
#' @param distMethod the method used for distance calculation
#' @param dc A numeric value specifying the distance cutoff.
#' @param gaussian Logical. Should a gaussian kernel be used to estimate the density (defaults to TRUE).
#' @param ifParallel A boolean decides if run parallelly
#' 
#' @return A vector of local density values with index matching the row names of data.
#' 
localDensity <- function(data, distMethod, dc, gaussian=FALSE, ifParallel = FALSE) {
    cat("calculate local density...")
    splitFactor <- splitFactorGenerator(nrow(data))
    dataFolds <- split.data.frame(data, splitFactor)

    rholist <- llply(dataFolds, function(datai, distMethod, data, dc, gaussian) {
        if(distMethod == "euclidean"){
            suppressWarnings(idist <- as.matrix(pdist(datai, data)))
        }else{
            suppressWarnings(idist <- ppdist(datai, data) ) }
        
        if(gaussian){
            apply((exp(-(idist/dc)^2)), 1, sum) - 1
        }else{
            apply(idist < dc, 1, sum) - 1 } 
        }, data = data, distMethod = distMethod, 
        dc = dc, gaussian = gaussian, .parallel = ifParallel)
    
    rho <- do.call(base::c, rholist)
    if(is.null(row.names(data))) {
        names(rho) <- NULL
    } else {
        names(rho) <- row.names(data)
    }
    cat("DONE!\n")
    rho
}


#' Calculate distance to closest observation of higher density
#' 
#' This function finds, for each observation, the minimum distance to an 
#' observation of higher local density. With a rowise implementation of 
#' the pairwise distance calculation, makes the local density estimation 
#' faster and memory efficient. A big benifit is the aviliability for big 
#' data. Parallel computing is supported for fast calculation.
#' 
#' @param data Numeric matrix of data or data frame.
#' @param rho A vector of local density values as outputted by \code{localDensity}
#' @param distMethod the method used for distance calculation
#' @param ifParallel A boolean decides if run parallelly
#' 
#' @return A list of distances to closest observation of higher density and the ID
#' 
minDistToHigher <- function(data, rho, distMethod, ifParallel = parallel) {
    cat("Search nearest neighobour with higher density...")
    splitFactor <- splitFactorGenerator(nrow(data))
    dataFolds <- split.data.frame(data, splitFactor)
    rhoFolds <- split(rho, splitFactor)
    
    dataRhoList <- mapply(function(datai, rhoi) {
        list(datai, rhoi)}, dataFolds, rhoFolds, SIMPLIFY = FALSE)
    
    deltaWidList <- llply(dataRhoList, function(x, data, distMethod, rho){
        datai <- x[[1]]
        rhoi <- x[[2]]
        
        if(distMethod == "euclidean"){
            suppressWarnings(datai2dataDist <- as.matrix(pdist(datai, data)))
        }else{
            suppressWarnings(datai2dataDist <- ppdist(datai, data)) }
        # datai2dataDist <- as.matrix(pdist::pdist(datai, data, rho))
        
        rhoi2rhoComp <- do.call(rbind, lapply(rhoi, function(x) x < rho)) ## <=
        sapply(seq_len(nrow(datai2dataDist)), function(i){
            distToAllHigherPoints <- datai2dataDist[i, rhoi2rhoComp[i, ]]
            if(length(distToAllHigherPoints) == 0) {
                c(max(datai2dataDist[i, ]), which.max(datai2dataDist[i, ]))
            } else {
                c(min(distToAllHigherPoints), 
                  which(rhoi2rhoComp[i, ] == TRUE)[which.min(distToAllHigherPoints)] )
            }} )
    }, data = data, distMethod = distMethod, rho = rho, .parallel = ifParallel)
    
    deltaWid <- do.call(cbind, deltaWidList)
    delta <- deltaWid[1, ]
    names(delta) <- names(rho)
    id <- deltaWid[2, ]
    names(id) <- names(rho)
    cat("DONE!\n")
    return(list(delta = delta, higherID = id))
}



#' Automatic peak detection
#' 
#' Automatic detect peaks by searching high denisty point with anomalous large distance to
#' higher denisty peaks. rho and delta are transformed to one index, and the anomalous peaks
#' are detected using generalized ESD method.
#' 
#' @param rho A vector of the local density, outout of \code{localDensity}
#' @param delta A vector of distance to closest observation of higher density
#' @param alpha The level of statistical significance for peak detection.
#' 
#' @return a vector containing the indexes of peaks
peakDetect <- function(rho, delta, alpha = 0.001){
    cat("Peak detection...")
    delta[is.infinite(delta)] <- max(delta[!(is.infinite(delta))])^2
    rdIndex <- scale01(rho) * delta   ## transform delta, important for big data
    #rdIndex <- log(rho + 1) * delta
    peakID1 <- detect_anoms_sd(rdIndex, direction = "pos", alpha = alpha)
    peakID2 <- detect_anoms_sd(delta, direction = "pos", alpha = alpha)
    peakID <- intersect(peakID1, peakID2)
    cat("DONE!\n")
    peakID
}

scale01 <- function(x){
    x <- (x-min(x))/(max(x)-min(x)) 
    x
}

truncTrans <- function(x, a, b=1){
    if(a < b)
        x[x<=a] <- b
    x
}

#' assign clusters to non-peak points
#' 
clusterAssign <- function(peakID, higherID, rho){
    cat("Cluster assigning...")
    runOrder <- order(rho, decreasing = TRUE)
    cluster <- rep(NA, length(rho))
    for(i in runOrder) {
        cluster[i] <- ifelse(i %in% peakID, match(i, peakID), cluster[higherID[i]] ) }
    cat("DONE!\n")
    cluster
}

## differentiate halo form cores
haloDetect <- function(data, rho, distMethod = "euclidean", cluster, peakID, dc){
    cat("diffenentiate halos from cores ...")
    clusterRhoThres <- sapply(1:length(peakID), function(clusteri){
        dataOfClusteri <- data[cluster == clusteri, ]
        otherData <- data[cluster != clusteri, ]
        rhoOfClusteri <- rho[cluster == clusteri]
        splitFactor <- splitFactorGenerator(nrow(dataOfClusteri), nrow(otherData))
        dataFolds <- split.data.frame(dataOfClusteri, splitFactor)
        haveColseNeighbourList <- lapply(dataFolds, function(x){
            if(distMethod == "euclidean"){
                suppressWarnings(dataFoldi2otherDataDist <- as.matrix(pdist(x, otherData)))
            }else{
                suppressWarnings(dataFoldi2otherDataDist <- ppdist(x, otherData) ) }
            #dataFoldi2otherDataDist <- as.matrix(pdist(x, otherData))
            apply(dataFoldi2otherDataDist < dc/2, 1, any)
        })
        checkRes <- do.call(base::c, haveColseNeighbourList)
        rhoThres <- ifelse(any(checkRes), max(rhoOfClusteri[checkRes]), min(rhoOfClusteri))
        rhoThres
    })
    halo <- rho < clusterRhoThres[cluster]
    cat("DONE!\n")
    halo
}



## generate split factors to split rowNum into folds, each of size around foldSize
splitFactorGenerator <- function(rowNum, colNum){
    if(missing(colNum)){
        colNum <- rowNum
    }
    foldSize <- round(32772663 / colNum)  ## each chunk with maxi 250Mb
    foldNum <- ceiling(rowNum / foldSize)
    lastfoldSize <- rowNum - (foldNum-1) * foldSize
    if(foldNum > 1){
        splitFactor <- c(rep(1:(foldNum-1), each = foldSize), rep(foldNum, lastfoldSize))
    }else{
        splitFactor <- rep(foldNum, lastfoldSize)
    }
    return(splitFactor)
}




#' Outlier detection
#' 
#' Using generialized ESD to detect outliers, iterate and remove point with ares higher than lamda 
#' in a univariate data set assumed to come from a normally distributed population.
#' 
#' @param data A vectors of boservations.
#' @param max_anoms Maximal percentile of anomalies.
#' @param alpha The level of statistical significance with which to accept or reject anomalies.
#' @param direction Directionality of the anomalies to be detected. Options are: 'pos' | 'neg' | 'both'.
#' 
#' @return A vector containing indexes of the anomalies (outliers).
#' 
detect_anoms <- function(data, max_anoms=0.1, alpha = 0.01, direction='pos') {
    
    num_obs <- length(data)
    names(data) <- 1:num_obs
    data <- na.omit(data)
    max_outliers <- trunc(num_obs*max_anoms)
    
    anoms <- NULL
    for (i in 1L:max_outliers){
        ares <- switch(direction, 
                       pos = data - median(data),
                       neg = median(data) - data,
                       both = abs(data - median(data)) )
        
        p <- switch(direction, 
                    pos = 1 - alpha/(num_obs-i+1),
                    neg = 1 - alpha/(num_obs-i+1),
                    both = 1 - alpha/(2*(num_obs-i+1)) )
        
        data_sigma <- mad(data)
        if(data_sigma == 0) break
        ares <- ares/data_sigma
        maxAres <- max(ares)
        
        t <- qt(p, (num_obs-i-1))
        lam <- t*(num_obs-i) / sqrt((num_obs-i-1+t**2)*(num_obs-i+1))
        
        if(maxAres > lam){
            maxAres_id <- which(ares == maxAres)
            anoms <- c(anoms, names(maxAres_id))
            data <- data[-maxAres_id ]
        }else
            break
    }
    
    return(as.numeric(anoms))
}



detect_anoms_sd <- function(data, max_anoms=0.1, alpha = 0.01, direction='pos') {
    
    num_obs <- length(data)
    names(data) <- 1:num_obs
    data <- na.omit(data)
    max_outliers <- trunc(num_obs*max_anoms)
    
    anoms <- NULL
    for (i in 1L:max_outliers){
        ares <- switch(direction, 
                       pos = data - mean(data),
                       neg = mean(data) - data,
                       both = abs(data - mean(data)) )
        
        p <- switch(direction, 
                    pos = 1 - alpha/(num_obs-i+1),
                    neg = 1 - alpha/(num_obs-i+1),
                    both = 1 - alpha/(2*(num_obs-i+1)) )
        
        data_sigma <- sd(data)
        if(data_sigma == 0) break
        ares <- ares/data_sigma
        maxAres <- max(ares)
        
        t <- qt(p, (num_obs-i-1))
        lam <- t*(num_obs-i) / sqrt((num_obs-i-1+t**2)*(num_obs-i+1))
        
        if(maxAres > lam){
            maxAres_id <- which(ares == maxAres)
            anoms <- c(anoms, names(maxAres_id))
            data <- data[-maxAres_id ]
        }else
            break
    }
    
    return(as.numeric(anoms))
}

