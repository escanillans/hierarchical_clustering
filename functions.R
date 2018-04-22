# This function runs bottom-up hierarchical clustering
# until a specified number, k, of clusters remain.
# Input:
# dat = data
# file = data that contains the identifier, the name, and 
# description of each gene.
# linkType = link type (single, complete, or average)
# k = number of clusters
# Output:
# output a list of genes for each k cluster
hierarchicalClutersting <- function(dat, file, linkType, k)
{
  # keep track of genes in each cluster k
  result = vector("list", k)
  
  # each instance is in its own cluster
  clusterList = as.list(1:nrow(dat))
  
  # compute initial distance matrix
  currDistMatrix = distMatrix(dat)
  
  while(length(clusterList) > k)
  {
    # find two closest clusters based on linkType
    if(linkType == "S")
    {
      # find two closest clusters in clusterList
      singleLinkResult = singleLink(currDistMatrix)
      clustersIdx = singleLinkResult$ind
      
      # add new cluster (merged) to clusterList
      mergedCluster = c()
      for(i in clustersIdx)
      {
        # get ith clusterIdx list of elements and attach to mergedCluster
        mergedCluster = c(mergedCluster, clusterList[[i]])
      }
      clusterList[[length(clusterList)+1]] = mergedCluster
      
      # remove the two closest clusters from clusterList
      clusterList = clusterList[-clustersIdx] 
      
      # update distance matrix with new clusters
      currDistMatrix = updateMatrix_SingleLink(currDistMatrix, clustersIdx)
      
    } else if(linkType == "C")
    {
      # find two closest clusters in clusterList
      completeLinkResult = completeLink(currDistMatrix)
      clustersIdx = completeLinkResult$ind
      
      # add new cluster (merged) to clusterList
      mergedCluster = c()
      for(i in clustersIdx)
      {
        # get ith clusterIdx list of elements and attach to mergedCluster
        mergedCluster = c(mergedCluster, clusterList[[i]])
      }
      clusterList[[length(clusterList)+1]] = mergedCluster
      
      # remove the two closest clusters from clusterList
      clusterList = clusterList[-clustersIdx] 
      
      # update distance matrix with new clusters
      currDistMatrix = updateMatrix_CompleteLink(currDistMatrix, clustersIdx)
    } else if(linkType == "A")
    {
      # find two closest clusters in clusterList
      avgLinkResult = averageLink(currDistMatrix)
      clustersIdx = avgLinkResult$ind
      
      # add new cluster (merged) to clusterList
      mergedCluster = c()
      for(i in clustersIdx)
      {
        # get ith clusterIdx list of elements and attach to mergedCluster
        mergedCluster = c(mergedCluster, clusterList[[i]])
      }
      clusterList[[length(clusterList)+1]] = mergedCluster
      
      # remove the two closest clusters from clusterList
      clusterList = clusterList[-clustersIdx] 
      
      # update distance matrix with new clusters
      currDistMatrix = updateMatrix_AverageLink(currDistMatrix, clustersIdx)
    }else{
      cat("There is an error in your code.", "\n")
    }
  }
  
  return(list(ListOfClusters = clusterList, Matrix = currDistMatrix))
}

# This function calculates the euclidian distance of
# two vectors
euclidianDistance <- function(x1, x2)
{
  return(sqrt(sum((x1 - x2)^2)))
}

# This function creates a distance matrix 
distMatrix <- function(dat)
{
  finalMatrix = matrix(rep(0, nrow(dat)*nrow(dat)), nrow = nrow(dat), ncol = nrow(dat))
  
  # For each row in dat, compute distance between all clusters
  for(i in 1:nrow(finalMatrix))
  {
    for(j in 1:nrow(finalMatrix))
    {
      # compute euclidian distance between i and j
      #distances = euclidianDistance(dat[i,], dat[j,])
      #distances = sqrt(sum((dat[i,] - dat[j,])^2))
      distances = norm(dat[i,] - dat[j,], type = "2")
      finalMatrix[i,j] = distances
    }
  }
  
  return(finalMatrix)
}

# This function returns the closest clusters based on single link
singleLink <- function(distanceMatrix)
{
  # set diagonal to Inf
  singleLinkMatrix = distanceMatrix
  diag(singleLinkMatrix) = Inf
  
  # find min value
  id = which.min(singleLinkMatrix)
  
  # index will return (row,col) location of min value
  ind <- as.numeric(arrayInd(id, dim(singleLinkMatrix)))
  
  return(list(min = singleLinkMatrix[id], ind = ind))
}

updateMatrix_SingleLink <- function(distanceMatrix, clusters)
{
  # remove clusters from distanceMatrix
  cluster1_vals = distanceMatrix[clusters[1],-1*clusters]
  cluster2_vals = distanceMatrix[clusters[2],-1*clusters]
  minClusterVals = pmin(cluster1_vals, cluster2_vals)
  
  # Create a new distance matrix and add minClusterVals to end of matrix
  newDistanceMatrix = distanceMatrix[-1*clusters,-1*clusters]
  newDistanceMatrix = rbind(cbind(newDistanceMatrix, minClusterVals, 
                                  deparse.level = 0), c(minClusterVals, 0), 
                            deparse.level = 0)
  
  return(newDistanceMatrix)
}

# This function returns the closest clusters based on complete link
completeLink <- function(distanceMatrix)
{
  # set diagonal to Inf
  completeLinkMatrix = distanceMatrix
  diag(completeLinkMatrix) = Inf
  
  # find max value
  id = which.min(completeLinkMatrix)
  
  # index will return (row,col) location of min value
  ind <- as.numeric(arrayInd(id, dim(completeLinkMatrix)))
  
  return(list(min = completeLinkMatrix[id], ind = ind))
}

updateMatrix_CompleteLink <- function(distanceMatrix, clusters)
{
  # remove clusters from distanceMatrix
  cluster1_vals = distanceMatrix[clusters[1],-1*clusters]
  cluster2_vals = distanceMatrix[clusters[2],-1*clusters]
  maxClusterVals = pmax(cluster1_vals, cluster2_vals)
  
  # Create a new distance matrix and add maxClusterVals to end of matrix
  newDistanceMatrix = distanceMatrix[-1*clusters,-1*clusters]
  newDistanceMatrix = rbind(cbind(newDistanceMatrix, maxClusterVals, 
                                  deparse.level = 0), c(maxClusterVals, 0), 
                            deparse.level = 0)
  
  return(newDistanceMatrix)
}

# This function returns the closest clusters based on single link
averageLink <- function(distanceMatrix)
{
  # set diagonal to Inf
  avgLinkMatrix = distanceMatrix
  diag(avgLinkMatrix) = Inf
  
  # find min value
  id = which.min(avgLinkMatrix)
  
  # index will return (row,col) location of min value
  ind <- as.numeric(arrayInd(id, dim(avgLinkMatrix)))
  
  return(list(min = avgLinkMatrix[id], ind = ind))
}

updateMatrix_AverageLink <- function(distanceMatrix, clusters)
{
  # remove clusters from distanceMatrix
  cluster1_vals = distanceMatrix[clusters[1],-1*clusters]
  oneNorm_Cluster1Vals = norm(as.matrix(cluster1_vals), type = "o")
  
  cluster2_vals = distanceMatrix[clusters[2],-1*clusters]
  oneNorm_Cluster2Vals = norm(as.matrix(cluster2_vals), type = "o")
  
  avgClusterVals =(oneNorm_Cluster1Vals * cluster1_vals + oneNorm_Cluster2Vals * cluster2_vals)/(oneNorm_Cluster1Vals + oneNorm_Cluster2Vals)
  
  # Create a new distance matrix and add avgClusterVals to end of matrix
  newDistanceMatrix = distanceMatrix[-1*clusters,-1*clusters]
  newDistanceMatrix = rbind(cbind(newDistanceMatrix, avgClusterVals, 
                                  deparse.level = 0), c(avgClusterVals, 0), 
                            deparse.level = 0)
  
  return(newDistanceMatrix)
}

# This function prints a description of the clusters (in ascending order)
printClusterInfo <- function(clusterList, file)
{
  # First compute average for each cluster
  avgVals = c()
  
  # For each cluster, calculate the average val of the measurements for
  # the genes contained in the cluster
  for(i in 1:length(clusterList))
  {
    # compute avg
    avg = mean(rowMeans(file[clusterList[[i]],3:ncol(file)])[order(rowMeans(file[clusterList[[i]],3:ncol(file)]))])
    avgVals = c(avgVals, avg)  
  }
  
  # arrange each cluster according to avgVals
  ascendingOrderIdx = order(avgVals)
  
  # print each cluster by printing each gene in each cluster in order (ascending)
  for(i in ascendingOrderIdx)
  {
    for(geneIdx in clusterList[[i]][sort(rowMeans(file[clusterList[[i]],3:ncol(file)]), index.return = TRUE)$ix])
    {
      # get identifier for current cluster gene
      identifier = file[geneIdx,1]
      
      # get name and description for current cluster gene
      description = file[geneIdx,2]
      
      # get avg value of the measurements for genes contained in that cluster
      val = as.numeric(rowMeans(file[geneIdx, 3:ncol(file)]))
      
      # print out above information separated by tabs
      info = paste(identifier, description, sprintf("%.3f", val), sep = "\t")
      cat(info, "\n")
      
    }
    
    # print avg value of all measurements in cluster i
    cat(sprintf("%.3f", avgVals[i]), "\n")
    cat("\n")
    
  }
}










