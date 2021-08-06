################################
### selectBestFeatureLinear Method ###
################################
#' selectBestFeatureLinear
#' @name selectBestFeatureLinear
#' @rdname selectBestFeatureLinear-RFTree
#' @description Find the best `splitfeature`, `splitValue` pair where
#' `splitfeature` is one of the features specified by `featureList`. 
#' @param x A data frame of all training predictors.
#' @param y A vector of all training responses.
#' @param featureList_spit A list of candidate variables to split on
#' @param featureList_lin A list of variables which should be used to fit the 
#' linear model on.
#' @param sampleIndex A list of index of dataset used in this node and its
#' children. `sampleIndex` contains two keys `averagingSampleIndex`
#' and `splittingSampleIndex`. `averagingSampleIndex` is used to generate
#' aggregated prediction for the node. `splittingSampleIndex` is used for
#' `honestRF` which stores the splitting data when creating the tree. In
#' default, `splittingSampleIndex` is the same as `averagingSampleIndex`.
#' @param nodesize The minimum observations contained in terminal nodes. This
#' parameter is actually a list containing the values for both
#' `splittingNodeSize` and `averagingNodeSize`.
#' @param categoricalFeatureCols A list of index for all categorical data. Used
#' for trees to detect categorical columns.
#' @export selectBestFeatureLinear
setGeneric(
  name = "selectBestFeatureLinear",
  def = function(x,
                 y,
                 featureList_spit,
                 featureList_lin,
                 lambda,
                 sampleIndex,
                 nodesize,
                 categoricalFeatureCols) {
    standardGeneric("selectBestFeatureLinear")
  }
)



#' @rdname selectBestFeatureLinear-RFTree
#' @aliases selectBestFeatureLinear, selectBestFeatureLinear-method
#' @return A list of two outputs: "splitFeature" is the best feature to split
#' in order to minimize the split loss, "splitValue" is its corresponding split
#' value.
#' @note This function is currently depreciated. It has been replaced by the
#' C++ version in the package. Although the functionality and parameters are
#' exactly the same.

selectBestFeatureLinear <- function(x,
                                    y,
                                    featureList_spit,
                                    featureList_lin,
                                    lambda,
                                    sampleIndex = list(
                                      "averagingSampleIndex" = 1:length(y),
                                      "splittingSampleIndex" = 1:length(y)
                                    ),
                                    nodesize = list("splittingNodeSize" = 5,
                                                    "averagingNodeSize" = 5),
                                    categoricalFeatureCols = list()) {
  
  if (!all(sampleIndex$averagingSampleIndex == sampleIndex$splittingSampleIndex)) {
    stop("honesty is not implemented yet")
  }
  
  # Get the number of total features
  mtry <- length(featureList_spit)
  
  # Initialize the minimum loss for each feature
  bestSplitLossAll <- rep(Inf, mtry)
  bestSplitValueAll <- rep(NA, mtry)
  bestSplitFeatureAll <- rep(NA, mtry)
  bestSplitCountAll <- rep(0, mtry)
  
  # Iterate each selected features
  for (i in 1:mtry) {
    currentFeature <- featureList_spit[i]
    allUniqueValues <-
      sort(as.numeric(union(unique(x[sampleIndex$splittingSampleIndex, currentFeature]),
                            unique(x[sampleIndex$averagingSampleIndex, currentFeature]))))
    # Test if the all the values for the feature are the same, then proceed
    if (length(allUniqueValues) == 1)
      next()
    
    # Test if the current feature is categorical
    if (currentFeature %in% categoricalFeatureCols) {
      # -- Categorical Feature -------------------------------------------------
      stop("Splitting on Categorical Features is not yet implemented")
    } else{
      # -- Continuous Feature --------------------------------------------------
      
      warning("Find_ridge_best_split does not take the min node size into account")
      best_split <- Find_ridge_best_split(
        feat = x,
        y = y,
        linear.idx = featureList_lin,
        current.splitting.idx = i,
        lambda = lambda,
        nodesize = nodesize
      )
      
      
      bestSplitLossAll[i] <- best_split["shifted_RSS_best"]
      bestSplitValueAll[i] <- best_split["split_val_best"]
      bestSplitFeatureAll[i] <- currentFeature
      bestSplitCountAll[i] <- NA
      # This used to contained, how often an optimal split happened, but it is
      # not implemented nay more.
    }
  }
  
  # Get the best split values among all features
  bestSplitLoss <- min(bestSplitLossAll)
  bestFeatureIndex <- which(bestSplitLossAll == bestSplitLoss)
  # If we found a feasible splitting point
  if (bestSplitLoss < Inf) {
    # A split is possible
    return(list("bestSplitFeature" = bestSplitFeatureAll[bestFeatureIndex],
                "bestSplitValue" = bestSplitValueAll[bestFeatureIndex]))
  } else{
    # If none of the features are possible, return NA
    return(list("bestSplitFeature" = NA,
                "bestSplitValue" = NA))
  }
}
