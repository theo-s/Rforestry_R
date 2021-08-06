#' @include R_honestRFTree.R

setClass(Class = "linearTree",
         contains = "RFTree")

# ------------------------------------------------------------------------------

#' @title findSplittingPoint_linear
#' @description Finds the best splitting value and feature to split on
#' @aliases selectBestFeature, selectBestFeature-method
#' @param x feature data set
#' @param y dependent outcome
#' @param featureList 
#' @param sampleIndex a list containing two integer vecotrs 
#' `averagingSampleIndex` and `splittingSampleIndex` which specify which of the 
#' features is used for splitting and which is used for estimating the beta in 
#' the leaves
#' @param nodesize a list containing `splittingNodeSize` and `averagingNodeSize`
#' which specify the nodesizes for the splitting and the aceraging sets.
#' @param splitrule = "variance",
#' @param categoricalFeatureCols = list()
#' @return A list of two outputs: "splitFeature" is the best feature to split
#' in order to minimize the split loss, "splitValue" is its corresponding split
#' value.
findSplittingPoint_linear <- function(x,
                                      y,
                                      featureList,
                                      sampleIndex = list(
                                        "averagingSampleIndex" = 1:length(y),
                                        "splittingSampleIndex" = 1:length(y)
                                      ),
                                      nodesize = list("splittingNodeSize" = 5,
                                                      "averagingNodeSize" = 5),
                                      splitrule = "variance",
                                      categoricalFeatureCols = list()) {
  # ----------------------------------------------------------------------------
  stop("To be implemented")
  # ----------------------------------------------------------------------------
  # Return the best splitFeature and splitValue
  return(
    list(
      "bestSplitFeature" = bestSplitFeatureAll[bestFeatureIndex],
      "bestSplitValue" = bestSplitValueAll[bestFeatureIndex]
    )
  )
}



# ------------------------------------------------------------------------------
#' @title honstRFTree Constructor
#' @name linearTree-constructor
#' @rdname linearTree-class
#' @description Create a linearTree by making specifc observatios as splitting
#' and averaging dataset.
#' @param x A data frame of all training predictors.
#' @param y A vector of all training responses.
#' @param mtry The number of variables randomly selected at each split point.
#' The default value is set to be one third of total number of features of the
#' training data.
#' @param nodesize The minimum observations contained in terminal nodes. The
#' default value is 5.
#' @param sampleIndex A list of the index of observations that are used as
#' averaging dataset. The index are based on the original dataset `x` and `y`
#' from forest. Essentially, `x[sampleIndex]` generates the whole splitting
#' dataset.
#' @param splitrule A string to specify how to find the best split among all
#' candidate feature values. The current version only supports `variance` which
#' minimizes the overall MSE after splitting. The default value is `variance`.
#' @param categoricalFeatureCols A list of index for all categorical data. Used
#' for trees to detect categorical columns.
#' @export linearTree
setGeneric(
  name = "linearTree",
  def = function(x,
                 y,
                 mtry,
                 nodesize,
                 sampleIndex,
                 splitrule,
                 categoricalFeatureCols) {
    standardGeneric("linearTree")
  }
)

linearTree <- function(x,
                       y,
                       mtry = max(floor(ncol(x) / 3), 1),
                       nodesize = list("averagingNodeSize" = 5,
                                       "splittingNodeSize" = 5),
                       sampleIndex = list(
                         "averagingSampleIndex" = 1:length(y),
                         "splittingSampleIndex" = 1:length(y)
                       ),
                       splitrule = "variance",
                       categoricalFeatureCols = list()) {
  # ----------------------------------------------------------------------------
  # Catch Errors
  honestRFTree_input_checker(nodesize, sampleIndex)
  
  # ----------------------------------------------------------------------------
  tree <- new("linearTree",
              sampleIndex = sampleIndex,
              root = list())
  
  # Grow the tree.
  root <- recursivePartition(
    x = x,
    y = y,
    mtry = mtry,
    nodesize = nodesize,
    sampleIndex = sampleIndex,
    splitrule = splitrule,
    categoricalFeatureCols = categoricalFeatureCols
  )
  
  tree@root <- list("node" = root)
  return(tree)
}