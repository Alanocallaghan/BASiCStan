#' Thin wrapper for BASiCS models in \linkS4class{stanfit} form before conversion.
#' 
#' @slot FeatureNames,ObservationNames The feature (gene) and observation (cell) names of the data used to fit the model.
setClass(
    "BASiCStan",
    representation = list(
        FeatureNames = "character", 
        ObservationNames = "character",
        SizeFactors = "numeric"
    ),
    contains = "stanfit",
    validity = function(object) {
        if (length(featureNames(object)) != nfeats(object)) {
            stop("FeatureNames vector isn't the right length")
        }
        if (length(sampleNames(object)) != ncells(object)) {
            stop("CellNames vector isn't the right length")
        }
        if (!all(sizeFactors(object) > 0)) {
            stop("SizeFactors should be > 0")
        }
        if (length(sizeFactors(object)) != ncells(object)) {
            stop("length(SizeFactors) should match the number of cells")
        }
        TRUE
    }
)
fnames <- function(object) object@FeatureNames
cnames <- function(object) object@ObservationNames

setfnames <- function(object, value) {
    object@FeatureNames <- value
    object
}

setcnames <- function(object, value) {
    object@ObservationNames <- value
    object
}

nfeats <- function(x) {
    nom <- names(x)
    length(grep("delta", nom))
}

ncells <- function(x) {
    nom <- names(x)
    if (any(grepl("nu", nom))) {
        length(grep("nu", nom))
    } else {
        length(sizeFactors(x))
    }
}

#' Feature and sample names for \linkS4class{BASiCStan} objects.
#' 
#' Getters and setters for names of \linkS4class{BASiCStan} objects.
#' @param object An object of class \linkS4class{BASiCStan}.
#' @param value The value for feature or sample names to be set to.
#' @return
#' For the getters, a character vector.
#' For the setters, an object of class \linkS4class{BASiCStan}.
#' @importFrom Biobase featureNames sampleNames featureNames<- sampleNames<-
#' @importFrom BiocGenerics sizeFactors
#' @rdname sample-features-sizefactors
#' @export
setMethod("featureNames", "BASiCStan", fnames)
#' @rdname sample-features-sizefactors
#' @export
setMethod("sampleNames", "BASiCStan", cnames)
#' @rdname sample-features-sizefactors
#' @export
setMethod("featureNames<-", signature("BASiCStan", "ANY"), setfnames)
#' @rdname sample-features-sizefactors
#' @export
setMethod("sampleNames<-", signature("BASiCStan", "ANY"), setcnames)
#' @rdname sample-features-sizefactors
#' @export
setMethod("sizeFactors", signature("BASiCStan"), function(object) object@SizeFactors)

# #' @export
# setMethod("nrow", "BASiCStan", nfeats)
# #' @export
# setMethod("ncol", "BASiCStan", ncells)
