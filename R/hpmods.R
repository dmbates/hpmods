##' Package to determine hierarchy-preserving models 
##' 
##' @name hpmods-package
##' @docType package
##' @author Douglas Bates \email{bates@@stat.wisc.edu}
##' @keywords package
NA

##' Determine hierarchy-preserving models
##'
##' Determine the hierarchy-preserving models from a formula
##' @title Hierarchy-preserving models
##' @param formula as in \code{"\link{lm}"}
##' @param data as in \code{"\link{lm}"}
##' @param subset as in \code{"\link{lm}"}
##' @param weights as in \code{"\link{lm}"}
##' @param na.action as in \code{"\link{lm}"}
##' @param offset as in \code{"\link{lm}"}
##' @param \dots optional, additional arguments.  At present, none are used.
##' @return an object of class \code{"hpmods"}
##' @author Douglas Bates
##' @keywords models
##' @examples
##' set.seed(12321)
##' fr <- data.frame(y = rnorm(900), x1 = rnorm(900), x2 = rnorm(900),
##'                  x3 = gl(3,10,900), x4 = gl(10,3,900))
##' (hpm <- hpmods(y ~ (x1 + x2 + x3 + x4)^4, fr))
##' @export
hpmods <- function(formula, data, subset, weights, na.action, offset, ...) {
    cl   <- match.call()
    mf   <- match.call(expand.dots = FALSE)
    m    <- match(c("formula", "data", "subset", "weights", "na.action", 
                    "offset"), names(mf), 0L)
    mf   <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf   <- eval(mf, parent.frame())
    atts <- attributes(terms(mf))
    inci <- crossprod(atts$factor) == atts$order # incidence matrix for terms
    mods <- array(FALSE, c(nrow(inci), 1))
    rownames(mods) <- rownames(inci)
    for (j in 1:ncol(inci))
        mods <- cbind(mods, t(unique(t(array(inci[,j], dim(mods)) | mods))))
    mods <- t(mods)
    rownames(mods) <- mods %*% 2^((seq_len(ncol(inci)))-1)
    res  <- list(call=cl, incidence=inci, models=mods,
                 frame=mf, X=model.matrix(terms(mf), mf))
    class(res) <- "hpmods"
    res
}

##' @S3method print hpmods
print.hpmods <- function(x, ...) {
    require("Matrix", quietly=TRUE, character.only=TRUE)
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(as(x$models, "sparseMatrix"))
    invisible(x)
}

##' Extract model subsets from an hpmods object
##'
##' After fitting a model with several terms the residual sums of
##' squares (and the coefficient estimates, if needed) for any leading
##' subset of those terms can be determined directly.  This function
##' extracts distinct subsets of the rank-preserving models that can
##' where each subset can be analyzed from a single model fit.
##' @title Extract model subsets
##' @param x an object of class \code{"hpmods"}
##' @param \dots optional, additional arguments.  At present, none are used.
##' @return a list of model subsets
##' @author Douglas Bates
##' ##' set.seed(12321)
##' fr <- data.frame(y = rnorm(900), x1 = rnorm(900), x2 = rnorm(900),
##'                  x3 = gl(3,10,900), x4 = gl(10,3,900))
##' msubs(hpmods(y ~ (x1 + x2 + x3 + x4)^4, fr))
##' @export
msubs <- function(x, ...) {
    stopifnot(class(x) == "hpmods")
    mm  <- x$models
    ans <- vector(mode="list", length=nrow(mm)-ncol(mm))
    i   <- 1L
    p2  <- 2^(0:(ncol(mm)-1))
    while (nrow(mm) > 0) {
        cr <- as.character(cumsum(p2[mm[nrow(mm),]]))
        if (i == 1L) cr <- c("0", cr)   # null model goes in first group
        ans[[i]] <- mm[intersect(rownames(mm), cr),,drop=FALSE]
        i <- i + 1L
        mm <- mm[setdiff(rownames(mm), cr),,drop=FALSE]
    }
    ans[!sapply(ans, is.null)]
}
