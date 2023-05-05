suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(SingleCellExperiment)
})

.df <- \(x, margin = 2, features = NULL, assay = "logcounts") {
    fun <- list(rowData, colData)[[margin]]
    df <- data.frame(fun(x))
    if (margin == 2) {
        if (is(x, "SpatialExperiment")) {
            xy <- spatialCoords(x)
            xy <- data.frame(xy)
            df <- cbind(df, xy)
        }
        if (length(reducedDims(x)) > 0) {
            for (. in reducedDimNames(x)) {
                y <- reducedDim(x, .)
                n <- seq(ncol(y))
                colnames(y) <- paste0(., n)
                reducedDim(x, .) <- y
            }
            dr <- do.call(cbind, reducedDims(x))
            df <- cbind(df, dr)
        }
        if (!is.null(features)) {
            fs <- intersect(features, rownames(x))
            y <- t(assay(x, assay)[fs, ])
            y <- data.frame(y, check.names = FALSE)
            df <- cbind(df, y)
        }
    }
    return(df)
}

.scale01 <- \(x, q = 0.01) {
    if (!is(x, "matrix")) 
        x <- as.matrix(x)
    if (any(dim(x)) == 1) {
        qs <- quantile(x, c(q, 1 - q))
        x <- (x - qs[1])/(qs[2] - qs[1])
    } else {
        qs <- c(rowQuantiles, colQuantiles)[[margin]]
        qs <- qs(x, probs = c(q, 1 - q))
        qs <- matrix(qs, ncol = 2)
        x <- switch(margin, 
            `1` = (x - qs[, 1])/(qs[, 2] - qs[, 1]), 
            `2` = t((t(x) - qs[, 1])/(qs[, 2] - qs[, 1])))
    }
    x[x < 0 | is.na(x)] <- 0
    x[x > 1] <- 1
    return(x)
}
