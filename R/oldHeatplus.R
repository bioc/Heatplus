## oldHeatplus.R
##
## R code for the original version of the annotated heatmaps (columns only)
## Apart from some helper functions, this is historical code 
##
## 2011-08-11   Alexander.Ploner@ki.se



#' Display Data as Heatmap
#' 
#' This function displays an expression data matrix as a heatmap. It is based
#' on an old version of \code{heatmap} in the \code{stats} package, but offers
#' more flexibility (e.g. skipping dendrograms, skipping row/column labelling,
#' adding a legend).
#' 
#' This function is just about to be deprecated. Please use \code{regHeatmap}
#' for new projects.
#' 
#' With all parameters at their default, this gives the same result as a very
#' old version of \code{heatmap} that was the base for the modifications. All
#' parameters of the same name have the same function as in \code{heatmap},
#' though \code{add.expr}, which can be used for adding graphical elements
#' after the call to \code{image}, will probably not produce useful results.
#' Note also that row- and column labels are optional, i.e. if the
#' corresponding \code{dimname} of \code{x} is \code{NULL}, no labels are
#' displayed.
#' 
#' Note that setting \code{Rowv} or \code{Colv} to \code{NA} completely
#' suppresses re-ordering of rows or columns as well as the corresponding
#' dendrogram. Setting both to \code{NA} works basically like \code{image}
#' (though you can still add a legend).
#' 
#' Setting \code{trim} to a number between 0 and 1 uses equidistant classes
#' between the (\code{trim})- and (1-\code{trim})-quantile, and lumps the
#' values below and above this range into separate open-ended classes. If the
#' data comes from a heavy-tailed distribution, this can save the display from
#' putting too many values into to few classes.
#' 
#' @param x the numerical data matrix to be displayed.
#' @param Rowv either a dendrogram or a vector of reordering indexes for the
#' rows.
#' @param Colv either a dendrogram or a vector of reordering indexes for the
#' columns.
#' @param distfun function to compute the distances between rows and columns.
#' Defaults to \code{dist}.
#' @param hclustfun function used to cluster rows and columns. Defaults to
#' \code{hclust}.
#' @param add.expr Expression to be evaluated after the call to \code{image}.
#' See Details.
#' @param scale indicates whether values should be scaled by either by row,
#' column, or not at all. Defaults to \code{row}.
#' @param na.rm logical indicating whther to remove NAs.
#' @param do.dendro logical vector of length two, indicating (in this order)
#' whether to draw the row and column dendrograms.
#' @param legend integer between 1 and 4, indicating on which side of the plot
#' the legend should be drawn, as in \code{mtext}.
#' @param legfrac fraction of the plot that is taken up by the legend; larger
#' values correspond to smaller legends.
#' @param col the color scheme for \code{image}. The default sucks.
#' @param trim Percentage of values to be trimmed. This helps to keep an
#' informative color scale, see Details.
#' @param \dots extra arguments to \code{image}.
#' @return Same as \code{heatmap} with \code{keep.dendro=FALSE}: an invisible
#' list giving the reordered indices of the row- and column-elements as
#' elements \code{rowInd} and \code{colInd}.
#' @author Original by Andy Liaw, with revisions by Robert Gentleman and Martin
#' Maechler.
#' 
#' Alexander Ploner for this version.
#' @seealso \code{\link{heatmap}}, \code{\link{hclust}},
#' \code{\link{heatmap_plus}}, \code{\link{regHeatmap}},
#' \code{\link{annHeatmap}}
#' @keywords hplot
#' @examples
#' 
#' # create data
#' mm = matrix(rnorm(1000, m=1), 100,10)
#' mm = cbind(mm, matrix(rnorm(2000), 100, 20))
#' mm = cbind(mm, matrix(rnorm(1500, m=-1), 100, 15))
#' mm2 = matrix(rnorm(450), 30, 15)
#' mm2 = cbind(mm2, matrix(rnorm(900,m=1.5), 30,30))
#' mm=rbind(mm, mm2)
#' colnames(mm) = paste("Sample", 1:45)
#' rownames(mm) = paste("Gene", 1:130)
#' 
#' # similar to base heatmap
#' heatmap_2(mm)
#' 
#' # remove column dendrogram
#' heatmap_2(mm, do.dendro=c(TRUE, FALSE))
#' 
#' # add a legend under the plot
#' heatmap_2(mm, legend=1)
#' # make it smaller
#' heatmap_2(mm, legend=1, legfrac=10)
#' # ... on the left side
#' heatmap_2(mm, legend=2, legfrac=10)
#' 
#' # remove the column labels by removing the column names
#' colnames(mm)=NULL
#' heatmap_2(mm, legend=1, legfrac=10)
#' 
#' # truncate the data drastically
#' heatmap_2(mm, legend=1, legfrac=10, trim=0.1)
#' 
#' @export heatmap_2
heatmap_2 = function (x, Rowv, Colv, distfun = dist, hclustfun = hclust, 
                      add.expr, scale = c("row", "column", "none"), na.rm = TRUE, 
                      do.dendro=c(TRUE,TRUE), legend=0, legfrac=8, 
                      col=heat.colors(12), trim, ...) 
#
# Name: heatmap_2
# Desc: an extended version of the heatmap()-function in library(mva)
# Auth: ?heatmap for original authors
#       Modifications by Alexander.Ploner@meb.ki.se 271003 onwards
#
# Chng: 031120 AP added trim
#       040129 AP added na.rm to quantiles
#       060303 AP changed name 
#       060602 AP Rowv/Colv=NA suppresses re-ordering completely,
#                 margins are now wider if no dendro is plotted
{
    # scaling & clustering, unchanged
    scale <- match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    r.cex <- 0.2 + 1/log10(nr)
    c.cex <- 0.2 + 1/log10(nc)
    if (missing(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (missing(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    ## Set up dendrograms otherwise told not to
    if (!identical(Rowv, NA)) {
        if (!inherits(Rowv, "dendrogram")) {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
        } else ddr <- Rowv
        rowInd <- order.dendrogram(ddr)    
    } else {
        rowInd = 1:nr
        do.dendro[1] = FALSE
    }
    if (!identical(Colv, NA)) {        
        if (!inherits(Colv, "dendrogram")) {
            hcc <- hclustfun(distfun(t(x)))
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, Colv)
        } else ddc <- Colv
        colInd <- order.dendrogram(ddc)
    } else {
        colInd = 1:nc
        do.dendro[2] = FALSE
    }
    x <- x[rowInd, colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sd <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sd, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sd <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sd, "/")
    }
    ## Only margin played with here; using the comprehensive list from
    ## no.readonly=FALSE crashed R CMD check
    op <- par()["mar"]
    on.exit(par(op))
    # slightly dirty: trim of extremes to make the legend more reasonable
    if (!missing(trim)) {
        trim = min(trim[1], 1-trim[1])
        lo = quantile(x, trim, na.rm=na.rm) 
        hi = quantile(x, 1-trim, na.rm=na.rm)
        x[x<lo] = lo
        x[x>hi] = hi
    }
    # Modifications: set up the plot properly
    do.xaxis = !is.null(colnames(x))
    do.yaxis = !is.null(rownames(x))
    margin = rep(0, 4)
    margin[1] = if (do.xaxis) 5 else 2
    margin[2] = if (do.dendro[1]) 0  else 2
    margin[3] = if (do.dendro[2]) 0  else 2    
    margin[4] = if (do.yaxis) 5 else 2
    # Choose layout according to dendrograms
    if (do.dendro[1] & do.dendro[2]) {
        ll = matrix(c(0, 3, 2, 1), 2, 2, byrow = TRUE)
        ll.width  = c(1,4)
        ll.height = c(1, 4)
    } else if (do.dendro[1]) {
        ll = matrix(c(2, 1), 1, 2, byrow = TRUE)
        ll.width  = c(1,4)
        ll.height = 4
    } else if (do.dendro[2]) {
        ll = matrix(c(2, 1), 2, 1, byrow = FALSE)
        ll.width  = 4
        ll.height = c(1,4)
    } else {
        # it's okay to plot nothing!
        ll = matrix(1, 1, 1)
        ll.width  = 1
        ll.height = 1
    }
    # add the legend at the required side
    if (legend %in% 1:4) {
        plotnum = max(ll) + 1
        nc = ncol(ll)
        nr = nrow(ll)
        if (legend == 1) {
            ll = rbind(ll, if (nc==1) plotnum else c(0,plotnum))
            ll.height = c(ll.height, sum(ll.height)/(legfrac-1))
            leg.hor =TRUE
        } else if (legend==2) {
            ll = cbind(if (nr==1) plotnum else c(0,plotnum), ll)
            ll.width = c(sum(ll.width)/(legfrac-1), ll.width)
            leg.hor = FALSE
        } else if (legend == 3) {
            ll = rbind(if (nc==1) plotnum else c(0,plotnum), ll)
            ll.height = c(sum(ll.height)/(legfrac-1), ll.height)
            leg.hor =TRUE            
        } else if (legend==4) {
            ll = cbind(ll, if (nr==1) plotnum else c(0,plotnum))
            ll.width = c(ll.width, sum(ll.width)/(legfrac-1))
            leg.hor =FALSE
        }
    }
    # now go
    layout(ll, width=ll.width, height=ll.height, respect=TRUE)            
    # the actual plot            
    par(mar=margin)
    image(1:ncol(x), 1:nrow(x), t(x), axes = FALSE, xlim = c(0.5, 
        ncol(x) + 0.5), ylim = c(0.5, nrow(x) + 0.5), xlab = "", 
        ylab = "", col=col, ...)
    # only do the axis if need be
    if (do.xaxis) {
        axis(1, 1:ncol(x), las = 2, line = -0.5, tick = 0, 
            labels = colnames(x), cex.axis = c.cex)
    }
    if (do.yaxis) {
        axis(4, 1:nrow(x), las = 2, line = -0.5, tick = 0, 
        labels = rownames(x), cex.axis = r.cex)
    }
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (do.dendro[1]) {
        mm = margin
        mm[2] = 3
        mm[4] = 0
        par(mar = mm)                
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    if (do.dendro[2]) {
        mm = margin
        mm[1] = 0
        mm[3] = 3
        par(mar = mm)            
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    # do the legend, if required
    # FIXME: the vertical legends do not work, the axis appears always below
    if (legend %in% 1:4) {    
        dummy.x <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length=length(col))
        dummy.z <- matrix(dummy.x,ncol=1)
        if (leg.hor) {
            par(mar=c(2, margin[2], 2, margin[4]))
            image(x=dummy.x, y=1, z=dummy.z, yaxt="n", col=col)
        } else {
            par(mar=c(margin[1], 2, margin[3], 2))            
            image(x=1, y=dummy.x, z=t(dummy.z), xaxt="n", col=col)
        }        
    }
    
    invisible(list(rowInd = rowInd, colInd = colInd))
}



#' Display an Annotated Heatmap
#' 
#' This function displays an expression data matrix as a heatmap with a column
#' dendrogram. A given clustering will be shown in color.  Additionally, a
#' number of binary and interval scaled covariates can be added to characterize
#' these clusters.
#' 
#' This function is just about to be deprecated. Please use functions
#' \code{annHeatmap} or \code{annHeatmap2} for new projects.
#' 
#' This is a heavily modified version of \code{heatmap_2}, which is a heavily
#' modfied version of an old version of \code{heatmap} in package \code{stats},
#' so some of the arguments are described in more detail there. The main
#' distinguishing feature of this routine is the possibility to color a cluster
#' solution, and to add a covariate display.
#' 
#' Covariates are assumed to be binary, coded as 0 and 1 (or \code{FALSE} and
#' \code{TRUE} respectively). One of the covariates can be interval scaled, the
#' column index of this variable is supplied via argument \code{covariate}. The
#' details of the added display are handled by the function \code{picketplot}.
#' 
#' Setting \code{trim} to a number between 0 and 1 uses equidistant classes
#' between the (\code{trim})- and (1-\code{trim})-quantile, and lumps the
#' values below and above this range into separate open-ended classes. If the
#' data comes from a heavy-tailed distribution, this can save the display from
#' putting too many values into to few classes. Alternatively, you can set
#' \code{equal=TRUE}, which uses an equidistant color scheme for the ranks of
#' the values.
#' 
#' @param x the numerical data matrix to be displayed.
#' @param addvar data frame with (mostly binary) covariates.
#' @param covariate integer indicating the one column in \code{addvar} that is
#' interval scaled.
#' @param picket.control list of option for drawing the covariates, passed to
#' \code{oldPicketplot}.
#' @param h height at which to cut the dendrogram, as in \code{oldCutree};
#' overrides \code{clus}.
#' @param clus an explicit vector of cluster memberships for the columns of
#' \code{x}, if no dendrogram is used; ignored if \code{do.dendro=TRUE} and
#' \code{h} is specified.
#' @param cluscol a vector of colors used to indicate clusters.
#' @param cluslabel labels to designate cluster names.
#' @param Rowv either a dendrogram or a vector of reordering indexes for the
#' rows.
#' @param Colv either a dendrogram or a vector of reordering indexes for the
#' columns.
#' @param reorder logical vector of length two, indicating whether the rows and
#' columns (in this order) should be reordered using \code{order.dendrogram}.
#' @param distfun function to compute the distances between rows and columns.
#' Defaults to \code{dist}.
#' @param hclustfun function used to cluster rows and columns. Defaults to
#' \code{hclust}.
#' @param scale indicates whether values should be scaled by either by row,
#' column, or not at all. Defaults to \code{row}.
#' @param na.rm logical indicating whther to remove NAs.
#' @param do.dendro logical indicating whether to draw the column dendrogram.
#' @param col the color scheme for \code{image}. The default sucks.
#' @param trim Percentage of values to be trimmed. This helps to keep an
#' informative color scale, see Details.
#' @param equalize logical indicating whther to use the ranks of the data for
#' setting the color scheme; alternative to \code{trim}, see Details.
#' @param \dots extra arguments to \code{image}.
#' @return A list with components \item{rowInd}{indices of the rows of the
#' display in terms of the rows of \code{x}.} \item{colInd}{ditto for the
#' columns of the display.} \item{clus}{the cluster indices of the columns of
#' the display.}
#' @author Original by Andy Liaw, with revisions by Robert Gentleman and Martin
#' Maechler.
#' 
#' Alexander Ploner for the modifications documented here.
#' @seealso \code{\link{heatmap_2}}, \code{\link{heatmap}},
#' \code{\link{oldPicketplot}}, \code{\link{oldCutplot.dendrogram}},
#' \code{\link{RGBColVec}}, \code{\link{annHeatmap}}, \code{\link{annHeatmap2}}
#' @keywords hplot
#' @examples
#' 
#' # create data
#' mm = matrix(rnorm(1000, m=1), 100,10)
#' mm = cbind(mm, matrix(rnorm(2000), 100, 20))
#' mm = cbind(mm, matrix(rnorm(1500, m=-1), 100, 15))
#' mm2 = matrix(rnorm(450), 30, 15)
#' mm2 = cbind(mm2, matrix(rnorm(900,m=1.5), 30,30))
#' mm=rbind(mm, mm2)
#' colnames(mm) = paste("Sample", 1:45)
#' rownames(mm) = paste("Gene", 1:130)
#' addvar = data.frame(Var1=rep(c(0,1,0),c(10,20,15)),
#'                     Var2=rep(c(1,0,0),c(10,20,15)),
#'                     Var3=rep(c(1,0), c(15,30)),
#'                     Var4=rep(seq(0,1,length=4), c(10,5,15,15))+rnorm(45, sd=0.5))
#' addvar[3,3] = addvar[17,2] = addvar[34,1] =NA
#' colnames(addvar) = c("Variable X","Variable Y", "ZZ","Interval")
#' 
#' 
#' # the lame default, without clustering
#' # Labels do not look too hot that way
#' heatmap_plus(mm)
#' 
#' # without labels, but with cluster
#' dimnames(mm)=NULL
#' heatmap_plus(mm, h=40)
#' 
#' # add some covariates, with nice names
#' heatmap_plus(mm, addvar=addvar, cov=4)
#' 
#' # covariates and clustering
#' heatmap_plus(mm, addvar=addvar, cov=4, h=20, col=RGBColVec(64), equal=TRUE)
#' 
#' # Clustering without the dendrogram
#' cc = cutree(hclust(dist(t(mm))), k=5)
#' heatmap_plus(mm, addvar=addvar, cov=4, clus=cc, do.dendro=FALSE)
#' 
#' @export heatmap_plus
heatmap_plus = function (x, addvar, covariate=NULL, picket.control=list(),
                         h, clus, cluscol, cluslabel=NULL, Rowv, Colv, 
                         reorder=c(TRUE, TRUE), distfun = dist, hclustfun = hclust, 
                         scale = c("row", "column", "none"), na.rm = TRUE, 
                         do.dendro=TRUE, col=heat.colors(12), trim, 
                         equalize=FALSE, ...) 
#
# Name: heatmap_plus
# Desc: a very much revised version of heatmap_2, which is 
#       a very much revised version of heatmap() in library(mva)
#       shows a dendrogram, the heatmap, barplots for binary covariates,
#       a dotplot for one intervalscaled varaible, plus color-coding for
#       a chosen cut through the dendrogram
# Auth: ?heatmap for original authors
#       Modifications by Alexander.Ploner@meb.ki.se 211203 onwards
#
# Chng: 221203 AP added more flexibility for layout (do.dendro etc.)
#       260104 AP added return value grp
#       040129 AP added na.rm to quantile, rank
#       060303 AP changed name
#
{
    # initial check up
    scale <- match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    r.cex <- 0.2 + 1/log10(nr)
    c.cex <- 0.2 + 1/log10(nc)
    if (!missing(h) & !missing(clus)) {
        warning("specified both h and clus - ignoring clus")
    }
    # clustering
    if (missing(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (missing(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (!inherits(Rowv, "dendrogram")) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
    }
    else ddr <- Rowv
    if (!inherits(Colv, "dendrogram")) {
        hcc <- hclustfun(distfun(t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
    }
    else ddc <- Colv
    if (reorder[1]) 
        rowInd <- order.dendrogram(ddr)
    else 
        rowInd = 1:nr
    if (reorder[2])         
        colInd <- order.dendrogram(ddc)
    else 
        colInd = 1:nc
    x <- x[rowInd, colInd]
    # Scaling for display
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sd <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sd, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sd <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sd, "/")
    }
    ## Problems downstream with R CMD check re-setting of margins, pin
    ##op <- par()["mar"]
    ##on.exit(par(op))
    # slightly dirty: trim of extremes to make the legend more reasonable
    if (!missing(trim)) {
        trim = min(trim[1], 1-trim[1])
        lo = quantile(x, trim, na.rm=na.rm) 
        hi = quantile(x, 1-trim, na.rm=na.rm)
        x[x<lo] = lo
        x[x>hi] = hi
    }
    # slightly dirty too - should not be used together with trim
    if (equalize) {
        att = attributes(x)
        x = rank(x, na.last=if (na.rm) TRUE else NA)
        attributes(x) = att
    }
    # Set up the plot, what to do etc.
    do.xlabels = !is.null(colnames(x))        
    do.ylabels = !is.null(rownames(x))        
    if (missing(addvar)) {
        extraplot = 0
    } else if (is.null(covariate) | ncol(addvar)<2) {
        extraplot = 1
    } else {
        extraplot = 2
    }
    ll.width = c(1,6)
    ll = ll.height = NULL
    if (do.dendro) {
        ll = rbind(ll, c(0,1), c(0,2))
        ll.height = c(ll.height, 2, 5)
    } else {
        ll = rbind(ll, c(0,1))
        ll.height = c(ll.height, 7)
    }
    cnt = max(ll)
    if (extraplot==1) {
        ll = rbind(ll, c(0, cnt+1))
        ll.height = c(ll.height, 1)
    } else if (extraplot==2) {
        ll = rbind(ll, c(0, cnt+1), c(0, cnt+2))
        ll.height = c(ll.height, 1, 1)
    }
    # now go
    layout(ll, width=ll.width, height=ll.height, respect=TRUE)            
    # fill in from the top, starting with the dendrogram
    mm = c(0,1,3,2)
    par(mar = mm)            
    # our own, unter Schmerzen
    if (do.dendro) {
        grp = oldCutplot.dendrogram(ddc, h, cluscol=cluscol, axes = FALSE, xaxs = "i", leaflab = "none")
    } 
    # the actual plot            
    mm = c(1,1, if (do.xlabels) 3 else 0, 2)
    par(mar=mm)
    image(1:ncol(x), 1:nrow(x), t(x), axes = FALSE, xlim = c(0.5, 
        ncol(x) + 0.5), ylim = c(0.5, nrow(x) + 0.5), xlab = "", 
        ylab = "", col=col, ...)
    box()
    if (do.xlabels) {
        axis(3, 1:ncol(x), las = 2, line = -0.5, tick = 0, labels = colnames(x), cex.axis = c.cex)
    }
    if (do.ylabels) {
        axis(2, 1:nrow(x), las = 2, tick = FALSE, labels = rownames(x), cex.axis = c.cex)
    }
    ## the strip plot
    # if we cluster, we want to sort by cluster number!
    if (!missing(h)) {    
        grp = cutree(hcc, h=h)[colInd]
    } else if (!missing(clus)) {
        grp = clus[colInd]
    } else {
        grp = NULL
    }
    if (!missing(addvar)) {
        par(mar=c(0,1,0,2))        
        oldPicketplot(addvar[colInd, , drop=FALSE], covariate=covariate, add=TRUE, control=picket.control, grp=grp, grpcol=cluscol, grplabel=cluslabel)
    }
    
        
    invisible(list(rowInd = rowInd, colInd = colInd, clus=grp))
}



#' Plot Subtrees of a Dendrogram in Different Colors
#' 
#' Plot a dendrogram, cut the tree at a given height, and draw the resulting
#' subtrees in different colors (OLD version, to be deprecated)
#' 
#' This routine makes use of the functions \code{plot.dendrogram} and
#' \code{plotNode} in package \code{stats}.
#' 
#' @param x a dendrogram.
#' @param h the height at which the dendrogram is cut.
#' @param cluscol the colors used for the subtrees; defaults to \code{rainbow}.
#' @param leaflab indicates how leaf labels are to be drawn< defaults to
#' 'perpendicular'.
#' @param horiz logical indicating whether to plot the dendrogram horizontally
#' or vertically.
#' @param lwd the line width used for the color subtrees.
#' @param \dots arguments to \code{plot.dendrogram}.
#' @author Alexander Ploner <Alexander.Ploner@@ki.se>
#' @seealso \code{\link{as.dendrogram}}
#' @keywords hplot aplot cluster
#' @examples
#' 
#' data(swiss)
#' cc = as.dendrogram(hclust(dist(swiss)))
#' oldCutplot.dendrogram(cc, h=80)
#' 
#' @export oldCutplot.dendrogram
oldCutplot.dendrogram = function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=3, 
                              ...)
#
# Name: OldCutplot.dendrogram
# Desc: takes a dendrogram as described in library(mva), cuts it at level h,
#       and plots the dendrogram with the resulting subtrees in different 
#       colors
# Auth: obviously based on plot.dendrogram in library(mva)
#       modifications by Alexander.Ploner@meb.ki.se  211203
#
# Chng: 050204 AP 
#       changed environment(plot.hclust) to environment(as.dendrogram) to
#       make it work with R 1.8.1
#       250304 AP added RainbowPastel() to make it consistent with oldPicketplot
#       030306 AP slightly more elegant access of plotNode
#       100811 AP name change to avoid collisions
#
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    
    # Not nice, but necessary
    pn  = stats:::plotNode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]], leaflab="none", ...)
    
    x = x[[2]]
    K = length(x)
    if (missing(cluscol)) {
       cluscol = RainbowPastel(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=horiz)
        x1 = x2 + 1
   }
   
    
}



#' Barplots for Several Binary Variables
#' 
#' Display one or more binary variables by using black bars for
#' presence/validity of a condition, empty space for absence/invalidity, and an
#' extra color for missing values. Additionally, an index plot for one interval
#' scaled variable can be added, possibly with a smoothing function (OLD
#' version, to be deprecated).
#' 
#' This routine is primarily intended for augmenting heatmaps. It might be
#' useful in other contexts, but misses most frills for using it comfortably.
#' 
#' The following named list elements can be set to change the appearance of the
#' plot: \describe{ \item{boxw}{the relative width of a marking box.}
#' \item{boxh}{the relative height of a marking box.} \item{hbuff}{the
#' horizontal separation around marking boxes; equals half the horizontal
#' distance between two marking boxes.} \item{vbuff}{ditto for vertical
#' separation.} \item{span}{passed on to \code{loess} used for the smoothing
#' curve.} \item{nacol}{color for missing values of binary variables.}
#' \item{degree}{if 0, no smoothing line is drawn; otherwise passed on to
#' \code{loess} used for the smoothing curve.} \item{cex.label}{the character
#' size for \code{grplabel}.} }
#' 
#' @param x a matrix or data frame containing the data.
#' @param covariate the index of the column in \code{x} that contains the
#' interval scaled variable, if any.
#' @param grp cluster indices for the rows of \code{x}, used for assigning
#' background color.
#' @param grpcol colors corresponding to the clusters.
#' @param grplabel cluster names.
#' @param add logical indicating whether to start a new plot, or whether to add
#' the plot to the existing one.
#' @param control a list of parameters controlling the appearance of the plot,
#' see Details.
#' @note The plot looks like a more or less derelict picket fence, and
#' 'picketplot' sounds somewhat like the 'pocketplot' used in geostatistics.
#' @author Alexander Ploner <Alexander.Ploner@@ki.se>
#' @seealso \code{\link{heatmap_plus}}
#' @keywords hplot aplot
#' @examples
#' 
#' # without covariate
#' mm = cbind(sample(0:1, 42, rep=TRUE), sample(0:1, 42, rep=TRUE))
#' mm[sample(42, 5), 1] = NA
#' oldPicketplot(mm)
#' 
#' # with clustering
#' cl = rep(1:3, c(10,22,10))
#' cn = c("Cluster I","Cluster II","Cluster III")
#' cc = c("lightblue","lightgreen","lightpink") # windows palette
#' oldPicketplot(mm, grp=cl, grplabel=cn, grpcol=cc)
#' 
#' # add a covariate; setting the colnames makes the variable labels
#' mm = cbind(mm, rnorm(42) + cl/2)
#' colnames(mm) = c("State A","State B", "X")
#' oldPicketplot(mm, covariate=3,grp=cl, grplabel=cn, grpcol=cc)
#' 
#' # using extra controls
#' oldPicketplot(mm, covariate=3,grp=cl, grplabel=cn, grpcol=cc, control=list(nacol="white", degree=0))
#' 
#' 
#' @export oldPicketplot
oldPicketplot = function (x, covariate=NULL, grp=NULL, grpcol, grplabel=NULL, 
                       add=FALSE, control=list()) 
#
# Name: oldPicketplot (looks like a picket fence with holes, and sounds like the
#                   pocketplot in geostatistics)
# Desc: visualizes a pattern of 0/1/NAs by using bars, great for annotating a 
#       heatmap
# Auth: Alexander.Ploner@meb.ki.se  181203
#
# Chng: 221203 AP loess() with degree < 2
#       260104 AP 
#       - made loess() optional
#       - better use of space if no covariate
#       030304 AP
#       - added RainbowPastel() as default colors
#       - check grplabel before passing it to axis
#       108011 AP name change to avoid collisions
#
{
    
    # deal with the setup
    cc = list(boxw=1, boxh=4, hbuff=0.1, vbuff=0.1, span=1/3, nacol=gray(0.85), 
              degree=1, cex.label=1.5)
    cc[names(control)] = control
    
    # if we do not add this, we have to set up the plot ourselves
    margin = par()$mar
    if (!add) {
        if (!is.null(covariate)) {
            layout(matrix(c(1,2), nc=1))
            mm = margin
            mm[1] = 0
            par(mar=mm)
        }
    }
    # we have leeway for exactly one covariate 
    x = as.matrix(x)    
    if (!is.null(covariate)) {
        if ((covariate < 1) | (covariate > ncol(x))) {
            stop("wrong index for covariate")
        }
        bindata = x[ , -covariate, drop=FALSE]
        covar   = x[ , covariate]
    } else {
        bindata = x
        covar =NULL
    }
    # here we could check for binary data 
        
    # we assume rows=cases, cols=variables
    k = ncol(bindata)    
    n = nrow(bindata)
    if (k > 0) {
        ww = n*(cc$boxw+2*cc$hbuff)
        hh = k*(cc$boxh+2*cc$vbuff)
        # the data matrix is read columnwise; the rectangles are 
        # plotted row-wise    
        x0 = rep(seq(cc$hbuff, by=cc$boxw+2*cc$hbuff, length=n),k)
        x1 = x0 + cc$boxw
        y0 = rep(seq(cc$vbuff, by=cc$boxh+2*cc$vbuff, length=k), rep(n,k))
        y1 = y0 + cc$boxh
        fill = ifelse(bindata==1, "black", "transparent")
        fill[is.na(fill)] = cc$nacol
        plot(c(0,ww), c(0, hh), type="n", ann=FALSE, xaxs="i", yaxs="i", xaxt="n", yaxt="n")
        # if we get a cluster variable, we have to set differently colored 
        # backgrounds; this assumes that the grp variable is sorted in the 
        # way it appears on the plot
        if (!is.null(grp)) {
            # this should take care of shit
            grp = as.integer(factor(grp, levels=unique(grp)))
            tt = table(grp)
            gg = length(tt)
            xx.grp = c(0,cumsum(tt/sum(tt))*ww)
            if (missing(grpcol)) {
                grpcol=RainbowPastel(gg)
            }
            for (i in 1:gg) {
                rect(xx.grp[i], 0, xx.grp[i+1], hh, col=grpcol[i], border="transparent")
            }
        }
       
        # now the real variables    
        rect(x0,y0,x1,y1, col=fill, border="transparent")
        box()
        # now do the annotations
        label = colnames(bindata)
        if (!is.null(label)) {
            yy = sort(unique((y0+y1)/2))
            axis(2, at=yy, label=label, las=TRUE, font=2, col=par("bg"), col.axis=par("fg"), tick=FALSE)
        }
    }
    # the extra plot
    if (!is.null(covar)) {
        if (!add) {
            mm = margin
            mm[3] = 0
            par(mar=mm)
        }
        plot(1:n, covar, ann=FALSE, xaxs="i", xlim=c(0.5, n+0.5), xaxt="n", yaxt="n", type="n")        
        # the same story for the background
        if (!is.null(grp)) {
            tt = table(grp)
            gg = length(tt)
            xx.grp = c(0.5,0.5+cumsum(tt/sum(tt))*n)
            uu = par("usr")
            for (i in 1:gg) {
                rect(xx.grp[i], uu[3], xx.grp[i+1], uu[4], col=grpcol[i], border="transparent")
            }
        }
        xx = 1:n
        points(xx, covar)
        if ((cc$degree>0) & (cc$span>0)){
            yy = predict(loess(covar~xx, span=cc$span, degree=cc$degree))
            lines(xx, yy)
        }
        axis(1, 1:n, labels=FALSE)
        axis(4)
        box()
        label = colnames(x)[covariate]
        if (!is.null(label)) {
            yy = mean(range(covar))
            axis(2, at=yy, label=label, las=TRUE, tick=FALSE, font=2)
        }
    }
    
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        axis(1, xx.grp, label=FALSE, tcl=-1.5)
        mids = (xx.grp[1:gg] + xx.grp[2:(gg+1)])/2
        # Is the grplabel ok?
        labelnum = length(grplabel)
        if (labelnum < gg) {
            warning("more groups than labels (filling up with blanks)")
            grplabel = c(grplabel, rep(" ", gg-labelnum))
        } else if (gg < labelnum) {
            warning("more labels than groups (ignoring the extras)")
            grplabel = grplabel[1:gg]
        }
        # Go
        axis(1, mids, label=grplabel, font=2, cex.axis=cc$cex.label, tick=FALSE)
    }
            
}



#' Alternative color schemes
#' 
#' \code{RGBColVec} returns a vector of colors that is equally spaced from red
#' through black to green, suitable for heatmaps.
#' 
#' \code{RainbowPastel} returns a vector of colors like \code{rainbow}, but
#' more pastelly.
#' 
#' 
#' @aliases RGBColVec RainbowPastel
#' @param nrgcols,n desired number of colors
#' @param blanche the amount of whiteness added; value between 0 and 255
#' @param \dots extra arguments to \code{rainbow}
#' @return A character vector of length \code{nrgcols} or \code{n} giving the
#' RGB codes for the colors.
#' @author \code{RGBColVec} is based on function \code{rgcolors.func} in
#' package \code{sma} by Sandrine Dudoit and Jane Fridlyand.
#' 
#' \code{RGBColVec} as documented and \code{RainbowPastel} by Alexander Ploner
#' @seealso \code{\link{heat.colors}}
#' @keywords color
#' @examples
#' 
#'  # A Color Wheel
#' pie(rep(1,12), col=RGBColVec(12))
#' 
#' # A color wheel in the original rainbow
#' pie(rep(1,6), col=rainbow(6))
#' 
#' # Pastel
#' pie(rep(1,6), col=RainbowPastel(6))
#' 
#' # Less whiteness
#' pie(rep(1,6), col=RainbowPastel(6, blanche=127))
#' 
#' # More steps require less whiteness
#' pie(rep(1,12), col=RainbowPastel(12, blanche=60))
#' 
#' # Test your screen & eyes: any differences?
#' pie(rep(1,12), col=RainbowPastel(12, blanche=80))
#' 
#' @export RGBColVec
RGBColVec = function (nrgcols=12)
{
# Name: RGBColVec
# Desc: construct a colorvector from red to green
# Auth: Alexander.Ploner@mep.ki.se      110303
# Note: this is an adaption of code found in rgcolors.func() of library(sma)
#
# TODO:
#
# Chng:
#

    k <- trunc(nrgcols/2)
    if (2*k == nrgcols) {   # the even case
        r <- c(rev(seq(0, 1, length = k)),rep(0, k))
        g <- c(rep(0, k),seq(0, 1, length = k))
        colvec <- rgb(r, g, rep(0, 2 * k))
    } else {                # the odd case
        r <- c(rev(seq(1/(2*k-1), 1, length = k)),rep(0, k+1))
        g <- c(rep(0, k+1),seq(1/(2*k-1), 1, length = k))
    colvec <- rgb(r, g, rep(0, 2 * k+1))
    }

    colvec
}

#' @rdname RGBColVec
#' @export RainbowPastel
RainbowPastel =  function (n, blanche=200, ...)
    #
    # Name: RainbowPastel
    # Desc: constructs a rainbow clolr vector, but more pastelly
    # Auth: Alexander.Ploner@mep.ki.se      030304
    #
    # Chng:
    #
    
{
    cv = rainbow(n, ...)
    rgbcv = col2rgb(cv)
    rgbcv = pmin(rgbcv+blanche, 255)
    rgb(rgbcv[1,], rgbcv[2,], rgbcv[3, ], maxColorValue=255)
}
