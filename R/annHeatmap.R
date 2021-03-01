## Imports (very fine-grained, could prbably be simplified)
#' @importFrom graphics axis
#' @importFrom graphics box
#' @importFrom graphics image
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics rect
#' @importFrom graphics layout.show
#' @importFrom graphics segments
#' 
#' @importFrom grDevices col2rgb
#' @importFrom grDevices gray
#' @importFrom grDevices heat.colors
#' @importFrom grDevices rainbow
#' @importFrom grDevices rgb
#' 
#' @importFrom stats as.dendrogram
#' @importFrom stats cutree
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats loess
#' @importFrom stats order.dendrogram
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats reorder
#' @importFrom stats dendrapply
#' @importFrom stats model.matrix
#' @importFrom stats na.exclude
#' @importFrom stats naresid
#' 
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
NULL


#' Generate a layout for an (annotated) heatmap
#' 
#' Generate a layout for an (annotated) heatmap. This function will generally
#' not be called directly, but only via \code{annHeatmap2}.
#' 
#' Space for plots is reserved via the \code{layout} mechanism. The function
#' starts with an empty maximum layout, fills in the plot, dendrograms,
#' annotation plots and legend as required, and compresses the resulting layout
#' by removing empty slots.
#' 
#' @param dendrogram A list with named entries \code{Row} and \code{Col}. Each
#' of these is a list with a named entry \code{status}. If the value of
#' \code{status} is the string \code{"yes"}, space will be set aside for
#' drawing a row- and/or column dendrogram.
#' @param annotation A list with named entries \code{Row} and \code{Col}. Each
#' of these is a list with a named entry \code{data}. If the value of
#' \code{data} is not \code{NULL}, space will be set aside for a picket plot
#' showing the row- and/or column annotation.
#' @param leg.side An integer indicating on where to reserve space for the
#' legend: values 1-4 correspond to below, to the left, above and to the right,
#' as in e.g. \code{axis}. For a value of \code{NULL}, the function provides a
#' reasonable default where there is space left in the layout. For any other
#' value, no space for a legend is put aside.
#' @param show A logical value; if \code{TRUE}, the layout defined by the
#' arguments is displayed graphically.
#' @return A list with the following entries: \item{plot }{A matrix describing
#' the plot layout; see \code{layout}} \item{width }{relative widths of plots
#' (i.e. columns)} \item{height }{relative heights of plots (i.e. rows)}
#' \item{legend.side }{side where to draw the legend}
#' @seealso \code{\link{annHeatmap2}}, \code{\link{picketPlot}},
#' \code{\link{layout}}
#' @keywords utilities
#' @examples
#'     def.par = par(no.readonly = TRUE) # save default, for resetting
#' 
#'     ## Heatmap with column dendrogram, column annotation, default legend
#'     dnd = list(Row=list(status="no"), Col=list(status="yes"))
#'     ann = list(Row=list(data=NULL), Col=list(data=1))
#'     ## 1 = heatmap, 2=dendrogram, 3=annotation, 4=legend
#'     ll = heatmapLayout(dendrogram=dnd, annotation=ann, leg.side=NULL, show=TRUE)
#'     ll
#'     
#'     par(def.par)  #- reset to default
#' 
#' @export heatmapLayout
heatmapLayout = function(dendrogram, annotation, leg.side=NULL, show=FALSE)
{
    ## Start: maximum matrix, 5 x 5, all zero
    ## Names for nice display post ante
    ll = matrix(0, nrow=5, ncol=5)
    ll.width = ll.height = rep(0, 5)
    cnt = 1
    rownames(ll) = c("leg3", "colDendro","image", "colAnn", "leg1")
    colnames(ll) = c("leg2", "rowDendro","image", "rowAnn", "leg4")
 
    ## The main plot
    ll[3,3] = cnt
    ll.width[3] = ll.height[3] = 5
    cnt = cnt+1
    ## The column dendrogram
    if (dendrogram$Col$status=="yes") {
        ll[2, 3] = 2
        ll.width[3]  = 5
        ll.height[2] = 2
        cnt = cnt+1
    }
    ## The row dendrogram
    if (dendrogram$Row$status=="yes") {
        ll[3, 2] = cnt
        ll.width[2]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    # Column annotation
    if (!is.null(annotation$Col$data)) {
        ll[4, 3] = cnt
        ll.width[3] = 5
        ll.height[4] = 2
        cnt = cnt+1
    }
    ## Row annotation
    if (!is.null(annotation$Row$data)) {
        ll[3, 4] = cnt
        ll.width[4]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    ## Legend: if no pref specified, go for empty, if possible
    if (is.null(leg.side)) {
        if (dendrogram$Row$status != "yes") {
            leg.side = 2
        } else if (is.null(annotation$Row$data)) {
            leg.side = 4
        } else if (is.null(annotation$Col$data)) {
            leg.side = 1
        } else if (dendrogram$Col$status != "yes") {
            leg.side = 3
        } else {
            leg.side = 4
        }
    }
    ## Add the legend space
    if (leg.side==1) {
        ll[5,3] = cnt
        ll.width[3] = 5
        ll.height[5] = 1
    } else if (leg.side==2) {
        ll[3,1] = cnt
        ll.width[1] = 1
        ll.height[3] = 5
    } else if (leg.side==3) {
        ll[1,3] = cnt
        ll.width[3]  = 5
        ll.height[1] = 1
    } else if (leg.side==4) {
        ll[3,5] = cnt
        ll.width[5]  = 1
        ll.height[3] = 5
    }
    
    ## Compress
    ndx = rowSums(ll)!=0
    ll  = ll[ndx, , drop=FALSE]
    ll.height = ll.height[ndx]
    ndx = colSums(ll)!=0
    ll  = ll[, ndx, drop=FALSE]
    ll.width = ll.width[ndx]
    ## Do it - show it
    if (show) {
        layout(ll, widths=ll.width, heights=ll.height, respect=TRUE)            
        layout.show(max(ll))
    }
    return(list(plot=ll, width=ll.width, height=ll.height, legend.side=leg.side))
}

#' Override existing list entries
#' 
#' Override existing list entries and extract arguments that are specified as
#' named lists
#' 
#' \code{modifyExistingList} is a general function that recursively overwrites
#' named items in \code{x} with the value of items of \code{val} with the same
#' name. Items in \code{val} that have no name, or do not correspond to an item
#' in \code{x} with the same name, are ignored.
#' 
#' \code{extractArg} is a specific helper function for setting default values
#' for the \code{annHeatmap2}-family of functions, where arguments are given as
#' a list with two named items, \code{Row} and \code{Col}. Each of these items
#' is again a named list of actual parameters. At the same time, all items with
#' other names than \code{Row} and \code{Col} at the top level are assumed to
#' be shared items with the same value for both sub-lists. \code{extractArg}
#' uses \code{modifyExistingList} to overwrite the default values specified in
#' \code{deflist} with the actual values specified in \code{arglist}, see
#' Examples.
#' 
#' @aliases modifyExistingList extractArg
#' @param x a named list, the target for replacing with entries with the same
#' name from \code{val}
#' @param val a named list that serves as template for filling in values in
#' \code{x}
#' @param arglist a named list; these are the specified arguments that override
#' the defaults.
#' @param deflist a named list whose entries are all possible slots (with
#' default values) that can be filled.
#' @return \code{modifyExistingList} returns \code{x}, with values replaced
#' from \code{val} where names match. \code{extractArg} returns a list with
#' items \code{Row} and \code{Col} fully specified according to both
#' \code{deflist} and \code{arglist}.
#' @seealso \code{\link{annHeatmap2}}
#' @keywords utilities
#' @examples
#' 
#'     ## Replace items with matching names recursively
#'     x   = list(a=1, b=2, c=list(a=31, b=32), 135)
#'     val = list(a=2, c=list(b=1114), d=92)
#'     modifyExistingList(x, val)
#'     
#'     ## Same defaults for rows/columns, no arguments specified
#'     defs = list(a="A", b="B", c="C")
#'     extractArg(NULL, defs)
#' 
#'     ## Shared and non-shared defaults
#'     defs = list(common.1=134, common.2=72, Row=list(row.only=14), Col=list(col.only=134))
#'     args = list(common.1 = -1, Row=list(row.only=94, common.2=-15))
#'     extractArg(args, defs)
#' @export modifyExistingList
modifyExistingList = function(x, val)
{
    if (is.null(x)) x = list()
    if (is.null(val)) val = list()
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    for (v in intersect(xnames, vnames)) {
        x[[v]] <- if (is.list(x[[v]]) && is.list(val[[v]])) 
            modifyExistingList(x[[v]], val[[v]])
        else val[[v]]
    }
    x
}


#' @rdname modifyExistingList
#' @export extractArg
extractArg = function(arglist, deflist)
{
    if (missing(arglist)) arglist = NULL
    al2 = modifyExistingList(deflist, arglist)
    row = col = al2
    row = modifyExistingList(row, arglist[["Row"]])    
    col = modifyExistingList(col, arglist[["Col"]])    
    list(Row=row, Col=col)  
}


#' Display a data frame of annotation information
#' 
#' Displays a data frame of both factor and numerical variables in parallel
#' panels. Factors levels are indicated by black rectangles, using dummy
#' variables for more than two levels. Numerical variables are shown as simple
#' index plots with an optional loess smoother. Panels can be arranged
#' horizontally or vertically, and different groups of subjects can be
#' indicated through different background colors.
#' 
#' Missing values are indicated by a box marking in \code{nacol} for factor
#' values.
#' 
#' @param x usually a data frame, which is passed to \code{convAnnData} to be
#' converted to a numerical matrix with dummy coding for the factor levels.
#' Alternatively, such a numerical matrix can be constructed manually and
#' passed in as \code{x}, see argument \code{asIs}
#' @param grp an optional vector of cluster memberships, in the same order as
#' the rows of \code{x}
#' @param grpcol an optional vector of background colors for the clusters
#' specified in \code{grp}
#' @param grplabel an optional vector of names for the clusters specified in
#' \code{grp}
#' @param horizontal logical value whether to plot variables horizontally
#' (default) or vertically
#' @param asIs a logical value indicating whether \code{x} should be passed to
#' \code{convAnnData} for pre-processing or not. Defaults to \code{FALSE}.
#' @param control a named list of control parameters that determines the visual 
#' appearance of the plot; see \code{picketPlotControl} for details.
#' 
#' @return Invisibly, a list containing the data and parameters used for plotting
#' each binary indicator and numerical variable, respectively. This is an internal
#' data structure, mostly useful for debugging. Irrelevant, as the main desired
#' effetc is a plot to the current graphical device. 
#' 
#' @seealso \code{\link{annHeatmap2}}, \code{\link{convAnnData}},
#' \code{\link{par}}, \code{\link{picketPlotControl}}
#' @keywords hplot
#' @examples
#' 
#'     ## Standard call
#'     data(mtcars)
#'     picketPlot(mtcars)
#'     
#'     ## Pre-process the data for display
#'     mm = convAnnData(mtcars, inclRef=FALSE)
#'     picketPlot(mm, asIs=TRUE)
#'     
#'     ## Higher panels for continous traits
#'     picketPlot(mm, asIs=TRUE, control=list(numfac=3))
#'  
#'     ## With clusters
#'     picketPlot(mtcars, grp = rep(1:2, c(16, 16)), grpcol = c("pink","lightblue"), grplabel=c("Cluster 1", "Cluster 2"))
#' 
#' @export picketPlot
picketPlot = function (x, grp=NULL, grpcol, grplabel=NULL, horizontal=TRUE, asIs=FALSE, control=list()) 
#
# Name: picketPlot (looks like a picket fence with holes, and sounds like the
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
#       2010-07-08 AP
#       - complete re-write
#       2010-08-28
#       - re-arranged code for vertical/horizontal drawing
{    
    # deal with the setup
    cc = picketPlotControl()
    cc[names(control)] = control
    
    ## Convert/check the data
    x = convAnnData(x, asIs=asIs)
    
    # Count variables, panels, types
    nsamp  = nrow(x)
    npanel = ncol(x)
    bpanel = apply(x, 2, function(y) all(y[is.finite(y)] %in% c(0,1)) )
    
    # Compute panel heights, widths
    panelw = nsamp*(cc$boxw+2*cc$hbuff)
    panelh = cc$boxh+2*cc$vbuff
    totalh = sum(panelh * ifelse(bpanel, 1, cc$numfac))
    LL = cbind(0, 0)
    UR = cbind(panelw, totalh)
     
    # Set up the x-values for a single panel
    xbase = seq(cc$hbuff, by=cc$boxw+2*cc$hbuff, length=nsamp)
    xcent = xbase + cc$boxw/2
    
    # if we get a cluster variable, we have to set differently colored 
    # backgrounds; this assumes that the  grp variable is sorted in the 
    # way it appears on the plot
    if (!is.null(grp)) {
        grp = as.integer(factor(grp, levels=unique(grp)))
        tt = table(grp)
        gg = length(tt)
        grpcoord = c(0,cumsum(tt/sum(tt))*panelw)
        grp0 = cbind(grpcoord[1:gg], rep(0, gg))
        grp1 = cbind(grpcoord[2:(gg+1)], rep(totalh, gg))
        if (missing(grpcol)) {
            grpcol=BrewerClusterCol
        }
        if (is.function(grpcol)) grpcol = grpcol(gg)
        ## In case of manually specified group colors, we only check/use the
        ## relevant colors and ignore the rest
        grpcol = grpcol[1:gg]
        if (gg > 1) {
            if ( any(grpcol[-1] == grpcol[-gg]) ) warning("neighboring clusters with same color, potentially misleading")
        }
    }

    # Loop over vars and fill in the panels
    panels = list()
    voff = 0
    for (i in 1:npanel) {
        if (bpanel[i]) {
            ## Coordinates
            x0 = xbase
            x1 = x0+cc$boxw
            y0 = voff + cc$vbuff
            y1 = y0   + cc$boxh
            ## Set fill
            fill = ifelse(x[, i, drop=FALSE]==1, "black", "transparent")
            fill[is.na(fill)] = cc$nacol
            label = colnames(x)[i]
            labcc = if (!is.null(label)) (y0+y1)/2 else NULL 
            panels[[i]] = list(ll=cbind(x0, y0), ur=cbind(x1, y1), fill=fill, label=label, labcc=labcc)
            voff = voff + panelh
        } else {
            xv = x[,i]
            rr = range(xv, na.rm=TRUE)
            yval = voff + cc$vbuff*cc$numfac + ((xv - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            if ((cc$degree>0) & (cc$span>0)){
                yy = predict(loess(yval~xcent, span=cc$span, degree=cc$degree))
            } else {
                yy = rep(NA, length(xcent))
            }
            label = colnames(x)[i]
            labcc = if (!is.null(label)) mean(range(yval, na.rm=TRUE)) else NULL
            ## Shrink the axis range for the labels, if specified
            if ( (cc$label_axis_shrink >= 0) & (cc$label_axis_shrink <= 1) ) {
                axis_offset <- cc$label_axis_shrink * (rr[2] - rr[1])/2
            } else {
                axis_offset <- 0
            }
            axis_range <- c(rr[1] + axis_offset, rr[2] - axis_offset)
            axlab = pretty(axis_range)
            axcc  = voff + cc$vbuff*cc$numfac + ((axlab - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            panels[[i]] = list(raw = cbind(xcent, yval), smo = cbind(xcent, yy), 
                               label = label, labcc = labcc, 
                               axlab = axlab, axcc = axcc, 
                               low_hi = c(voff, voff + panelh*cc$numfac))
            voff = voff + panelh*cc$numfac
        }
    }
    
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        mids = (grpcoord[1:gg] + grpcoord[2:(gg+1)])/2
        # Is the grplabel ok?
        labelnum = length(grplabel)
        if (labelnum < gg) {
            warning("more groups than labels (filling up with blanks)")
            grplabel = c(grplabel, rep(" ", gg-labelnum))
        } else if (gg < labelnum) {
            warning("more labels than groups (ignoring the extras)")
            grplabel = grplabel[1:gg]
        }
    }
            
    ## Switch coordinates, if you have to
    h2v = function(cc) cbind(cc[,2]-totalh, cc[,1])
    if (horizontal) {
        grpaxis = 1
        labaxis = 2
        covaxis = 4
        las = 1        
        draw_box = function(x) abline(h = x)
    } else {
        grpaxis = 4
        labaxis = 3
        covaxis = 1
        las = 3
        draw_box = function(x) abline(v = x - totalh)
        ## Rotate
        LL = h2v(LL)
        UR = h2v(UR)
        if (!is.null(grp)) {
            grp0 = h2v(grp0)
            grp1 = h2v(grp1)
        }
        for (i in 1:npanel) {
            panels[[i]][[1]] = h2v(panels[[i]][[1]])
            panels[[i]][[2]] = h2v(panels[[i]][[2]])
            panels[[i]]$labcc = panels[[i]]$labcc - totalh 
            panels[[i]]$axcc  = panels[[i]]$axcc - totalh
        }
    }
    
    # Set up the plot
    plot(rbind(LL, UR), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
    # Add the colored rectangles, if required
    if (!is.null(grp)) {
        rect(grp0[,1], grp0[,2], grp1[,1], grp1[,2], col=grpcol, border="transparent")
    }
    # Loop over vars and fill in the panels
    for (i in 1:npanel) {
        if (bpanel[i]) {
            ## Do the rectangles
            with(panels[[i]], rect(ll[,1], ll[,2], ur[,1], ur[,2], col=fill, border="transparent") )
        } else {
            with(panels[[i]], points(raw[,1], raw[,2], pch=cc$pch, cex=cc$cex.pch, col=cc$col.pch))
            if ((cc$degree>0) & (cc$span>0)){
                with(panels[[i]], lines(smo[,1], smo[,2]))
            }
            with(panels[[i]], axis(covaxis, at=axcc, labels=axlab, las = las))
            ## Draw the base reference line
            if (cc$plot_baseline) draw_box(panels[[i]]$low_hi)
        }
        ## Name panel (regardless of type)
        if (!is.null(panels[[i]]$label)) {
            axis(labaxis, at=panels[[i]]$labcc, labels=panels[[i]]$label, las=las, tick=FALSE, font=2, col=par("bg"), col.axis=par("fg"))
        }
    }
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        axis(grpaxis, grpcoord, labels=FALSE, tcl=-1.5)
        axis(grpaxis, mids, labels=grplabel, font=2, cex.axis=cc$cex.label, tick=FALSE)
    }                                
    invisible(panels)
}


#' Default parameter settings for picketPlot
#' 
#' This function returns a named list of parameters that affect how a
#' picketPlot is generated. This list can be used as a template for 
#' overriding the defaults partially or completely.
#' 
#' @details The following parameter affects the overall appearance of the plot:
#' \itemize{ 
#'   \item\code{cex.label} is the expansion factor for the size of the cluster labels at
#'     the bottom of the plot; default is 1.5.
#' }
#' 
#' The following parameters directly affect how binary indicator variables are 
#' displayed:
#' \itemize{ 
#'   \item\code{boxw} is the relative length of the short side of a box marking 
#'      (width for a horizontal plot); default is 1.
#'   \item\code{boxh} is the relative length of the long side of a box marking
#'     (default: 4) 
#'   \item\code{hbuff} is the relative distance between two box markings for the 
#'     same variable (horizontal buffer for a horizontal plot); default is 0.1 
#'   \item\code{vbuff} is the relative distance between two box
#'      markings for the same subject, but different variables (default: 0.1)
#'   \item\code{nacol} is the color for box markings indicating missing values 
#'     (default: \code{gray(0.85)}) 
#'  }
#'  Note that \code{boxh} and \code{vbuff} also affect the display of numerical
#'  variables as a scatter plot: as the amount of vertical space allowed for a 
#'  single numerical variable (see also \code{numfac} below) and the vertical 
#'  space between two neighboring variable panels (binary or ornumerical), 
#'  respectively.
#'
#' The following parameters only affect the display of a numerical variable:
#' \itemize{
#'   \item\code{numfac} is the expansion factor indicating how much higher (for a
#'     horizontal plot) or wider (for a vertical plot) panels with numerical
#'     variables are than a panels for a single binary indicator
#'   \item\code{span} is the span argument for the loess smoother. Default is 1/3; 
#'     setting this to zero switches off smoothing.
#'   \item\code{degree} is the degree of loess smoothing. Default is 1; setting
#'     this to zero switches off smoothing 
#'   \item\code{pch} is the plotting character for numerical variables; uses the 
#'     device default. 
#'   \item\code{cex.pch} is the size of the plotting character for numerical 
#'     variables; uses the device default. 
#'   \item\code{col.pch} is the color of the plotting character for numerical 
#'     variables; uses the device default.
#'   \item\code{label_axis_shrink} controls the range of the axis for which axis 
#'     ticks are labeled: by default, labels covering the whole observed range 
#'     are defined (via a call to \code{pretty}); if set to a number between zero
#'     and one, the range covered by labels is shortened by that fraction and 
#'     centered within the observed range (and the fed to \code{pretty}); this can 
#'     be used to avoid overlapping labels for multiple adjacent panes with 
#'     numerical variables. 
#'   \item\code{plot_baseline} is a logical value indicating whether to draw 
#'     a baseline for panes showing numerical variables: FALSE by default, this
#'     can be useful to visually separate multiple adjacent panes with numerical
#'     variables. 
#' }
#'  
#' @return A named list
#' @seealso \code{\link{picketPlot}}, \code{\link{par}}, \code{\link{pretty}}
#' @export picketPlotControl
picketPlotControl = function()
{
    list(cex.label = 1.5,
         boxw = 1, boxh = 4, hbuff = 0.1, vbuff = 0.1, nacol = gray(0.85), 
         span = 1/3, degree = 1, numfac = 2, 
         pch = par("pch"), cex.pch = par("cex"), col.pch = par("col"),
         label_axis_shrink = 0, plot_baseline = FALSE 
    )
}

#' Get nice (symmetric) breaks for an interval
#' 
#' Given a minimum and a maximum, this function returns a vector of equidistant
#' breaks that covers this interval, and has a pretty interval length (1, 2, or
#' 5 times a power of 10). If the interval contains zero, it will be one of the
#' breaks, so that the intervals are arranged somewhat symmetrically around it.
#' 
#' The number of desired breaks is honored as far as possible, which is not
#' actually that often in practice. However, major deviations of three or more
#' are reasonably rare.
#' 
#' The functiona allows the specification of a set of breaks instead of the
#' desired number of breaks, somewhat like in \code{cut}. However, if the
#' length of \code{breaks} is greater than one, the function just sorts the
#' values and returns them otherwise unchanged.
#' 
#' @param xr the range to be covered, as \code{c(min, max)}
#' @param breaks either the desired number of breaks, or a pre-specified vector
#' of breaks
#' @return A vector of pretty breaks covering the specified interval, more or
#' less of the desired length.
#' @seealso \code{\link{pretty}}
#' @keywords utilities
#' @examples
#' 
#'     
#'     ## Niceness overrules specified number
#'     niceBreaks(c(-1,1), 5)
#'     niceBreaks(c(-1,1), 6)
#' 
#'     ## Zero appears always as break
#'     niceBreaks(c(-2.75, 1.12), 8)
#' 
#'     ## Not invariant to translation (of course)
#'     niceBreaks(3.27 + c(-2.75, 1.12), 8)
#' 
#'     
#' 
#' @export niceBreaks
niceBreaks = function(xr, breaks)
{
    ## If you want it, you get it
    if (length(breaks) > 1) {
        return(sort(breaks))
    }
    ## Ok, so you proposed a number
    ## Neg and pos?
    if ( (xr[1] < 0) & (xr[2] > 0) ) {
        xminAbs = abs(xr[1])
        xmax    = xr[2]
        nneg = max(round(breaks * xminAbs/(xmax+xminAbs)), 1)
        npos = max(round(breaks * xmax/(xmax+xminAbs)), 1)
        nbr  = pretty(c(xr[1], 0), nneg)
        pbr  = pretty(c(0, xr[2]), npos)
        ## Average of the proposed interval lengths,
        ##  nice enough for us
        diff = ( (nbr[2]-nbr[1]) + (pbr[2] - pbr[1]) ) / 2
        nbr  = diff * ( (xr[1] %/% diff)  : 0 ) 
        pbr  = diff * ( 1 : (xr[2] %/% diff + 1) )
        breaks = c(nbr, pbr)
    } else { ## only pos or negs
        breaks = pretty(xr, breaks)
    }
    breaks
}


#' Color palette for (symmetric) breaks
#' 
#' Given a vector of breaks specifying a set of intervals, this function
#' provides a vector of colors for the indicating the intervals graphically. If
#' the intervals are arranged symmetrically around a specified value, the
#' colors try to reflect this.
#' 
#' The meaning of symmetrical is rather generous here: it is enough that the
#' intervals specified by \code{breaks} are of equal length and that
#' \code{center} is one of the breaks. This means we allow for more or less
#' intervals on one side of \code{center}.
#' 
#' This really only works well if \code{colors} is specified as
#' \code{g2r.colors}, which returns a symmetrical color vector (from green to
#' red) if an even number of colors is requested. The whole point is then that
#' if there are more classes to one side of \code{center} than to the other,
#' this will be reflected by deeper shades of red or green on the appropriate
#' side.
#' 
#' @param breaks a vector of breaks
#' @param colors either an explicit vector of colors, or a palette function
#' that takes a number and returns a vector of colors
#' @param center optional center around which to check for symmetry
#' @param tol tolerance (as relative error) for deviation from mathematically
#' exact symmetry
#' @return A vector of colors, of length one less than the number of breaks.
#' @seealso \code{\link{g2r.colors}}
#' @keywords utilities
#' @examples
#' 
#'     ## Fully symmetrical breaks
#'     br1 = (-3) : 3
#'     co1 = breakColors(br1, g2r.colors)
#'     co1
#'     doLegend(br1, co1, 1)
#'     
#'     ## Truncated on one side
#'     br2 = (-2) : 4
#'     co2 = breakColors(br2, g2r.colors)
#'     co2
#'     doLegend(br2, co2, 1)
#'     
#'     ## Does not work with other color schemes
#'     co3 = breakColors(br2, heat.colors)
#'     co3
#'     doLegend(br2, co3, 1)    
#' 
#' @export breakColors
breakColors = function(breaks, colors, center=0, tol=0.001)
{
    ## In case of explicit color definitions
    nbreaks = length(breaks)
    nclass  = nbreaks - 1
    if (!is.function(colors)) {
        ncolors = length(colors)
        if (ncolors > nclass) {
            warning("more colors than classes: ignoring ", ncolors-nclass, " last colors")
            colors = colors[1:nclass]
        } else if (nclass > ncolors) {
            stop(nclass-ncolors, " more classes than colors defined")
        }
    } else {
        ## Are the classes symmetric and of same lengths?
        clens = diff(breaks)
        aclen = mean(clens)
        if (aclen==0) stop("Dude, your breaks are seriously fucked up!")
        relerr = max((clens-aclen)/aclen)
        if ( (center %in% breaks) & (relerr < tol) ) { ## yes, symmetric
            ndxcen = which(breaks==center)
            kneg = ndxcen -1 
            kpos = nbreaks - ndxcen
            kmax = max(kneg, kpos)
            colors = colors(2*kmax)
            if (kneg < kpos) {
                colors = colors[ (kpos-kneg+1) : (2*kmax) ]
            } else if (kneg > kpos) {
                colors = colors[ 1 : (2*kmax - (kneg-kpos)) ]
            }
        } else {                                      ## no, not symmetric
            colors = colors(nclass)
        }
    }
    colors
}




#' Palette from green to red via black
#' 
#' Returns a color vector of the requested length, ranging from pure red to
#' pure green via slighlty tinted black.
#' 
#' If \code{n} is even, the colors range from pure green to green-tinted black
#' to red-tinted black to pure red. If \code{n} is odd, the colors range from
#' pure red to pure green, with full black for the median class.
#' 
#' @param n the number of requested colors
#' @param min.tinge the proportion of red/green added to black to make it
#' recognizably green or red
#' @return A vector of (RGB-) colors of the specified length
#' @seealso \code{\link{breakColors}}
#' @keywords utilities
#' @examples
#' 
#'     ## Even number: residual tint shows left/right of center
#'     co_even = g2r.colors(10)
#'     co_even
#'     doLegend(1:11, co_even, 1)
#'     
#'     ## Odd number: central class all black
#'     co_odd = g2r.colors(9)
#'     co_odd
#'     doLegend(1:10, co_odd, 1)
#' 
#'     ## Lighter tint in the middle
#'     co_light = g2r.colors(10, min.tinge=0.50)
#'     co_light
#'     doLegend(1:11, co_light, 1)
#' 
#' @export g2r.colors
g2r.colors = function(n=12, min.tinge = 0.33)
{
    k <- trunc(n/2)
    if (2 * k == n) {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, k))
        r <- c(rep(0, k), seq(min.tinge, 1, length = k))
        colvec <- rgb(
        r, g, rep(0, 2 * k))
    }
    else {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, 
            k + 1))
        r <- c(rep(0, k + 1), seq(min.tinge, 1, length = k))
        colvec <- rgb(r, g, rep(0, 2 * k + 1))
    }
    colvec
}



#' A simple legend
    #' 
#' Add a simple legend in form of a color bar to a plot.
#' 
#' This is an extremely simple way of giving a visual impression of what
#' numerical values correspond to a given color. The actual plot is done via a
#' call to \code{\link{image}} and \code{\link{axis}}.
#' 
#' @param breaks a vector of breaks defining a set of intervals for the data
#' @param col a vector of colors corresponding to the intervals.
#' @param side integer between 1 and 4, indicating on which side of the main
#' plot the legend is supposed to be drawn. Standard interpretation: 1 = below,
#' continuing clock-wise.
#' @return The locations of the ticks returned by the call to
#' \code{\link{axis}}
#' @seealso \code{\link{plot.annHeatmap}}, \code{\link{niceBreaks}},
#' \code{\link{g2r.colors}}
#' @keywords utilities
#' @examples
#' 
#'     ## Set up data
#'     doLegend(1:9, g2r.colors(8), 2)
#' 
#' @export doLegend
doLegend = function(breaks, col, side)
{
    zval = ( breaks[-1] + breaks[-length(breaks)] ) / 2
    z  = matrix(zval, ncol=1)
    if (side %in% c(1,3)) {
        image(x=zval, y=1, z=z, xaxt="n", yaxt="n", col=col, breaks=breaks , xaxs="i", xlab="", ylab="")
    } else {
        image(x=1, y=zval, z=t(z), xaxt="n", yaxt="n", col=col, breaks=breaks, yaxs="i", xlab="", ylab="")
    }        
    axis(side, las=1)
}




#' Converting data frames for display as annotation
#' 
#' Converts a data frames for display as annotation in a heatmap. This is
#' mostly intended as an internal function, but might be useful for finetuning
#' an annotation data frame manually.
#' 
#' Logical variables are converted to factors. So are numerical variables with
#' less than \code{nval.fac} unique values.
#' 
#' @param x the data frame to be converted
#' @param nval.fac lower limit for unique values in numerical variables
#' @param inclRef logical value indicating whether to include the reference
#' level among the dummy variables for factors
#' @param asIs logical value indicating whether to perform a conversion; if
#' \code{TRUE}, the input \code{x} is simply returned, provided it is a
#' numerical matrix (otherwise, the function stops with an error message)
#' @return \code{convAnnData} returns the converted data frame, which is a
#' numerical matrix
#' @seealso \code{\link{annHeatmap2}}
#' @keywords utilities
#' @examples
#' 
#' data(mtcars)
#' summary(mtcars)
#' summary(convAnnData(mtcars))
#' summary(convAnnData(mtcars, nval.fac=2))    
#' summary(convAnnData(mtcars, nval.fac=2, inclRef=FALSE))    
#' 
#' 
#' @export convAnnData
convAnnData = function(x, nval.fac=3, inclRef=TRUE, asIs=FALSE)
{
    if (is.null(x)) return(NULL)
    if (asIs) {
        if (is.matrix(x) & is.numeric(x)) return(x)
        else stop("argument x not a numerical matrix, asIs=TRUE does not work")
    }
    
    x = as.data.frame(x)
    if (!is.null(nval.fac) & nval.fac>0) doConv = TRUE
    vv = colnames(x)
    for (v in vv) {
        if (is.logical(x[,v])) {
            x[,v] = factor(as.numeric(x[,v]))
        }
        if (doConv & length(unique(x[is.finite(x[,v]),v])) <= nval.fac) {
            x[,v] = factor(x[,v])    
        }
    }
    ret  = NULL
    ivar = 0
    for (v in vv) {
         xx = x[, v]
         if (is.factor(xx)) {
            nandx = is.na(xx)
            if (length(unique(xx[!nandx])) > 1) {
                naAction = attr(na.exclude(x[, v, drop=FALSE]), "na.action")
                modMat   = model.matrix(~xx-1)
                if (!inclRef) modMat = modMat[ , -1, drop=FALSE]
                binvar   = naresid(naAction, modMat)
                colnames(binvar) = paste(v, "=", levels(xx)[if (!inclRef) -1 else TRUE], sep="")
            } else {
                nlev = length(levels(xx))
                ilev = unique(as.numeric(xx[!nandx]))
                if (length(ilev)==0) {
                    binvar = matrix(NA, nrow=length(xx), ncol=nlev)
                } else {
                    binvar = matrix(0, nrow=length(xx), ncol=nlev)
                    binvar[, ilev] = 1
                    binvar[nandx, ] = NA
                }
                colnames(binvar) = paste(v, "=", levels(xx), sep="")
            }
            ret = cbind(ret, binvar)
            ivar = ivar + ncol(binvar)
         } else {
            ret = cbind(ret, x[,v])
            ivar = ivar + 1
            colnames(ret)[ivar] = v
         }
    } 
    ret
}


###################################################
### code chunk number 12: cut.dendrogram_Def
###################################################
cutree.dendrogram = function(x, h)
{
    # Cut the tree, get the labels
    cutx = cut(x, h)
    cutl = lapply(cutx$lower, getLeaves)
    # Set up the cluster vector as seen in the plot
    nclus = sapply(cutl, length)
    ret   = rep(1:length(nclus), nclus)
    # Return cluster membership in the order of the original data, if possible
    ord = order.dendrogram(x) 
    # Is the order a valid permutation of the data?
    if (!all(sort(ord)==(1:length(ret)))) {
        stop("dendrogram order does not match number of leaves - is this a subtree?")
    }
    # Ok, proceed
    ret[ord] = ret
    ret = as.integer(factor(ret, levels=unique(ret))) # recode for order of clus
    names(ret)[ord] = unlist(cutl)
    ret
}



#' Undocumented functions
#' 
#' These functions are currently undocumented. Please refer to the source and
#' source comments if you feel you need to use them.
#' 
#' @aliases getLeaves
#' @rdname Undocumented
#' @keywords internal
getLeaves = function(x)
{
    unlist(dendrapply(x, function(x) attr(x, "label")))
}


#' Printing information about annotated heatmaps
#' 
#' Printing method for annotated heatmaps
#' 
#' A very simple printing method, displaying a minimum of information about
#' dendrograms and annotation
#' 
#' @param x an object of class \code{annHeatmap}
#' @param \dots extra arguments, currently ignored
#' @return \code{x} is returned invisibly
#' @seealso \code{\link{annHeatmap}}, \code{\link{annHeatmap2}},
#' \code{\link{plot.annHeatmap}}
#' @keywords hplot
#' @examples
#' 
#'     set.seed(219)
#'     mat = matrix(rnorm(100), ncol=5)
#'     ann = data.frame(Class=c("A","A","B","A","B"))
#'     map1 = annHeatmap(mat, ann)
#'     map1
#' 
#' @export print.annHeatmap
#' @rawNamespace S3method(print, annHeatmap) 
print.annHeatmap = function(x, ...)
{
    cat("annotated Heatmap\n\n")
    cat("Rows: "); print(x$dendrogram$Row$dendro)
    cat("\t", if (is.null(x$annotation$Row$data)) 0 else ncol(x$annotation$Row$data), " annotation variable(s)\n")
    cat("Cols: "); print(x$dendrogram$Col$dendro)
    cat("\t", if (is.null(x$annotation$Col$data)) 0 else ncol(x$annotation$Col$data), " annotation variable(s)\n")
    invisible(x)
}


#' Color scheme for clusters
#' 
#' This function returns a color vector based on one of the qualitative
#' paletters supported by \code{RColorBrewer}. This allows visually distinct
#' coloring of clusters and ensures sure that adjacent clusters have different
#' colors.
#' 
#' This is just a wrapper for \code{\link{brewer.pal}} that checks that the
#' specified palette is qualitative, and allows for an arbitrary number of
#' colors: for less than three colors, it just returns the first and second
#' colors of the palette; for more than \code{maxcolors} colors, it recycles
#' the basic palette as often as required. This is ok, because the main point
#' is to have different colors for neighboring clusters.
#' 
#' @param n desired number of colors
#' @param name name of the qualitative palette from which colors are taken, see
#' \code{\link{brewer.pal.info}}
#' @return A character vector of length \code{n} of hexadecimal color codes.
#' @seealso \code{\link{brewer.pal}}
#' @keywords color
#' @examples
#' 
#' ## A Color Wheel: default palette with maximum number of colors
#' pie(rep(1,9), col=BrewerClusterCol(9))
#' 
#' ## Double the number of colors 
#' pie(rep(1,18), col=BrewerClusterCol(18))
#' 
#' ## Only two clusters/colors
#' pie(rep(1,2), col=BrewerClusterCol(2))
#' 
#' ## Different qualitative palette: stronger colors
#' pie(rep(1,12), col=BrewerClusterCol(12, "Paired"))
#' 
#' 
#' @export BrewerClusterCol
BrewerClusterCol = function(n, name="Pastel1")
{
    ## Check the name of the palette
    qualpal = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category=="qual", ]
    name = match.arg(name, rownames(qualpal))
    nmax = qualpal[name, "maxcolors"]

    ## Get the full color vector of the palette
    cols = RColorBrewer::brewer.pal(nmax, name)

    ## Build the (shortened or recycled) index vector
    ndx  = rep(1:nmax, length=n)

    cols[ndx]
}


cutplot.dendrogram = function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=3, ...)
#
# Name: cutplot.dendrogram
# Desc: takes a dendrogram as described in library(mva), cuts it at level h,
#       and plots the dendrogram with the resulting subtrees in different 
#       colors
# Auth: obviously based on plot.dendrogram in library(mva)
#       modifications by Alexander.Ploner@meb.ki.se  211203
#
# Chng: 050204 AP 
#       changed environment(plot.hclust) to environment(as.dendrogram) to
#       make it work with R 1.8.1
#       250304 AP added RainbowPastel() to make it consistent with picketplot
#       030306 AP slightly more elegant access of plotNode
#       220710 AP also for horizontal plots
#       120811 AP use edgePar instead of par() for col and lwd
#
{
    ## If there is no cutting, we plot and leave
    if (missing(h) | is.null(h)) {
        return(plot(x, leaflab=leaflab, horiz=horiz, edgePar=list(lwd=lwd), ...))
    }
    ## If cut height greater than tree, don't cut, complain and leave
    treeheight = attr(x, "height")
    if (h >= treeheight) {
        warning("cutting height greater than tree height ", treeheight, ": tree uncut")
        return(plot(x, leaflab=leaflab, horiz=horiz, edgePar=list(lwd=lwd), ...))
    }
            
    ## Some param processing
    if (missing(cluscol) | is.null(cluscol)) cluscol = BrewerClusterCol
    
    # Not nice, but necessary
    pn  = stats:::plotNode
       
    x = cut(x, h)
    plot(x[[1]], leaflab="none", horiz=horiz, edgePar=list(lwd=lwd), ...)
    
    x = x[[2]]
    K = length(x)
    if (is.function(cluscol)) {
       cluscol = cluscol(K)
    }
    left = 1
    for (k in 1:K) {
        right = left + attr(x[[k]],"members")-1
        if (left < right) {         ## not a singleton cluster 
            pn(left, right, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(lwd=lwd, col=cluscol[k]), horiz=horiz)
        } else if (left == right) { ## singleton cluster
            if (!horiz) {
                segments(left, 0, left, h, lwd=lwd, col=cluscol[k])
            } else {
                segments(0, left, h, left, lwd=lwd, col=cluscol[k])
            }
        } else stop("this totally should not have happened")
        left = right + 1
   }

}



#' Annotated heatmaps
#' 
#' This function plots a data matrix as intensity heatmap, with optional
#' dendrograms, annotation panels and clustering for both rows and columns.
#' This is the actual working function called by numerous wrappers.
#' 
#' Arguments \code{scale}, \code{breaks}, \code{col} and \code{legend} control
#' different aspects of the whole plot directly as described. Arguments
#' \code{dendrogram}, \code{annotation}, \code{cluster} and \code{labels}
#' control aspects that may differ for the rows and columns of the central
#' heatmap and have a special structure: each is a named list with different
#' entries controling e.g. the look of a dendrogram, the data for annotation
#' etc. Additionally, they can contain two extra entries called simply
#' \code{Row} and \code{Col}; these are again named lists that can contain all
#' the same entries as the parent list. Entries specified directly in the list
#' apply to both rows and columns; entries specified as part of \code{Row} or
#' \code{Col} override these defaults for the rows or columns only.
#' 
#' Recognized parameters for argument \code{dendrogram}: \describe{
#' \item{clustfun}{the clustering function for generating the dendrogram;
#' defaults to \code{hclust} for rows and columns} \item{distfun}{a function
#' that returns the pairwise distances between samples/features as an object of
#' class \code{dist}; defaults to \code{dist} for rows and columns}
#' \item{status}{a string that controls the display of the dendrogram:
#' \code{yes} means use the dendrogram to re-order the rows/columns and display
#' the dendrogram; \code{hidden} means re-rorder, but do not display; \code{no}
#' means do not use the dendrogram at all.} \item{lwd}{the line width of the
#' branches of the dendrogram; defaults to 3.} \item{dendro}{an override
#' argument that allows to pass in a dendrogram directly, bypassing the
#' \code{clustfun} and \code{distfun} mechanism; defaults to \code{NULL} (i.e.
#' is not used)} }
#' 
#' Recognized entries for argument \code{annotation}: \describe{ \item{data}{a
#' data frame containing the annotation data; defaults to \code{NULL}, i.e. no
#' annotation is displayed} \item{control}{a list of fine-tuning parameters
#' that is passed directly to \code{picketPlot}; defaults to an empty list,
#' i.e. the default settings in \code{picketPlot}} \item{asIs}{logical value
#' indicating whether the annotation \code{data} needs to be pre-processed via
#' \code{convAnnData} or not; defaults to \code{TRUE}} \item{inclRef}{logical
#' value indicating whether to include all levels of factor variables in
#' \code{data}, or whether to drop the reference level (i.e. the first level).
#' Defaults to \code{TRUE}} }
#' 
#' Recognized entries for argument \code{cluster}: \describe{ \item{cuth}{the
#' height at which to cut through the dendrogram to define groups of similar
#' features/samples; defaults to \code{NULL}, i.e. no cutting}
#' \item{label}{labels for the clusters; defaults to \code{NULL}, i.e. no
#' labels} \item{col}{colors for the different clusters; the colors are used
#' for coloring both the sub-trees of the dendrogram and the corresponding area
#' in the annotation plot (if there is one). This is either a vector of colors,
#' or a palette function that takes a number and returns a vector of colors of
#' the specified length; defaults to \code{BrewerClusterCol}} \item{grp}{an
#' override argument that directly specifies group memberships for the
#' features/samples, completely bypassing the whole \code{dendrogram} and
#' \code{cuth} mechanism. This probably only works for
#' \code{dendrogram$status="no"}.} }
#' 
#' Recognized entries for argument \code{labels}: \describe{ \item{cex}{size of
#' the text for the labels; defaults to \code{NULL}, i.e. use a hard-coded
#' default guess} \item{nrow}{amount of space available for the labels between
#' the central heatmap and the dendrogram, expressed as lines of text; defaults
#' to 3.} \item{side}{side at which to draw the labels, coded as integer
#' between 1 and 4 in the usual way (1 = below the plot, continuing clockwise).
#' A common default for rows and columns does not make sense: rows only work
#' with 2 and 4, columns only with 1 and 3. Defaults try to make use of empty
#' space, depending on the presence of a dendrogram.} \item{labels}{labels for
#' rows and columns; defaults to \code{NULL}, i.e. using the row- and column
#' names of \code{x}. Note that these labels are applied \emph{after}
#' re-sorting rows and columns as per dendrogram, so these have to be already
#' sorted accordingly. If you want to change the labels \emph{before}
#' re-sorting, it is is easier to re-set the row- and/or column names of
#' \code{x}.} }
#' 
#' @param x the numerical matrix to be shown as heatmap
#' @param dendrogram a list that controls how row- and column diagrams are
#' determined and displayed
#' @param annotation a list that controls the data and the way it is shown in
#' row- and column annotation panels
#' @param cluster a list that controls how many clusters are chosen, and how
#' these clusters are labelled and colored
#' @param labels a list that controls the row- and column labels, as well as
#' their size and placement
#' @param scale a character string indicating how the matrix \code{x} is
#' standardized (by row, by column or not at all). This affects only display,
#' not dendrograms or clustering
#' @param breaks specifies the interval breaks for displaying the data in
#' \code{x}; either a vector of explicit interval breaks, or just the desired
#' number of intervals. See \code{niceBreaks} for details.
#' @param col specifies a palette of colors for the heatmap intensities; either
#' a vector of explicit color definitions (one less than breaks) or a palette
#' function. See \code{breakColors}.
#' @param legend whether and where to draw a legend for the colors/intervals in
#' the heatmap. If \code{TRUE}, a legend is placed in a position determined by
#' the function to be suitable; alternatively, integer values 1-4 indicate the
#' side where the legend is to be drawn; and \code{FALSE} indicates that no
#' legend should be drawn.
#' @return An object of class \code{annHeatmap}. Use \code{plot} to display it
#' graphically. 
#' @seealso \code{\link{heatmapLayout}}, \code{\link{niceBreaks}},
#' \code{\link{breakColors}}, \code{\link{g2r.colors}}, \code{\link{BrewerClusterCol}}
#' @keywords hplot
#' @examples
#' 
#' require(Biobase)
#' data(sample.ExpressionSet)
#' ex1 = sample.ExpressionSet[51:85,]
#' map1 = annHeatmap2(exprs(ex1), ann=list(Col=list(data=pData(ex1))),
#'                    cluster=list(Col=list(cuth=3000)))
#' plot(map1)
#' 
#' @export annHeatmap2
annHeatmap2 = function(x, dendrogram, annotation, cluster, labels, scale=c("row", "col", "none"), breaks=256, col=g2r.colors, legend=FALSE)
#
# Name: annHeatmap2
# Desc: a (possibly doubly) annotated heatmap
# Auth: Alexander.Ploner@ki.se 2010-07-12
#
# Chng: 
#
{
    ## Process arguments
    if (!is.matrix(x) | !is.numeric(x)) stop("x must be a numeric matrix")
    nc = ncol(x); nr = nrow(x)
    if (nc < 2 | nr < 2) stop("x must have at least two rows/columns")
    
    ## Process the different lists: dendrogram, cluster, annotation
    ## See lattice:::xyplot.formula, modifyLists, lattice:::construct.scales
    def = list(clustfun=hclust, distfun=dist, status="yes", lwd=3, dendro=NULL)
    dendrogram = extractArg(dendrogram, def)
    def = list(data=NULL, control=picketPlotControl(), asIs=FALSE, inclRef=TRUE)
    annotation = extractArg(annotation, def)
    def = list(cuth=NULL, grp=NULL, label=NULL, col=BrewerClusterCol)
    cluster = extractArg(cluster, def)
    def = list(cex=NULL, nrow=3, side=NULL, labels=NULL)
    labels = extractArg(labels, def)
    
    ## Check values for the different lists

    
    ## Generate the layout: TRUE means default, FALSE means none
    ## Otherwise, integer 1-4 indicates side
    if (is.logical(legend)) {
        if (legend) leg = NULL else leg = 0
    } else {
        if (!(legend %in% 1:4)) stop("invalid value for legend: ", legend)
        else leg=legend
    }
    layout = heatmapLayout(dendrogram, annotation, leg.side=leg)
    
    ## Copy the data for display, scale as required
    x2 = x
    scale = match.arg(scale)
    if (scale == "row") {
        x2 = sweep(x2, 1, rowMeans(x, na.rm = TRUE))
        sd = apply(x2, 1, sd, na.rm = TRUE)
        x2 = sweep(x2, 1, sd, "/")
    }
    else if (scale == "col") {
        x2 = sweep(x2, 2, colMeans(x, na.rm = TRUE))
        sd = apply(x2, 2, sd, na.rm = TRUE)
        x2 = sweep(x2, 2, sd, "/")
    }
    
    ## Construct the breaks and colors for display
    breaks = niceBreaks(range(x2, na.rm=TRUE), breaks)
    col    = breakColors(breaks, col)
    
    ## Generate the dendrograms, if required; re-indexes in any cases
    ## We could put some sanity checks on the dendrograms in the else-branches
    ## FIXME: store the names of the functions, not the functions in the object
    dendrogram$Row = within(dendrogram$Row, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(x))
            dendro = reorder(as.dendrogram(dendro), rowMeans(x, na.rm=TRUE))
        }
    )
    dendrogram$Col = within(dendrogram$Col, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(t(x)))
            dendro = reorder(as.dendrogram(dendro), colMeans(x, na.rm=TRUE))
        }
    )
    ## Reorder the display data to agree with the dendrograms, if required
    rowInd = with(dendrogram$Row, if (status!="no") order.dendrogram(dendro) else 1:nr)
    colInd = with(dendrogram$Col, if (status!="no") order.dendrogram(dendro) else 1:nc)
    x2 = x2[rowInd, colInd]
    
    ## Set the defaults for the sample/variable labels
    labels$Row = within(labels$Row, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nr)
        if (is.null(side)) side = if (is.null(annotation$Row$data)) 4 else 2
        if (is.null(labels)) labels = rownames(x2)
    })
    labels$Col = within(labels$Col, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nc)
        if (is.null(side)) side = if (is.null(annotation$Col$data)) 1 else 3
        if (is.null(labels)) labels = colnames(x2)
    })
    
    ## Generate the clustering, if required (cut, or resort the cluster var)
    ## FIXME: does not deal with pre-defined grp form outside
    cluster$Row = within(cluster$Row, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Row$dendro, cuth)[rowInd]
        })
    cluster$Col = within(cluster$Col, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Col$dendro, cuth)[colInd]
        })

    ## Process the annotation data frames (factor/numeric, re-sort?)
    annotation$Row = within(annotation$Row, {
        data = convAnnData(data, asIs=asIs, inclRef=inclRef)
    })
    annotation$Col = within(annotation$Col, {
        data = convAnnData(data, asIs=asIs, inclRef=inclRef)
    })

        
    ## Generate the new object
    
    ## print, return invisibly
    ret = list(data=list(x=x, x2=x2, rowInd=rowInd, colInd=colInd, breaks=breaks, col=col), dendrogram=dendrogram, cluster=cluster, annotation=annotation, labels=labels, layout=layout, legend=legend)
    class(ret) = "annHeatmap"
    ret

}


#' Plotting method for annotated heatmaps
#' 
#' Plotting method for annotated heatmaps
#' 
#' This function displays an annotated heatmap object that has been previously
#' generated by \code{annHeatmap2} or on of its wrappers. The arguments
#' \code{widths} and \code{heights} work as in \code{layout}.
#' 
#' @param x an object of class \code{annHeatmap}
#' @param widths a numerical vector giving the widths of the sub-plots
#' currently defined
#' @param heights a numerical vector giving the heights of the sub-plots
#' currently defined
#' @param \dots extra graphical parameters, currently ignored
#' @return \code{x}, invisibly returned. If \code{widths} or \code{heights}
#' have been specified, they overwrite the corresponding items
#' \code{x$layout$width} and \code{x$layout$height} in \code{x}.
#' @seealso \code{\link{annHeatmap2}}, \code{\link{heatmapLayout}},
#' \code{\link{layout}}
#' @keywords hplot
#' @examples
#' 
#'     ## Define the map
#'     require(Biobase)
#'     data(sample.ExpressionSet)
#'     ex1 = sample.ExpressionSet[51:85,]
#'     map1 = annHeatmap2(exprs(ex1), ann=list(Col=list(data=pData(ex1))),
#'                    cluster=list(Col=list(cuth=3000)))
#'     
#'     ## Plot it               
#'     plot(map1)
#'     
#'     ## More heatmap, smaller dendrogram/annotation
#'     map2 = plot(map1, heights = c(1,6,1))
#'     
#'     ## Compare layout before/after
#'     with(map1$layout, layout(plot, width, height))
#'     layout.show(4)
#'     with(map2$layout, layout(plot, width, height))
#'     layout.show(4)
#' 
#' @export plot.annHeatmap
#' @rawNamespace S3method(plot, annHeatmap)
plot.annHeatmap = function(x, widths, heights, ...)
{
	## Preserve parameters that are set explicitly below
	## Not doing this has lead to Issue 8: inconsistent distance
	## between dendrogram and heatmap after repeated calls on same device
	opar = par("oma", "mar", "xaxs", "yaxs")
	on.exit(par(opar))
    ## If there are cluster labels on either axis, we reserve space for them
    doRClusLab = !is.null(x$cluster$Row$label) 
    doCClusLab = !is.null(x$cluster$Col$label)
    omar = rep(0, 4)
    if (doRClusLab) omar[4] = 2
    if (doCClusLab) omar[1] = 2
    par(oma=omar)
    ## Set up the layout    
    if (!missing(widths)) x$layout$width = widths
    if (!missing(heights)) x$layout$height = heights    
    with(x$layout, layout(plot, width, height, respect=TRUE))
    
    ## Plot the central image, making space for labels, if required
    nc = ncol(x$data$x2); nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels) 
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
        par(mar=mmar)
        image(1:nc, 1:nr, t(x2), axes = FALSE, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), xlab = "", ylab = "", col=col, breaks=breaks, ...)    
    })
    with (x$labels, {
        if (doRlab) axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, labels = Row$labels, cex.axis = Row$cex)
        if (doClab) axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, labels = Col$labels, cex.axis = Col$cex)
    })

    ## Plot the column/row dendrograms, as required
    with(x$dendrogram$Col,
        if (status=="yes") {
            par(mar=c(0, mmar[2], 3, mmar[4]))
            cutplot.dendrogram(dendro, h=x$cluster$Col$cuth, cluscol=x$cluster$Col$col, horiz=FALSE, axes = FALSE, xaxs = "i", leaflab = "none", lwd=x$dendrogram$Col$lwd)
        })
    with(x$dendrogram$Row,
        if (status=="yes") {
            par(mar=c(mmar[1], 3, mmar[3], 0))
            cutplot.dendrogram(dendro, h=x$cluster$Row$cuth, cluscol=x$cluster$Row$col, horiz=TRUE, axes = FALSE, yaxs = "i", leaflab = "none", lwd=x$dendrogram$Row$lwd)
        })

    ## Plot the column/row annotation data, as required
    if (!is.null(x$annotation$Col$data)) {
        par(mar=c(1, mmar[2], 0, mmar[4]), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Col$data[x$data$colInd, ,drop=FALSE],
          grp=x$cluster$Col$grp, grplabel=x$cluster$Col$label, grpcol=x$cluster$Col$col,
          control=x$annotation$Col$control, asIs=TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
        par(mar=c(mmar[1], 0, mmar[3], 1), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Row$data[x$data$rowInd, ,drop=FALSE],
          grp=x$cluster$Row$grp, grplabel=x$cluster$Row$label, grpcol=x$cluster$Row$col,
          control=x$annotation$Row$control, asIs=TRUE, horizontal=FALSE)
    }

    ## Plot a legend, as required
    if (x$legend) {
        if (x$layout$legend.side %in% c(1,3)) {
            par(mar=c(2, mmar[2]+2, 2, mmar[4]+2))
        } else {
            par(mar=c(mmar[1]+2, 2, mmar[3]+2, 2))            
        }       
        doLegend(x$data$breaks, col=x$data$col, x$layout$legend.side)
    }    
    
    invisible(x)    
    
}


#' Regular heatmaps with a legend
#' 
#' Creating regular heatmaps, without annotation, but allowing for a legend
#' 
#' A gelded wrapper for \code{annHeatmap2} that allows for heatmaps without
#' annotation or clustering on the dendrograms, but still offer some control
#' over dendrograms, labels and legend.
#' 
#' These functions generate an object representing the heatmap; in order to
#' produce graphical output, you have to invoke the \code{plot} method, see
#' Examples.
#' 
#' @aliases regHeatmap regHeatmap.default
#' @param x a numerical matrix
#' @param dendrogram a list controlling the options for row- and column
#' dendrogram, see \code{annHeatmap2}
#' @param labels a list controlling the row- and column labels as well as their
#' location and size, see \code{annHeatmap2}
#' @param legend either a logical value, indicating whether to draw a legend at
#' the default location determined by the function, or one of the sides of the
#' plot (1-4), see \code{annHeatmap2}
#' @param \dots extra options passed to \code{annHeatmap2}
#' @return An object of class \code{annHeatmap}
#' @seealso \code{\link{annHeatmap}}, \code{\link{annHeatmap2}},
#' \code{\link{plot.annHeatmap}}
#' @keywords hplot
#' @examples
#' 
#'     
#'     ## Default
#'     set.seed(219)
#'     mat = matrix(rnorm(100), ncol=5)
#'     map1 = regHeatmap(mat)
#'     plot(map1)
#'   
#' 
#' @export regHeatmap
regHeatmap = function(x, ...) UseMethod("regHeatmap")


#' Annotated heatmaps
#' 
#' Creating heatmaps with annotated columns
#' 
#' These functions generate an object representing the heatmap; in order to
#' produce graphical output, you have to invoke the \code{plot} method, see
#' Examples.
#' 
#' @aliases annHeatmap annHeatmap.default annHeatmap.ExpressionSet
#' @param x either a numerical matrix with the data for the central heatmap
#' (for the default method) or an object of class \code{ExpressionSet}
#' @param annotation a data frame containing the annotation for the columns of
#' \code{x}
#' @param dendrogram a list controlling the options for row- and column
#' dendrogram, see \code{annHeatmap2}
#' @param cluster a list controlling the options for clustering rows and
#' columns of \code{x}, see \code{annHeatmap2}
#' @param labels a list controlling the row- and column labels as well as their
#' location and size, see \code{annHeatmap2}
#' @param legend either a logical value, indicating whether to draw a legend at
#' the default location determined by the function, or one of the sides of the
#' plot (1-4), see \code{annHeatmap2}
#' @param \dots extra options passed to \code{annHeatmap2}
#' @return An object of class \code{annHeatmap}
#' @section Warning: These are currently simple convenience functions that
#' allow quick plotting, but little control over the finer details. This may
#' change in the future, but for now, if you want to do anything fancy, you
#' should invoke \code{annHeatmap2} directly.
#' @seealso \code{\link{annHeatmap2}}, \code{\link{plot.annHeatmap}}
#' @keywords hplot
#' @examples
#' 
#' 
#'     ## Default method
#'     set.seed(219)
#'     mat = matrix(rnorm(100), ncol=5)
#'     ann = data.frame(Class=c("A","A","B","A","B"))
#'     map1 = annHeatmap(mat, ann)
#'     plot(map1)
#'     
#'     ## Expression set
#'     require(Biobase)
#'     data(sample.ExpressionSet)
#'     map2 = annHeatmap(sample.ExpressionSet)
#'     plot(map2)
#' 
#' @export annHeatmap
annHeatmap = function(x, ...) UseMethod("annHeatmap")


#' @rdname regHeatmap
#' @export regHeatmap.default
#' @rawNamespace S3method(regHeatmap, default)
regHeatmap.default = function(x, dendrogram=list(clustfun=hclust, distfun=dist, status="yes"), labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=NULL, cluster=NULL,  labels=labels, legend=legend, ...)
    ret
}


#' @rdname annHeatmap
#' @export annHeatmap.default
#' @rawNamespace S3method(annHeatmap, default)
annHeatmap.default = function(x, annotation, dendrogram=list(clustfun=hclust, distfun=dist, Col=list(status="yes"), Row=list(status="hidden")), cluster=NULL, labels=NULL, legend=TRUE, ...)
{
	if (!is.data.frame(annotation)) stop("Argument 'annoation' needs to be data frame")

    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=list(Col=list(data=annotation, fun=picketPlot)), cluster=cluster,  labels=labels, legend=TRUE, ...)
    ret
}


#' @rdname annHeatmap
#' @export annHeatmap.ExpressionSet
#' @rawNamespace S3method(annHeatmap, ExpressionSet) 
annHeatmap.ExpressionSet = function(x, ...)
{
    expmat = Biobase::exprs(x)
    anndat = Biobase::pData(x)
    annHeatmap(expmat, anndat, ...)
}


