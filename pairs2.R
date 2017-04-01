
.panelhist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE, breaks=8)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan3", ...)
#    if(length(unique(x)) == 2){
#     text(breaks[c(1,nB-1)]+.05, y[c(1,nB-1)]+.2, c("F","M"))
#    }
}


.panelcor <- function(x, y, digits=2, prefix="", cex.cor)
{
    if(length(unique(y)) == 2) return(NULL)
    else{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    # r <- abs(cor(x, y, use="complete.obs", method='spearman'))
    tmp <- cor.test(x,y, exact=F, method="spearman")
    r <- tmp$estimate; pvalue <- tmp$p.value
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep=""); tcol <- 1
    if(r >= 0.7) tcol <- 'red'  #(tmp$p.value < 0.05)
    else if(r <= -0.7) tcol <- 'blue2'
    text(0.5, 0.5, txt, cex = .8/strwidth(txt), col = tcol, font=(tmp$p.value < 0.05)+1)
    }
}  

.panelmain <- function(x,y,...){
   points(x,y,pch=1,cex=.5, col='#225ea8')
   abline(lm(y~x), col='red')
}


pairs2 <- 
  function (x, labels, panel = points, ..., lower.panel = panel, 
            upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
            label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
            row1attop = TRUE, gap = 1) 
  {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                          oma, ...) {
      if (side%%2 == 1) 
        {par(mgp=c(.5,0,0),tck=-0.1); 
        Axis(x, side = side, xpd = NA, ...)}
#        {vec <- seq(1,length(x),by=2); 
#         print(y[vec])
#         text(x=x[vec],y=rep(min(y),length(vec)),labels=round(x[vec]), pos=side,xpd=TRUE,adj=c(0,0), cex=.4, col='black',srt=30)
#        }
      else {
       par(mgp=c(1.1,.1,0), tck=-0.1)
       Axis(y, side = side, xpd = NA, ...)
      }
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
      x <- as.data.frame(x)
      for (i in seq_along(names(x))) {
        if (is.factor(x[[i]]) || is.logical(x[[i]])) 
          x[[i]] <- as.numeric(x[[i]])
        if (!is.numeric(unclass(x[[i]]))) 
          stop("non-numeric argument to 'pairs'")
      }
    }
    else if (!is.numeric(x)) 
      stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
      lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
      upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
      diag.panel <- match.fun(diag.panel)
    if (row1attop) {
      tmp <- lower.panel
      lower.panel <- upper.panel
      upper.panel <- tmp
      tmp <- has.lower
      has.lower <- has.upper
      has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
      stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
      dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
      dots$main
    else NULL
    if (is.null(oma)) {
      oma <- c(2.3,2.3,0,0)+.1 #c(4, 4, 4, 4)
      if (!is.null(main)) 
        oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), oma = oma, mar=c(0,0,0,0))  #mar = rep.int(gap/2, 4), 
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    for (i in if (row1attop) 
      1L:nc
         else nc:1L) for (j in 1L:nc) {
           localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                     type = "n", ...)
           if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
             box()
             # edited here...
             #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
             #           localAxis(1 + 2 * row1attop, x[, j], x[, i], 
             #                       ...)
             # draw x-axis
             if (i == nc & j != nc & j%%2) 
               localAxis(1, x[, j], x[, i], 
                         ...)
             # draw y-axis
             if (j == 1 & i != 1 & (!i%%2|i==nc)) 
               localAxis(2, x[, j], x[, i], ...)
             #           if (j == nc && (i%%2 || !has.upper || !has.lower)) 
             #             localAxis(4, x[, j], x[, i], ...)
             mfg <- par("mfg")
             if (i == j) {
               if (has.diag) 
                 localDiagPanel(as.vector(x[, i]), ...)
               if (has.labs) {
                 par(usr = c(0, 1, 0, 1))
                 if (is.null(cex.labels)) {
                   l.wid <- strwidth(labels, "user")
                   cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                 }
                 text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                            font = font.labels)
               }
             }
             else if (i < j) 
               localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                              i]), ...)
             else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                                 i]), ...)
             if (any(par("mfg") != mfg)) 
               stop("the 'panel' function made a new plot")
           }
           else par(new = FALSE)
         }
    if (!is.null(main)) {
      font.main <- if ("font.main" %in% nmdots) 
        dots$font.main
      else par("font.main")
      cex.main <- if ("cex.main" %in% nmdots) 
        dots$cex.main
      else par("cex.main")
      mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
  }
  
  
  
addAxisLabels <- function(){
x.coords = par('usr')[1:2]
y.coords = par('usr')[3:4]

# Offset is estimated distance between edge of plot area and beginning of actual plot
x.offset = 0.0 * (x.coords[2] - x.coords[1])  
xrng =  (x.coords[2] - x.coords[1]) - 2*x.offset
x.width = xrng/p*1.9

y.offset = 0.0 * (y.coords[2] - y.coords[1])
yrng =  (y.coords[2] - y.coords[1]) - 2*y.offset
y.width = yrng/p*2.5

# seq function below calculates the location of the midpoint of each panel

# x-axis labels

# # seting for PDF
#text(seq(x.coords[1] + x.offset + 0.5*x.width, x.coords[2] - x.offset - 0.5*x.width,length.out=p)+.015, .01, Xnams,xpd=TRUE,adj=c(0,0), cex=.5, col='black',srt=30)

text(seq(x.coords[1] + x.offset + 0.5*x.width, x.coords[2] - x.offset - 0.5*x.width,length.out=p)+.015, .01, Xnams,xpd=TRUE,adj=c(0,0), cex=.5, col='black',srt=30)

# y-axis labels
# # seting for PDF
#text(.017, seq(y.coords[1] + y.offset + 0.5*y.width, y.coords[2] - 3*y.offset - 0.5*y.width, length.out=p)+.032,rev(Xnams), xpd=TRUE, adj=c(0.5, 0.5), srt=45, cex=.5,col='black')  

text(.017, seq(y.coords[1] + y.offset + 0.5*y.width, y.coords[2] - 3*y.offset - 0.5*y.width, length.out=p)+.032,rev(Xnams), xpd=TRUE, adj=c(0.5, 0.5), srt=45, cex=.5,col='black')  

}  