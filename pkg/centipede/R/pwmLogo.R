## Adapted from seqLogo for using regular polygons instead of the grid thing.
letterA <- function (x.pos, y.pos, ht, wt, id = NULL) 
{
    x <- c(0, 4,  6, 10, 8, 6.8, 3.2, 2, 0,NA, 3.6,   5, 6.4, 3.6)
    y <- c(0,10, 10,  0, 0,   3,   3, 0, 0,NA,   4, 7.5,   4, 4)
    x <- 0.1 * x
    y <- 0.1 * y
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    if (is.null(id)) {
        id <- c(rep(1, 9),NA, rep(2, 4))
    }
    else {
        id <- c(rep(id, 9),NA, rep(id + 1, 4))
    }
    fill <- c("green", "white")
    list(x = x, y = y, id = id, fill = fill)
}
letterC <- function (x.pos, y.pos, ht, wt, id = NULL) 
{
    angle1 <- seq(0.3 + pi/2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    if (is.null(id)) {
        id <- rep(1, length(x))
    }
    else {
        id <- rep(id, length(x))
    }
    fill <- "blue"
    list(x = x, y = y, id = id, fill = fill)
}
letterG <- function (x.pos, y.pos, ht, wt, id = NULL) 
{
    angle1 <- seq(0.3 + pi/2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    h1 <- max(y.l1)
    r1 <- max(x.l1)
    h1 <- 0.4
    x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
    y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
    if (is.null(id)) {
        id <- c(rep(1, length(x)),NA, rep(2, length(x.add)))
    }
    else {
        id <- c(rep(id, length(x)),NA, rep(id + 1, length(x.add)))
    }
    x <- c(rev(x),NA, x.add)
    y <- c(rev(y),NA, y.add)
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    fill <- c("orange", "orange")
    list(x = x, y = y, id = id, fill = fill)
}
letterT <- function (x.pos, y.pos, ht, wt, id = NULL) 
{
    x <- c(0, 10, 10, 6, 6, 4, 4, 0)
    y <- c(10, 10, 9, 9, 0, 0, 9, 9)
    x <- 0.1 * x
    y <- 0.1 * y
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    if (is.null(id)) {
        id <- rep(1, 8)
    }
    else {
        id <- rep(id, 8)
    }
    fill <- "red"
    list(x = x, y = y, id = id, fill = fill)
}
addLetter <- function (letters, which, x.pos, y.pos, ht, wt) 
{
    if (which == "A") {
        letter <- letterA(x.pos, y.pos, ht, wt)
    }
    else if (which == "C") {
        letter <- letterC(x.pos, y.pos, ht, wt)
    }
    else if (which == "G") {
        letter <- letterG(x.pos, y.pos, ht, wt)
    }
    else if (which == "T") {
        letter <- letterT(x.pos, y.pos, ht, wt)
    }
    else {
        stop("which must be one of A,C,G,T")
    }
    letters$x <- c(letters$x,NA,letter$x)
    letters$y <- c(letters$y,NA,letter$y)
    lastID <- ifelse(is.null(letters$id), 0, max(letters$id))
    letters$id <- c(letters$id,NA, lastID + letter$id)
    letters$fill <- c(letters$fill, letter$fill)
    letters
}
pwm2ic <- function (pwm) 
{
    npos <- ncol(pwm)
    ic <- numeric(length = npos)
    for (i in 1:npos) {
        ic[i] <- 2 + sum(sapply(pwm[, i], function(x) {
            if (x > 0) {
                x * log2(x)
            } else {
                0
            }
        }))
    }
    ic
}

pwmLogo <- function (pwm, x.off = 0.5 , y.off = 0.0 , x.sc=1,y.sc=1,ic.scale = TRUE, add=FALSE, axis = TRUE,...)
{
    if (class(pwm) == "pwm") {
        pwm <- pwm@pwm
    }
    else if (class(pwm) == "data.frame") {
        pwm <- as.matrix(pwm)
    }
    else if (class(pwm) != "matrix") {
        stop("pwm must be of class matrix or data.frame")
    }
    if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) 
        stop("Columns of PWM must add up to 1.0")
    chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(pwm)
    if (ic.scale) {
        ylim <- 2
        ylab <- "Information content"
        facs <- pwm2ic(pwm)
    }
    else {
        ylim <- 1
        ylab <- "Probability"
        facs <- rep(1, npos)
    }
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        column <- pwm[, j]
        hts <- 0.95 * column * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0){
                ##letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
                letters <- addLetter(letters, letter, x.pos, y.pos, ht, wt)
              }
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }
##    browser()
    ##iid <- which(diff(letters$id)==1)
    ##letters$x[iid] <- NA
    ##letters$y[iid] <- NA
    letters$x <- (letters$x)*x.sc+x.off;
    letters$y <- (letters$y)*y.sc+y.off;
    if(add==FALSE){
      plot.new()
      plot.window(xlim=c(1,npos),ylim=c(0,2))
    }
    polygon(letters,col=letters$fill,...)
    str(letters)
    if(axis==TRUE){
      axis(1,at=1:npos)
      axis(2)
    }
}

## Sequence LOGO
## IdxPos <- which((R1$PostPr > 0.99))
## IdxNeg <- which((R1$LogRatios < 0.01))
## PwmBound <- ComputePWM(Sequence[IdxPos,11:51])
## pwmLogo(PwmBound,border=NA)

