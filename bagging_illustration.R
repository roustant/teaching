## Author: O. Roustant

# generate the data

N <- 30
set.seed(0)
x <- runif(N, min=-10, max=40)
fun <- function(x){
  0.5*pmax(x-20, 0) + pmax(16-x, 0)
}

y <- fun(x) + 2*rnorm(N)
par(mfrow = c(1,1))
plot(x, y)
xnew <- seq(-10, 40, length.out = 200)
lines(xnew, fun(xnew), lty="dotted")

cex.points <- 1
pdfPlot <- FALSE

## estimation d'un arbre de regression
library(tree)
mtree <- tree(y~x, data = data.frame(x=x, y=y))
p <- predict(mtree, data.frame(x = xnew))
if (pdfPlot) pdf("Tree1D.pdf", width = 8, height = 8)
par(mfrow = c(2,1))
plot(mtree); text(mtree)
plot(x, y)
lines(xnew, p, col = "blue")
lines(xnew, fun(xnew), lty="dotted")
if (pdfPlot) dev.off()

## effect of boostrapping on the tree
treeBoot <- function(x, y, control = tree.control(length(x), ...), 
                     plot = TRUE, plot.OOB = FALSE, size.adapt = TRUE, ...){
  N <- length(x)
  indBoot <- sample(1:N, N, replace = TRUE)
  # print(indBoot)
  xBoot <- x[indBoot]
  yBoot <- y[indBoot]
  mtreeBoot <- tree(y~x, data = data.frame(x=xBoot, y=yBoot),
                    control = control)
  pBoot <- predict(mtreeBoot, data.frame(x = xnew))
  if (plot) {
    if (size.adapt){
      OOB <- setdiff(1:N, indBoot)
      indBootSort <- sort(indBoot, index.return = TRUE)
      indBoot <- indBootSort$x
      xBoot <- xBoot[indBootSort$ix]
      yBoot <- yBoot[indBootSort$ix]
      repBoot <- rle(indBoot)$lengths
    } else {
      repBoot <- 1
    }
    points(x, y, cex = cex.points)
    points(unique(xBoot), unique(yBoot), pch = 19, cex = cex.points * repBoot, ...)
    repText <- as.character(repBoot)
    repText[repBoot == 1] <- ""
    text(unique(xBoot), unique(yBoot), repText, col = "white")
    # verification : rep(unique(xBoot), times = repBoot) - xBoot
    lines(xnew, pBoot, ...)
    if (plot.OOB) arrows(x0 = x[OOB], x1 = x[OOB], 
                         y0 = y[OOB], y1 = predict(mtreeBoot, data.frame(x = x[OOB])),
                         code = 0, col = "red", lwd = 3)
  }
  invisible(pBoot)
}

par(mfrow = c(1,1))

if (pdfPlot) pdf("Bagging0.pdf", width = 10, height = 6)
plot(x, y, cex = cex.points)
lines(xnew, fun(xnew), lty="dotted")
if (pdfPlot) dev.off()
#treeBoot(x, y)

set.seed(0)
if (pdfPlot) pdf("Bagging1.pdf", width = 10, height = 6)
plot(x, y, cex = cex.points)
lines(xnew, fun(xnew), lty="dotted")
treeBoot(x, y, control = tree.control(nobs = length(x), mincut = 1, minsize = 2), 
         col = "blue")
if (pdfPlot) dev.off()
set.seed(0)
if (pdfPlot) pdf("Bagging1OOB.pdf", width = 10, height = 6)
plot(x, y, cex = cex.points)
lines(xnew, fun(xnew), lty="dotted")
treeBoot(x, y, control = tree.control(nobs = length(x), mincut = 1, minsize = 2), 
         col = "blue", plot.OOB = TRUE)
if (pdfPlot) dev.off()

set.seed(100)
if (pdfPlot) pdf("Bagging2.pdf", width = 10, height = 6)
plot(x, y, cex = cex.points)
lines(xnew, fun(xnew), lty="dotted")
treeBoot(x, y, control = tree.control(nobs = length(x), mincut = 1, minsize = 2), col = "violet")
if (pdfPlot) dev.off()
set.seed(100)
if (pdfPlot) pdf("Bagging2OOB.pdf", width = 10, height = 6)
plot(x, y, cex = cex.points)
lines(xnew, fun(xnew), lty="dotted")
treeBoot(x, y, control = tree.control(nobs = length(x), mincut = 1, minsize = 2), col = "violet",
         plot.OOB = TRUE)
if (pdfPlot) dev.off()

# bagging trees
B <- 500
pBoot <- matrix(NA, length(xnew), B)

#par(mfrow = c(2,1))

if (pdfPlot) pdf("Bagging500.pdf", width = 10, height = 6)
plot(x, y)
for (b in 1:B){
  pBoot[, b] <- treeBoot(x, y, plot = TRUE, col = "lightgrey",
                         control = tree.control(nobs = length(x), 
                                                mincut = 1, minsize = 2),
                         size.adapt = FALSE,
                         lty = "dotted", lwd = 0.2)
}
pBagging <- rowMeans(pBoot)
lines(xnew, pBagging, lwd = 3, col = "black")
lines(xnew, fun(xnew), lty="dotted", col = "blue")
if (pdfPlot) dev.off()
#lines(xnew, fun(xnew), lwd = 3, col = "red", lty = "dashed")



library(randomForest)  # of course here, feature sampling is useless
d <- data.frame(x = x, y = y)
mRF <- randomForest(y~x, data = d, keep.inbag = TRUE)
plot(x, y)
pRF <- predict(mRF, data.frame(x = xnew))
lines(xnew, pRF, lwd = 3, col = "blue")
lines(xnew, fun(xnew), lwd = 3, col = "red", lty = "dashed")

