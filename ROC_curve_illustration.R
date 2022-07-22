set.seed(6)
n <- 15
x1 <- rnorm(n = n, mean = 5, sd = 1)
x2 <- rnorm(n = n, mean = 7, sd = 1)
x <- c(x1, x2)
y <- factor(c(rep("blue", n), rep("red", n)), levels = c("blue", "red"))

m <- glm(y ~ ., data = data.frame(x = x, y = y), family = "binomial")
summary(m)

pdfPlot <- FALSE

myplot <- function(m, s = 0.5, ROC = FALSE){
  x <- m$data$x
  y <- m$data$y
  col <- ifelse(y == "red", "red", "blue")
  pch <- ifelse(y == "red", 20, 4)
  cex <- 2
  xmin <- min(x) * 0.9
  xmax <- max(x) * 1.1
  
  #par(mfrow = c(1, ifelse(ROC, 2, 1)))
  plot(x, y == "red", col = col, pch = pch, cex = cex, 
       xlab = "x", xlim = c(xmin, xmax),
       ylab = "Probability", ylim = c(0, 1))  
  abline(h = c(0, 1), col = "black")
  prob <- predict(m, type = "response")
  t <- seq(xmin * 0.8, xmax * 1.2, length.out = 200)
  xdotb <- coef(m)[1] + t*coef(m)[2]
  lines(t, exp(xdotb)/(1 + exp(xdotb)), lwd = 2)
  points(x, prob, col = col, pch = pch, cex = cex)
  abline(h = s, lty = "dotted")
  p <- ifelse(prob > s, "red", "blue")
  
  if (ROC){
    tab <- table(p, y)
    TP <- tab["red", "red"] / sum(tab[, "red"])
    FP <- tab["red", "blue"] / sum(tab[, "blue"])
    plot(FP, TP, col = "brown", cex = cex, pch = 19, 
         xlab = "FALSE positive rate", xlim = c(0, 1), 
         ylab = "TRUE positive rate", ylim = c(0, 1))
    text(FP, TP, paste("(s = ", s, ")"), 
         pos = ifelse(TP > 0.8, 2, 4), offset = 1)
  }
  
  return(p = p)
}

if (pdfPlot) pdf("regLog_n30.pdf", width = 10, height = 5)
myplot(m, s = 0.5, ROC = FALSE)
if (pdfPlot) dev.off()

for (s in 0:9){
  if (pdfPlot) pdf(paste("regLog_ROC_s_", s, ".pdf", sep = ""), width = 10, height = 5)
  par(mfrow = c(1, 2))
  myplot(m, s = s/10, ROC = TRUE)
  if (pdfPlot) dev.off()
}

library(ROCR)
roclogit <- predict(m, type="response")
predlogit <- prediction(roclogit, y)
perflogit <- performance(predlogit, "tpr","fpr")
# Tracé de la courbe
if (pdfPlot) pdf("regLog_ROC_fullCurve.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
myplot(m, s = 1)
plot(perflogit, col = "brown", type = "o", pch = 19, cex = 2, lwd = 2)
if (pdfPlot) dev.off()
