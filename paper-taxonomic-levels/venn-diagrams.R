library(eulerr)

VennDiag1 <- euler(c("16S" = 33, "MGX" = 14,"16S&MGX" = 41))
plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

VennDiag2 <- euler(c("16S" = 62, "MGX" = 34,"16S&MGX" = 105))
plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

VennDiag3 <- euler(c("16S" = 0, "MGX" = 385,"16S&MGX" = 0))
plot(VennDiag3, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("MGX")))
