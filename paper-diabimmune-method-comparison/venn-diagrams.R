library(eulerr)

VennDiag1 <- euler(c("16S" = 16, "MGX" = 51,"16S&MGX" = 39))
plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

VennDiag2 <- euler(c("16S" = 14, "MGX" = 81,"16S&MGX" = 77))
plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

VennDiag3 <- euler(c("16S" = 0, "MGX" = 519,"16S&MGX" = 0))
plot(VennDiag3, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("MGX")))

