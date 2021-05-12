library(eulerr)

# family updated
VennDiag1 <- euler(c("16S" = 223, "MGX" = 4,"16S&MGX" = 66))
plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

# genera updated
VennDiag2 <- euler(c("16S" = 375, "MGX" = 40,"16S&MGX" = 154))
plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

# species
VennDiag3 <- euler(c("16S" = 146, "MGX" = 317,"16S&MGX" = 150))
plot(VennDiag3, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S","MGX")))
