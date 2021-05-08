library(eulerr)

# family updated
VennDiag1 <- euler(c("16S" = 337, "MGX" = 55,"16S&MGX" = 17))
plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

# genera updated
VennDiag2 <- euler(c("16S" = 197, "MGX" = 40,"16S&MGX" = 158))
plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S", "MGX")))

# species updated
VennDiag3 <- euler(c("16S" = 127, "MGX" = 383,"16S&MGX" = 128))
plot(VennDiag3, counts = TRUE, font=1, cex=1, alpha=0.5,
     fill=c("coral", "cyan", "lightgrey"),
     quantities = list(cex = 5),
     legend = list(labels = c("16S","MGX")))

