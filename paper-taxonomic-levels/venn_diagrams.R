library(eulerr)
set.seed(1)

family <- c(A = 1, B = 2)

plot(venn(family))
plot(euler(family), quantities = TRUE)


v <- venneuler(c(A=450, B=1800, "A&B"=230))


library(VennDiagram) 
venn.diagram(list(B = 1:74, A = 41:88), fill = c("pink", "blue"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_diagram.tiff")
