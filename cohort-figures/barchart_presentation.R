
v <- data.frame(taxa=c('A','B','C'), abund=c(6, 6, 3))

ggplot(data = v, aes(x=taxa, y=abund, fill=taxa)) + geom_bar(stat="identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("#a1d76a", "#56B4E9", "orange"))
  