
# making legend for itol plot

legend_df <- data.frame("method"= c("both", "16S", "mgx"), "color"= c("#9ebcda", "#fa9fb5", "#7fcdbb"))

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('both', '16S', 'mgx'), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('#9ebcda', '#fa9fb5', '#7fcdbb'))
mtext("profiling method", at=0.2, cex=2)
