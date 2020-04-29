library(vegan)
library(ggplot2)
library(gridExtra)
library(scales)

set.seed(3434)

setwd("/Users/danielle/Documents/thesis/theoretical")

mock <- read.csv("mock_communities2.csv", header=TRUE)

for (ii in 1:30) {
  
    mock[ii, 2:(mock[ii, "num_bugs"]+1)] <- (rpois(mock[ii,"num_bugs"], mock[ii,"parameter"]))
    
}

mock["total"] <- rowSums(mock[,2:101], na.rm = TRUE)
mock[is.na(mock)] <- 0
mock["shannon"] <- diversity(mock[,2:101], index="shannon")
mock["evenness"] <- diversity(mock[,2:101], "simpson")
mock["richness"] <- apply(mock[,2:101]>0,1,sum)
mock["min_abund"] <- (1.0/mock["total"])


genome_length <- 4000000
read_length <- 150

mock["reads"] <- genome_length/(read_length*mock$min_abund)

write.csv(mock, "mock_communities2_filled.csv" )

##### Make plot of estimated reads ~ eveness + richness ######

plot<- ggplot(mock, aes(log(richness), evenness))

plot1 <- plot + geom_point(aes(color = reads)) +
  scale_color_gradientn(colours = rainbow(6), labels= comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.1, 0.85)) + 
  geom_smooth(method='lm', formula= y~x)

plot1

# adding on actual data points

babies2 <- read.csv("/Users/danielle/Documents/thesis/theoretical/theoretical_babies_df.csv")

babies2["shannon"] <- diversity(babies2[,6:134], index="shannon")
babies2["evenness"] <- diversity(babies2[,6:134], "simpson")
babies2["richness"] <- apply(babies2[,6:134]>0,1,sum)

# creating linear model
log_reads <- log(mock$reads)
log_reads[is.infinite(log_reads)] <- NA
model <- lm(log_reads ~ evenness + richness, data=mock)
predict_data <-  babies2
predictions <- predict(model, predict_data)
babies2$read_predictions <- exp(predictions) #unlog transformation

plot2 <- plot + geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Shannon Index") + ylab("estimated reads")

babies2$dev_stage <- factor(babies2$dev_stage, levels = c("less than 15 months", "15 to 30 months", "older than 30 months"))


# calculating mean sequencing depth for each age group

plot3 <- ggplot(data=subset(babies2, !is.na(dev_stage)), aes(log(richness), evenness)) + 
  geom_point(aes(color = dev_stage), alpha=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.15, 0.55))+
  guides(color=guide_legend("age"))

plot3
  
grid.arrange(plot1, plot3, nrow=2)

