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

plot1 <- plot + geom_point(aes(color = reads), size = 1.0) +
  scale_color_gradientn(colours = rainbow(6), labels= comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.85, 0.25))

plot1

# adding on actual data points

babies2 <- read.csv("/Users/danielle/Documents/thesis/theoretical/theoretical_babies_df.csv")

babies2["shannon"] <- diversity(babies2[,6:134], index="shannon")
babies2["evenness"] <- diversity(babies2[,6:134], "simpson")
babies2["richness"] <- apply(babies2[,6:134]>0,1,sum)

# creating linear model
log_reads <- log(mock$reads)
log_reads[is.infinite(log_reads)] <- NA
model <- lm(log_reads ~ evenness + log(richness), data=mock)
predict_data <-  babies2
predictions <- predict(model, predict_data)
babies2$read_predictions <- exp(predictions) #unlog transformation

babies2$dev_stage <- factor(babies2$dev_stage, levels = c("less than 15 months", "15 to 30 months", "older than 30 months"))


# calculating mean sequencing depth for each age group

mean.under15 <- mean(na.omit(babies2[babies2$dev_stage == "less than 15 months",]$read_predictions))
var.under15 <- var(na.omit(babies2[babies2$dev_stage == "less than 15 months",]$read_predictions))
c(mean.under15-1.96*sqrt(var.under15), mean.under15+1.96*sqrt(var.under15))

mean.over15 <- mean(na.omit(babies2[babies2$dev_stage == "15 to 30 months" 
                |babies2$dev_stage == "older than 30 months",]$read_predictions))
var.over15 <- var(na.omit(babies2[babies2$dev_stage == "15 to 30 months" 
                                 |babies2$dev_stage == "older than 30 months",]$read_predictions))
c(mean.over15-1.96*sqrt(var.over15), mean.over15+1.96*sqrt(var.over15))

mean(babies2[babies2$dev_stage == "older than 30 months",]$read_predictions)

# boxplots
par(mfrow=c(1,2))
boxplot(read_predictions~dev_stage, babies2)
boxplot(read_depth~dev_stage, babies2)

# statistics
aov <- aov(read_predictions~dev_stage, babies2)
summary(aov)
TukeyHSD(aov)

# plotting read predictions

plot2 <- ggplot(data=subset(babies2, !is.na(dev_stage)), aes(log(richness), evenness))
plot2<- plot2+ geom_point(aes(color = dev_stage), alpha=0.4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.15, 0.85)) +
  guides(color=guide_legend("age")) + ylim(0, 1.0) + xlim(0,4)
  
plot2
  
grid.arrange(plot1, plot2, nrow=2)

