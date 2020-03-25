library(vegan)
library(ggplot2)
library(gridExtra)

set.seed(3434)

setwd("/Users/danielle/Documents/thesis/theoretical")

mock <- read.csv("mock_communities2.csv", header=TRUE)

for (ii in 1:30) {
  
    mock[ii, 2:(mock[ii, "num_bugs"]+1)] <- (rpois(mock[ii,"num_bugs"], mock[ii,"parameter"]))
    
}

mock["total"] <- rowSums(mock[,2:101], na.rm = TRUE)
mock[is.na(mock)] <- 0
mock["shannon"] <- diversity(mock[,2:101], index="shannon")
mock["min_abund"] <- (1.0/mock["total"])


genome_length <- 4000000
read_length <- 150

mock["reads"] <- genome_length/(read_length*mock$min_abund)

write.csv(mock, "mock_communities2_filled.csv" )


plot<- ggplot(mock, aes(shannon, reads)) 


plot1 <- plot + geom_point() + 
  geom_smooth(aes(color =as.factor(parameter)), se=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.22, 0.55)) +
  geom_jitter() +
  guides(color=guide_legend("sampling distribution parameters", ncol=2))+
  labs(color=mock$parameter)+
  xlab("Shannon Index") + ylab("estimated reads")

# adding on actual data points

babies2 <- read.csv("/Users/danielle/Documents/thesis/theoretical/sorted_babies.csv")

log_reads <- log(mock$reads)
log_reads[is.infinite(log_reads)] <- NA

model <- lm(log_reads ~ mock$shannon)
exp(predict(model, babies2))
babies2$read_predictions <- exp(babies2$shannon * 1.179 + 10.322) #unlog transformation

  
plot2 <- plot + geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Shannon Index") + ylab("estimated reads")

babies2$dev_stage <- factor(babies2$dev_stage, levels = c("less than 15 months", "15 to 30 months", "older than 30 months"))


# calculating mean sequencing depth for each age group

babies3 <- babies2[!is.na(babies2$dev_stage),] # remove nas

mean_less15 <- mean(babies3[babies3$dev_stage == "less than 15 months","read_predictions"])
mean_15to30 <- mean(babies3[babies3$dev_stage == "15 to 30 months","read_predictions"])
mean_greater30 <- mean(babies3[babies3$dev_stage == "older than 30 months","read_predictions"])

plot3 <- ggplot(data=subset(babies2, !is.na(dev_stage)), aes(shannon, read_predictions)) + 
  geom_point(aes(color = dev_stage), alpha=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.15, 0.55))+
  guides(color=guide_legend("age")) +
  xlab("Shannon Index") + ylab("estimated reads")+
  geom_hline(yintercept = mean_less15, color="red", alpha=0.4)+
  geom_hline(yintercept = mean_15to30, color="green", alpha=0.4)+
  geom_hline(yintercept = mean_greater30, color="blue", alpha=0.4)
  
  
grid.arrange(plot1, plot3, nrow=2)

