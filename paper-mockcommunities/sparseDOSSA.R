library(BiocManager)
library("sparseDOSSA")
library("dplyr")
library("vegan")

n.microbes <- 150
n.samples <- 2500

#sparseDOSSA::sparseDOSSA( number_features = n.microbes, 
                          #number_samples = n.samples)

setwd("/Users/danielle/Documents/thesis/paper-mockcommunities")

# read in file of mock communities
mock <- read.csv("SyntheticMicrobiome-Counts.csv", header=TRUE, 
                 stringsAsFactors=FALSE)

rownames(mock) <- mock$X.SampleID
mock$X.SampleID <- NULL

# only keep rows starting with Feature_Lognormal
mock <- mock[grep('^Feature_Lognormal', rownames(mock)),]

# tranpose
mock <- t(mock)

# convert to numeric
mock <- as.data.frame(   
  apply(mock, 2, as.numeric))
sapply(mock, class)

# calculate minimum read length
mock["total"] <- rowSums(mock, na.rm = TRUE)
mock["shannon"] <- diversity(mock, index="shannon")
mock["evenness"] <- diversity(mock, "simpson")
mock["richness"] <- apply(mock>0,1,sum)
mock["least_abund_bug"] <- apply(mock[,1:150], 1, FUN = function(x) {min(x[x > 0])})
mock["min_abund"] <- mock["least_abund_bug"]/mock["total"]


genome_length <- 4000000 # assuming that average genome length is 4 mil bp
read_length <- 150 # assuming average read length to be 150 bp

mock["reads"] <- genome_length/(read_length*mock$min_abund)

write.csv(mock, "mock_communities.csv")


plot<- ggplot(mock, aes(log(richness), evenness))
plot1 <- plot + geom_point(aes(color = reads), size = 1.0) +
  scale_color_gradientn(colours = rainbow(6), labels= comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position=c(0.85, 0.5))+
  labs(title = "", tag = "A")+
  stat_compare_means(method = "t.test", label = "p.signif", label.x = 2, label.y = 1)

+xlim(3.25,4.25)

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

babies2$dev_stage <- factor(babies2$dev_stage, levels = c("less than 15 months", 
                                                          "15 to 30 months", 
                                                          "older than 30 months"))


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
babies3 <- subset(babies2, select=c("sampleid", "dev_stage", "read_predictions"))
colnames(babies3)[3] <- "predicted necessary depth"
babies3 <- melt(babies3)

babies4 <- subset(babies2, select=c("sampleid", "dev_stage", "read_depth"))
colnames(babies4)[3] <- "sequenced depth"
babies4 <- melt(babies4)

babies5 <- rbind(babies3, babies4)


p3 <- ggplot(na.omit(babies5), aes(x = dev_stage, y = value)) + 
  geom_boxplot(aes(colour = variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("read depth") + xlab("developmental stage")+
  labs(title = "", tag = "C")+ theme(legend.title = element_blank()) +
  theme(legend.position = c(0.25, 0.85)) +
  scale_y_continuous(label=comma)+
  scale_fill_manual(values = c("#E7B800", "#FC4E07"))
p3

# statistics
aov <- aov(read_predictions~dev_stage, babies2)
summary(aov)
TukeyHSD(aov)

# plotting read predictions

babies2$dev_stage <- factor(babies2$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)

plot2 <- ggplot(data=subset(babies2, !is.na(dev_stage)), aes(log(richness), 
                                                             evenness))
plot2<- plot2+ geom_point(aes(color = dev_stage), alpha=0.4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = 
          element_line(colour = "black"),
        legend.position=c(0.25, 0.85)) +
  guides(color=guide_legend("age")) + ylim(0, 1.0) + xlim(0,5)+
  labs(title = "", tag = "B")

gl <- list(plot1, p3, plot2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 2,2))
)  

