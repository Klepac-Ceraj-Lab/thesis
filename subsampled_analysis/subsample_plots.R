library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(phyloseq)
library(ggpubr)

# index of where bug abundances columns start and end
start_bugs <- 8
end_bugs <- 107

setwd("/Users/danielle/Documents/thesis/subsampled_analysis")

df <- read.csv("subsampled_df.csv", header=TRUE)


df[is.na(df)] <- 0
df["shannon"] <- diversity(df[,start_bugs:end_bugs], "shannon")
df["evenness"] <- diversity(df[,start_bugs:end_bugs], "simpson")
df["richness"] <- apply(df[,start_bugs:end_bugs]>0,1,sum)

# subsample children statistics

# average age of children less than 15 months old
mean(df[df$dev_stage == "less than 15 months",]$AgeMonths)

# number of kids in each age group
nrow(df[df$dev_stage == "less than 15 months",])
nrow(df[df$dev_stage == "15 to 30 months",])
nrow(df[df$dev_stage == "older than 30 months",])

# plots evenness and richness by developmental stage and read depth

# color by devlopmental stage
# opacity by read depth
# darker colors are more reads

# remove 10k samples

df <- subset(df, sampling_cat != 10)

plot1 <- ggplot(df, aes(sampling_cat, shannon)) + 
  geom_point(aes(shape=dev_stage, color = sampling_cat)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x="sampling depth", y="Shannon Index", 
       color="read depth", shape="developmental stage") +
  labs(title = "", tag = "A") 
plot1

df$dev_stage <- factor(df$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)
df$sampling_cat <- factor(df$sampling_cat,
                          levels = c(10, 100, 250, 500, 750, 1000, "original depth"),
                          ordered = TRUE)

plot2 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(color=dev_stage, alpha = sampling_cat))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", 
       color="read depth", alpha="sampling depth") +
  labs(title = "", tag = "A")+ theme(legend.position = c(.85,.35)) +
  ylab("evenness (Pielou's measure)")
plot2

plot3 <- ggplot(df, aes(richness, evenness)) + 
  geom_point(aes(shape=sampling_cat, color = sampleid)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="richness", y="evenness", color="sample", shape = "sampling depth")

plot3


gl <- list(plot1, plot2, plot3)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1,2, 2),
                        c(3, 2,2))
)

# pcoa plots
# http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

# put RDA plots in supplement
# plot C, stratified by developmental stage
# grouped by developmental stage, sampling depth for each developmental stage
# tests: compare between read depths for all samples (what we have currently)
# in addition: compare read depth within each age group
# with younger kids, can you get away with getting fewer reads??

abund_table <- df[,start_bugs:end_bugs]
abund_table<-subset(abund_table,rowSums(abund_table)!=0)
meta_table <- subset(df,rowSums(df[,start_bugs:end_bugs])!=0)[,1:start_bugs-1]

sol<-rda(abund_table ~ ., data=meta_table)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
df_sites<-data.frame(scrs$sites,meta_table$sampleid, 
                     meta_table$sampling_cat, meta_table$dev_stage)
colnames(df_sites)<-c("x,","y","sampleid", "sampling_cat", "dev_stage")

# finding main axes
axis.expl <- function(mod, axes = 1:2) {
  
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
  
}

axises <- axis.expl(sol)
axis_1 <- axises[1]
axis_2 <- axises[2]

# coloring by dev stage and sampling depth



df_sites$dev_stage <- factor(df_sites$dev_stage,
                            levels = c("less than 15 months", 
                                       "15 to 30 months", 
                                       "older than 30 months"),ordered = TRUE)


plot4<-ggplot(df_sites, aes(df_sites$x, df_sites$y, colour=sampling_cat, shape = dev_stage)) +
  geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="RDA 1, 37.01%", y="RDA 2, 21.67%", 
                  color="sampling depth", shape = "developmental stage")+ 
  theme(legend.position = "none")+
  labs(title = "", tag = "B")
plot4

# calculating statistics

df_subsampled <- subset(df, sampling_cat != "original depth")

plot5 <- ggplot(df_subsampled, aes(x = sampling_cat, y = shannon, fill= dev_stage)) + 
  geom_boxplot()
plot5 <- plot5 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("alpha-diversity (Shannon Index)")+ xlab("sampling depth (k reads)") +
  labs(title = "", tag = "B") + theme(legend.position = "none") + ylim(0,3)
plot5

original_data <- subset(df, sampling_cat == "original depth")
original_data$dev_stage <- factor(original_data$dev_stage,
                             levels = c("less than 15 months", 
                                        "15 to 30 months", 
                                        "older than 30 months"),ordered = TRUE)

plot6 <- ggplot(original_data, aes(x = (read_depth/1000), y = shannon)) +  
  geom_point(aes(color = dev_stage))
plot6 <- plot6 +
  theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line.x = element_line(colour = "black"),
                 axis.text.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(), 
                axis.line.y = element_blank()) + 
  theme(legend.position = "none") + 
  ylim(0, 3) +
  xlab("original read depth")
plot6

original_data$read_depth <- original_data$read_depth/1000
alldata <- rbind(original_data, df_subsampled)


ggplot(alldata, aes(x = (read_depth), y = shannon)) +  
  geom_point(aes(color = dev_stage)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank()) + 
  theme(legend.position = "none") + 
  ylim(0, 3) +
  xlab("original read depth")


gl <- list(plot2, plot5, plot6)

lay <- rbind(c(1,1,2,2,3),
             c(1,1,2,2,3),
             c(1,1,2,2,3))

grid.arrange(grobs = gl,layout_matrix = lay)

mean((df[df$sampling_cat == "100" ,]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "100" ,]$shannon), na.rm = TRUE)

mean((df[df$sampling_cat == "250" ,]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "250" ,]$shannon), na.rm = TRUE)

mean((df[df$sampling_cat == "500" ,]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "500" ,]$shannon), na.rm = TRUE)

mean((df[df$sampling_cat == "750" ,]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "750" ,]$shannon), na.rm = TRUE)

mean((df[df$sampling_cat == "1000" ,]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "1000" ,]$shannon), na.rm = TRUE)

mean((df[df$sampling_cat == "original depth",]$shannon), na.rm = TRUE)
IQR((df[df$sampling_cat == "original depth",]$shannon), na.rm = TRUE)

anova_sampling <- aov(df$shannon~df$sampling_cat)
summary(anova_sampling)

anova <- aov(df$shannon~df$sampling_cat*df$dev_stage)
summary(anova)
posthoc <- TukeyHSD(anova,conf.level=0.95)
posthoc <- as.data.frame(posthoc$`df$sampling_cat:df$dev_stage`)

keep_cat_tukey <- scan("/Users/danielle/Documents/thesis/subsampled_analysis/keep_cat_tukey.txt", 
     what="character", sep=",", 
     strip.white = TRUE)

# only keep comparisons we care about and highlight significant p-values

posthoc$comparisons <- rownames(posthoc)
posthoc <- posthoc[posthoc$comparisons %in% keep_cat_tukey,]
posthoc$less_0.05 <- posthoc$`p adj`< 0.05
posthoc$less_0.005 <- posthoc$`p adj`< 0.005

# write.csv(posthoc,"subsampling_tukey.csv")

# range of shannon diversity by age

max((df[df$dev_stage == "less than 15 months" ,]$shannon), na.rm = TRUE) - min((df[df$dev_stage == "less than 15 months" ,]$shannon), na.rm = TRUE)
median((df[df$dev_stage == "less than 15 months" ,]$shannon), na.rm = TRUE)
max((df[(df$dev_stage == "15 to 30 months"  | df$dev_stage == "older than 30 months"),]$shannon), na.rm = TRUE) - min((df[(df$dev_stage == "15 to 30 months"  | df$dev_stage == "older than 30 months"),]$shannon), na.rm = TRUE)
median((df[(df$dev_stage == "15 to 30 months"  | df$dev_stage == "older than 30 months"),]$shannon), na.rm = TRUE)

par(mfrow=c(2,1))
hist((df[df$dev_stage == "less than 15 months" ,]$shannon), 
     na.rm = TRUE,
     xlim = c(0,3), xlab = "Shannon Index", main = "A. Children younger than 15 months")
hist((df[(df$dev_stage == "15 to 30 months"  | df$dev_stage == "older than 30 months"),]$shannon), 
     na.rm = TRUE,
     xlim = c(0,3), xlab = "Shannon Index", main = "B. Children older than 15 months")


# supplemental figures

# broken down by dev_stage

update_geom_defaults("point", list(colour = NULL))
m <- ggplot(movies, aes(y = votes, x = factor(round(rating)),
                        colour = factor(Animation)))

my_comparisons_1 <- list(c("100", "250"), c("250", "750"), c("750", "1000"))

# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
# define outlier function, beyound 5 and 95% percentiles
o <- function(x) {
  subset(x, x < quantile(x,probs=c(0.05))[1] | quantile(x,probs=c(0.95))[1] < x)
}

s1 <- ggplot(df, aes(x = sampling_cat, y = richness, fill= dev_stage)) + 
  facet_wrap(~dev_stage) +
  geom_boxplot()
s1 <- s1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("richness (genus counts)")+ xlab("sampling depth (k reads)")  + 
  labs(title = "", tag = "A")  

s1 <- s1 + stat_summary(fun.data=f, geom='boxplot')
s1 <- s1 + stat_summary(fun.y=o, geom='point', aes(colour=factor(dev_stage)))
s1 + scale_y_log10()
s1

labs(title = "", tag = "B") + theme(legend.position = "none")+ ylim(0,40)

# broken down by sampling depth

update_geom_defaults("point", list(colour = NULL))
s2 <- ggplot(df, aes(x = sampling_cat, y = richness, fill= dev_stage)) + 
  geom_boxplot(outlier.colour = NULL)
s2 <- s2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black")) +
  ylab("richness (genus counts)") + xlab("sampling depth (k reads)")  + 
  labs(title = "", tag = "B") 
s2

