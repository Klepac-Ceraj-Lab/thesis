library(ggplot2)
library(vegan)

# looking at shannnon data from cohort

setwd("/Users/danielle/Documents/thesis/theoretical")
babies <- read.csv("/Users/danielle/Documents/thesis/theoretical/shannon.csv")

# creating age categories

babies$AgeMonths <- babies$correctedAgeDays/30.0

babies$dev_stage[babies$AgeMonths <= 15] <- "less than 15 months"
babies$dev_stage[babies$AgeMonths > 15 & babies$AgeMonths <= 30] <- "15 to 30 months"
babies$dev_stage[babies$AgeMonths > 30] <- "older than 30 months"

babies$color[babies$dev_stage=="less than 15 months"] <- "red"
babies$color[babies$dev_stage=="15 to 30 months"] <- "blue"
babies$color[babies$dev_stage=="older than 30 months"] <- "yellow"

write.csv(babies, "sorted_babies.csv" )


# plotting correlation age vs Shannon diversity

shannon1 <- plot(log(babies$AgeMonths),babies$shannon,
     xlab = "Log Age(months)",
     ylab = "Shannon diversity",
     col = babies$color)

dev_stages <- c("less than 15 months", "15 to 30 months", "older than 30 months")
colors <- c("red", "blue", "yellow")

legend("bottomright", legend = dev_stages, fill = colors, title = 'Developmental Stages')

# boxplots of Shannon diversity given developmental stage
babies$dev_stage <- factor(babies$dev_stage, 
    levels=c("less than 15 months", "15 to 30 months", "older than 30 months"))
shannon2 <- boxplot(babies$shannon~babies$dev_stage, 
        xlab = "developmental stage",
        ylab = "Shannon diversity")

par(mfrow=c(1,2))
shannon1
shannon2

# calculating mean of shannon diversity for different developmental stages

mean(babies[babies$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
mean(babies[babies$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
mean(babies[babies$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)

IQR(babies[babies$dev_stage == "less than 15 months",]$shannon, na.rm=TRUE)
IQR(babies[babies$dev_stage == "15 to 30 months",]$shannon, na.rm=TRUE)
IQR(babies[babies$dev_stage == "older than 30 months",]$shannon, na.rm=TRUE)
