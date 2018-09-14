#libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(survival)
library(survminer)
library(survMisc)
library(pROC)
theme_set(theme_pubr())
library(ber)
library(rms)
library(readr)
library(readxl)
library(ROCR)


#read data file
df <- read.csv("D:/Study/BCL2/df_total_180121.csv", 
               stringsAsFactors = FALSE)
df$phenotype[which(df$phenotype == "CD20 -")] <- "CD20-"
df$phenotype[which(df$phenotype == "CD20+ Mb")] <- "CD20+"
df$phenotype[which(df$phenotype == "CD20- Mb")] <- "CD20-"
df.bcl2 <- dplyr::select(df, c(phenotype, x_position, y_position, nuc_DAPI_mean, nuc_myc_mean, cyt_bcl2_mean, slideID, rowTMA, colTMA, confidence, SN, batch))
df.bcl2 <- dplyr::filter(df.bcl2, confidence > 79.16)
df.bcl2$cyt_bcl2_mean[which(df.bcl2$cyt_bcl2_mean == "#N/A")] <- 0
df.bcl2$cyt_bcl2_mean <- as.numeric(df.bcl2$cyt_bcl2_mean)
df.cd20 <- dplyr::filter(df.bcl2, phenotype == "CD20+")
df.cd3 <- dplyr::filter(df.bcl2, phenotype == "CD3+")

#remove batch effect of BCL2 in tumor cells
nl.2 <- as.matrix(df.cd20$cyt_bcl2_mean)
batch2 <- as.factor(df.cd20$batch)
nl.2 <- as.data.frame(ber(nl.2, batch2))
nl.2 <- dplyr::rename(nl.2, bcl2 = V1)
df.cd20 <- cbind(df.cd20, nl.2)

# remove batch effect of BCL2 in CD3
nl.3 <- as.matrix(df.cd3$cyt_bcl2_mean)
batch3 <- as.factor(df.cd3$batch)
nl.3 <- as.data.frame(ber(nl.3, batch3))
nl.3 <- dplyr::rename(nl.3, bcl2 = V1)
df.cd3 <- cbind(df.cd3, nl.3)

# remove batch effect of MYC in tumor cells
nl.1 <- as.matrix(df.cd20$nuc_myc_mean)
batch1 <- as.factor(df.cd20$batch)
nl.1 <- as.data.frame(ber(nl.1, batch1))
nl.1 <- dplyr::rename(nl.1, myc = V1)
df.cd20 <- cbind(df.cd20, nl.1)

# normalize BCL2 and MYC data
df.cd20$nor.bcl2 <- 
        (df.cd20$bcl2 - min(df.cd20$bcl2))/(max(df.cd20$bcl2) - min(df.cd20$bcl2))

df.cd3$nor.bcl2 <- 
        (df.cd3$bcl2 - min(df.cd3$bcl2))/(max(df.cd3$bcl2) - min(df.cd3$bcl2))

df.cd20$nor.myc <- 
        (df.cd20$myc - min(df.cd20$myc)) / (max(df.cd20$myc) - min(df.cd20$myc))

###################################################################
### Conventional analysis using proportional cutoff ###############
###################################################################

# open medical data
med.07_09 <- read_excel("D:/Study/MYC analysis/clinical07_09_2.xlsx", sheet = 1, col_names = TRUE)
med.10_12 <- read_excel("D:/Study/MYC analysis/clinical10_12_2.xlsx", sheet = 1, col_names = TRUE)

# cleansing the clinical data
med.data <- rbind(med.07_09, med.10_12)
med.data$Site <- gsub("Stomcah", "Stomach", med.data$Site)
med.data <- dplyr::filter(med.data, Site != "Brain")
#med.data <- dplyr::filter(med.data, Site != "Stomach")
med.data <- dplyr::filter(med.data, Tx.regimen %in% c("R-CHOP", "R-CVP"))
med.data[which(duplicated(med.data$ID) == TRUE),]

med.data$OS[which(med.data$OS == 9)] <- 0
med.data$Event[which(med.data$Event == 9)] <- 0
med.data$yr5.OS[which(med.data$yr5.OS == 9)] <- 0
med.data$yr2.Event[which(med.data$yr2.Event == 9)] <- 0

# OS duration / EFS duration
med.data$Dx.Date <- ymd(med.data$Dx.Date)
med.data$yr5.Last.Date <- ymd(med.data$yr5.Last.Date)
med.data$yr2.Event.Date <- ymd(med.data$yr2.Event.Date)
med.data$Event.Date <- ymd(med.data$Event.Date)
med.data$Last.Date <- ymd(med.data$Last.Date)
med.data$CR.Date <- ymd(med.data$CR.Date)

mon.os <- as.period(interval(med.data$Dx.Date, med.data$Last.Date), unit = "months")
med.data$mon.os <- mon.os%/% months(1)
mon.efs <- as.period(interval(med.data$Dx.Date, med.data$Event.Date), unit = "months")
med.data$mon.efs <- mon.efs%/% months(1)
mon.5yr.os <- as.period(interval(med.data$Dx.Date, med.data$yr5.Last.Date), unit = "months")
med.data$mon.5yr.os <- mon.5yr.os %/% months(1)
mon.2yr.efs <- as.period(interval(med.data$Dx.Date, med.data$yr2.Event.Date), unit = "months")
med.data$mon.2yr.efs <- mon.2yr.efs %/% months(1)

# Distribution of BCL2 expression on tumor cells
summary(df.cd20$nor.bcl2)
# Min 0.00000 / 1Q 0.07177 / Median 0.09441 / 3Q 0.13310 / Max 1.00000

# Histogram of BCL2 expression on tumor cells
options(scipen = 999)
g <- qplot(nor.bcl2, data = df.cd20, geom = "histogram", color = "black",
           bins = 100)

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/bcl2_tumor.jpg")
print(g)
dev.off()

# Distribution of BCL2 expression on CD3+ T-cells
summary(df.cd3$nor.bcl2)
# Min 0.0000 / 1Q 0.1373 / Median 0.1926 / 3Q 0.2638 / Max 1.0000

# Histogram of BCL2 expression on T-cells
g <- qplot(nor.bcl2, data = df.cd3, geom = "histogram", color = "black",
           bins = 100)
print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/bcl2_Tcell.jpg")
print(g)
dev.off()

# Using intensity cutpoint (0.1)
df.cd20$Q1 <- ifelse(df.cd20$nor.bcl2 < 0.1, 0, 1)
Q1 <- aggregate(Q1 ~ SN, data = df.cd20, FUN = sum)
tot <- aggregate(Q1 ~ SN, data = df.cd20, FUN = length)
tot <- dplyr::rename(tot, total = Q1)
pos <- inner_join(Q1, tot, by = "SN")
pos$Q1.ratio <- (pos$Q1 / pos$total) * 100
pos <- dplyr::select(pos, c(1,4))

# merge cell data and medical data
df.amc <- inner_join(pos, df.amc, by = "SN")

# cox proportional hazard model for Q1 ratio
fit <- coxph(Surv(mon.os, OS == 1) ~ Q1.ratio, data = df.amc)
summary(fit) # p value = 0.00116 / exp(coef) = 1.012517

# 50% proportional cutpoint
df.amc$cut50 <- ifelse(df.amc$Q1.ratio < 50, 0, 1)
df.amc$cut50 <- factor(df.amc$cut50, levels = c(0,1), 
                       labels = c("Negative", "Positive"))

# K-M curve (OS)
fit <- survfit(Surv(mon.os, OS == 1) ~ cut50, data = df.amc)
g <- ggsurvplot(fit,
                surv.plot.height = 0.8,
                pval = TRUE, pval.size = 4,
                xlab = "Time (month)", font.x = 12,
                ylab = "OS", font.y = 12,
                risk.table = TRUE, risk.table.height = 0.35,
                risk.table.fontsize = 3.5,
                risk.table.y.text.col = FALSE,
                legend = "right",
                legend.tital = "BCL2",
                legend.labs = c("Negative", "Positive"),
                palette = c("blue", "red"))

g$table <- g$table +
        theme(plot.title = element_text(hjust = 0))

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/blc2_os_conventional.jpg")
print(g)
dev.off()

# K-M curve (2yr EFS)
fit <- survfit(Surv(mon.2yr.efs, yr2.Event == 1) ~ cut50, 
               data = df.amc)
g <- ggsurvplot(fit,
                surv.plot.height = 0.8,
                pval = TRUE, pval.size = 4,
                xlab = "Time (month)", font.x = 12,
                ylab = "2yr EFS", font.y = 12,
                risk.table = TRUE, risk.table.height = 0.35,
                risk.table.fontsize = 3.5,
                risk.table.y.text.col = FALSE,
                legend = "right",
                legend.tital = "BCL2",
                legend.labs = c("Negative", "Positive"),
                palette = c("blue", "red"))

g$table <- g$table +
        theme(plot.title = element_text(hjust = 0))

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/blc2_efs_conventional.jpg")
print(g)
dev.off()

# Distribution of MYC expression on tumor cells
summary(df.cd20$nor.myc)
# Min 0.00000 / 1Q 0.04862 / Median 0.06973 / 3Q 0.10430 / Max 1.00000

# Histogram of MYC expression on tumor cells
options(scipen = 999)
g <- qplot(nor.myc, data = df.cd20, geom = "histogram", color = "black",
           bins = 100)

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/myc.jpg")
print(g)
dev.off()

# Using intensity cutpoint (0.1)
df.cd20$myc.1 <- ifelse(df.cd20$nor.myc < 0.1, 0, 1)
myc.1 <- aggregate(myc.1 ~ SN, data = df.cd20, FUN = sum)
myc.pos <- inner_join(myc.1, tot, by = "SN")
myc.pos$myc.ratio <- (myc.pos$myc.1 / myc.pos$total) * 100
myc.pos <- dplyr::select(myc.pos, c(1,4))

# merge cell data and medical data
df.amc <- inner_join(myc.pos, med.data, by = "SN")

# cox proportional hazard model for myc ratio
fit <- coxph(Surv(mon.os, OS == 1) ~ myc.ratio, data = df.amc)
summary(fit) # p value = 0.0204 / exp(coef) = 1.011228

# 40% proportional cutpoint
df.amc$myc40 <- ifelse(df.amc$myc.ratio < 40, 0, 1)
df.amc$myc40 <- factor(df.amc$myc40, levels = c(0,1), 
                       labels = c("Negative", "Positive"))

# K-M curve (OS)
fit <- survfit(Surv(mon.os, OS == 1) ~ myc40, data = df.amc)
g <- ggsurvplot(fit,
                surv.plot.height = 0.8,
                pval = TRUE, pval.size = 4,
                xlab = "Time (month)", font.x = 12,
                ylab = "OS", font.y = 12,
                risk.table = TRUE, risk.table.height = 0.35,
                risk.table.fontsize = 3.5,
                risk.table.y.text.col = FALSE,
                legend = "right",
                legend.tital = "MYC",
                legend.labs = c("Negative", "Positive"),
                palette = c("blue", "red"))

g$table <- g$table +
        theme(plot.title = element_text(hjust = 0))

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/myc_os_conventional.jpg")
print(g)
dev.off()

# K-M curve (2yr EFS)
fit <- survfit(Surv(mon.2yr.efs, yr2.Event == 1) ~ myc40, 
               data = df.amc)
g <- ggsurvplot(fit,
                surv.plot.height = 0.8,
                pval = TRUE, pval.size = 4,
                xlab = "Time (month)", font.x = 12,
                ylab = "2yr EFS", font.y = 12,
                risk.table = TRUE, risk.table.height = 0.35,
                risk.table.fontsize = 3.5,
                risk.table.y.text.col = FALSE,
                legend = "right",
                legend.tital = "MYC",
                legend.labs = c("Negative", "Positive"),
                palette = c("blue", "red"))

g$table <- g$table +
        theme(plot.title = element_text(hjust = 0))

print(g)

jpeg(filename = "D:/Study/BCL2/bcl2 figure/myc_efs_conventional.jpg")
print(g)
dev.off()
###################################################################
############# BCL2 analysis using AQUA score ######################
###################################################################

