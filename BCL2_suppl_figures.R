###############################################################################
########Represent various intensity of BCL2 expression in tumor cells##########
###############################################################################

# label cells with quantile group
df.cd20$Q1 <- ifelse(df.cd20$nor.bcl2 < 0.07177, 1, 0)
df.cd20$Q2 <- ifelse(df.cd20$nor.bcl2 >= 0.07177 & df.cd20$nor.bcl2 < 0.09441, 1, 0)
df.cd20$Q3 <- ifelse(df.cd20$nor.bcl2 >= 0.09441 & df.cd20$nor.bcl2 < 0.13310, 1, 0)
df.cd20$Q4 <- ifelse(df.cd20$nor.bcl2 >= 0.13310, 1, 0)

# calculate the number of cells in each quantile group
Q1 <- aggregate(Q1 ~ SN, data = df.cd20, FUN = sum)
Q2 <- aggregate(Q2 ~ SN, data = df.cd20, FUN = sum)
Q3 <- aggregate(Q3 ~ SN, data = df.cd20, FUN = sum)
Q4 <- aggregate(Q4 ~ SN, data = df.cd20, FUN = sum)

pos <- inner_join(Q1, Q2, by = "SN")
pos <- inner_join(pos, Q3, by = "SN")
pos <- inner_join(pos, Q4, by = "SN")
pos <- inner_join(pos, tot, by = "SN")

pos$Q1.ratio <- pos$Q1 / pos$total * 100
pos$Q2.ratio <- pos$Q2 / pos$total * 100
pos$Q3.ratio <- pos$Q3 / pos$total * 100
pos$Q4.ratio <- pos$Q4 / pos$total * 100

pos <- dplyr::select(pos, c(1, 7, 8, 9, 10))

# stacked barplot according to the cell numbers of each intensity quantile
df.bar <- reshape2::melt(pos, id = "SN")

g <- ggplot(df.bar, aes(fill = variable, y = value, x = SN))
g + geom_bar(stat = "identity") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position = "bottom") +
        labs(x = "Cases", y = "Proportion (%)") +
        scale_fill_manual(values = alpha(c("light blue", "dark green", "yellow", "dark red"), 0.8),
                          name = "Intensity", 
                          labels = c("Negative", "Weak", "Moderate", "Strong"))

jpeg(filename = "D:/Study/BCL2/bcl2 figure/intensity.jpg", width = 24, height = 8, units = "cm", res = 300)
g <- ggplot(df.bar, aes(fill = variable, y = value, x = SN))
g + geom_bar(stat = "identity") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position = "bottom") +
        labs(x = "Cases", y = "Proportion (%)") +
        scale_fill_manual(values = alpha(c("light blue", "dark green", "yellow", "dark red"), 0.8),
                          name = "Intensity", 
                          labels = c("Negative", "Weak", "Moderate", "Strong"))
dev.off()

###############################################################################
#########Represent various intensity of MYC expression in tumor cells##########
###############################################################################
# label cells with quantile group
df.cd20$Q1 <- ifelse(df.cd20$nor.myc < 0.04862, 1, 0)
df.cd20$Q2 <- ifelse(df.cd20$nor.myc >= 0.04862 & df.cd20$nor.myc < 0.06973, 1, 0)
df.cd20$Q3 <- ifelse(df.cd20$nor.myc >= 0.06973 & df.cd20$nor.myc < 0.10430, 1, 0)
df.cd20$Q4 <- ifelse(df.cd20$nor.myc >= 0.10430, 1, 0)

# calculate the number of cells in each quantile group
Q1 <- aggregate(Q1 ~ SN, data = df.cd20, FUN = sum)
Q2 <- aggregate(Q2 ~ SN, data = df.cd20, FUN = sum)
Q3 <- aggregate(Q3 ~ SN, data = df.cd20, FUN = sum)
Q4 <- aggregate(Q4 ~ SN, data = df.cd20, FUN = sum)

pos <- inner_join(Q1, Q2, by = "SN")
pos <- inner_join(pos, Q3, by = "SN")
pos <- inner_join(pos, Q4, by = "SN")
pos <- inner_join(pos, tot, by = "SN")

pos$Q1.ratio <- pos$Q1 / pos$total * 100
pos$Q2.ratio <- pos$Q2 / pos$total * 100
pos$Q3.ratio <- pos$Q3 / pos$total * 100
pos$Q4.ratio <- pos$Q4 / pos$total * 100

pos <- dplyr::select(pos, c(1, 7, 8, 9, 10))

# stacked barplot according to the cell numbers of each intensity quantile
df.bar <- reshape2::melt(pos, id = "SN")

g <- ggplot(df.bar, aes(fill = variable, y = value, x = SN))
g + geom_bar(stat = "identity") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position = "bottom") +
        labs(x = "Cases", y = "Proportion (%)") +
        scale_fill_manual(values = alpha(c("light blue", "dark green", "yellow", "dark red"), 0.8),
                          name = "MYC intensity", 
                          labels = c("Negative", "Weak", "Moderate", "Strong"))

jpeg(filename = "D:/Study/BCL2/bcl2 figure/myc_intensity.jpg", 
     width = 24, height = 8, units = "cm", res = 300)
g <- ggplot(df.bar, aes(fill = variable, y = value, x = SN))
g + geom_bar(stat = "identity") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), legend.position = "bottom") +
        labs(x = "Cases", y = "Proportion (%)") +
        scale_fill_manual(values = alpha(c("light blue", "dark green", "yellow", "dark red"), 0.8),
                          name = "MYC intensity", 
                          labels = c("Negative", "Weak", "Moderate", "Strong"))
dev.off()
