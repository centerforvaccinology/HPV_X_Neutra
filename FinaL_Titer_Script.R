library(tidyverse)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggrepel)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library("gmodels")
library("DescTools")
library("qqplotr")
library("dplyr")
library(tidyplots)

titers <- read_excel("HPV_Total_Cohort_Titers.xlsx", 
                     sheet = "Sheet2",
                     col_types = c("numeric", "numeric", "numeric", 
                                   "date", "date", "date", "text", "numeric", 
                                   "numeric", "text", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", "numeric"), 
                     na = "na")
View(titers)            

titers$Vax <- as.factor(titers$Vax)
levels(titers$Vax) <- c("Cervarix", "Gardasil")
class(titers$Vax)
titers$Visit <- as.factor(titers$Visit)
levels(titers$Visit) <- c("Day_0", "Day_187")
titers$Assay <- as.factor(titers$Assay)
titers$HPV6 <- as.numeric(titers$HPV6)
titers$HPV16 <- as.numeric(titers$HPV16)
titers$HPV18 <- as.numeric(titers$HPV18)
titers$HPV31 <- as.numeric(titers$HPV31)
titers$HPV33 <- as.numeric(titers$HPV33)
titers$HPV45 <- as.numeric(titers$HPV45)
titers$HPV52 <- as.numeric(titers$HPV52)
titers$HPV58 <- as.numeric(titers$HPV58)

pbna <- titers[(titers$Assay == "PBNA"),]
pbna$HPV11 <- NULL

#GMT calculation and Plot PBNA
x <- pbna |> 
  dplyr::filter(Visit == "Day_187") |> 
  dplyr::filter(Vax == "Cervarix")

GMT_HPV6_Day187_C <- Gmean(x$HPV6)
SD_HPV6_D187_C <- Gsd(x$HPV6)
GMT_HPV16_Day187_C <- Gmean(x$HPV16)
SD_HPV16_D187_C <- Gsd(x$HPV16)
GMT_HPV18_Day187_C <- Gmean(x$HPV18)
SD_HPV18_D187_C <- Gsd(x$HPV18)
GMT_HPV31_Day187_C <- Gmean(x$HPV31)
SD_HPV31_D187_C <- Gsd(x$HPV31)
GMT_HPV33_Day187_C <- Gmean(x$HPV33)
SD_HPV33_D187_C <- Gsd(x$HPV33)
GMT_HPV45_Day187_C <- Gmean(x$HPV45)
SD_HPV45_D187_C <- Gsd(x$HPV45)
GMT_HPV52_Day187_C <- Gmean(x$HPV52)
SD_HPV52_D187_C <- Gsd(x$HPV52)
GMT_HPV58_Day187_C <- Gmean(x$HPV58)
SD_HPV58_D187_C <- Gsd(x$HPV58)

x1 <- pbna |> 
  filter(Visit == "Day_187") |> 
  filter(Vax == "Gardasil")

GMT_HPV6_Day187_G <- Gmean(x1$HPV6)
SD_HPV6_D187_G <- Gsd(x1$HPV6)
GMT_HPV16_Day187_G <- Gmean(x1$HPV16)
SD_HPV16_D187_G <- Gsd(x1$HPV16)
GMT_HPV18_Day187_G <- Gmean(x1$HPV18)
SD_HPV18_D187_G <- Gsd(x1$HPV18)
GMT_HPV31_Day187_G <- Gmean(x1$HPV31)
SD_HPV31_D187_G <- Gsd(x1$HPV31)
GMT_HPV33_Day187_G <- Gmean(x1$HPV33)
SD_HPV33_D187_G <- Gsd(x1$HPV33)
GMT_HPV45_Day187_G <- Gmean(x1$HPV45)
SD_HPV45_D187_G <- Gsd(x1$HPV45)
GMT_HPV52_Day187_G <- Gmean(x1$HPV52)
SD_HPV52_D187_G <- Gsd(x1$HPV52)
GMT_HPV58_Day187_G <- Gmean(x1$HPV58)
SD_HPV58_D187_G <- Gsd(x1$HPV58)

y <- pbna |> 
  filter(Visit == "Day_0") |> 
  filter(Vax == "Cervarix")

GMT_HPV6_Day0_C <- Gmean(y$HPV6)
SD_HPV6_D0_C <- Gsd(y$HPV6)
GMT_HPV16_Day0_C <- Gmean(y$HPV16)
SD_HPV16_D0_C <- Gsd(y$HPV16)
GMT_HPV18_Day0_C <- Gmean(y$HPV18)
SD_HPV18_D0_C <- Gsd(y$HPV18)
GMT_HPV31_Day0_C <- Gmean(y$HPV31)
SD_HPV31_D0_C <- Gsd(y$HPV31)
GMT_HPV33_Day0_C <- Gmean(y$HPV33)
SD_HPV33_D0_C <- Gsd(y$HPV33)
GMT_HPV45_Day0_C <- Gmean(y$HPV45)
SD_HPV45_D0_C <- Gsd(y$HPV45)
GMT_HPV52_Day0_C <- Gmean(y$HPV52)
SD_HPV52_D0_C <- Gsd(y$HPV52)
GMT_HPV58_Day0_C <- Gmean(y$HPV58)
SD_HPV58_D0_C <- Gsd(y$HPV58)

y1 <- pbna |> 
  filter(Visit == "Day_0") |> 
  filter(Vax == "Gardasil")

GMT_HPV6_Day0_G <- Gmean(y1$HPV6)
SD_HPV6_D0_G <- Gsd(y1$HPV6)
GMT_HPV16_Day0_G <- Gmean(y1$HPV16)
SD_HPV16_D0_G <- Gsd(y1$HPV16)
GMT_HPV18_Day0_G <- Gmean(y1$HPV18)
SD_HPV18_D0_G <- Gsd(y1$HPV18)
GMT_HPV31_Day0_G <- Gmean(y1$HPV31)
SD_HPV31_D0_G <- Gsd(y1$HPV31)
GMT_HPV33_Day0_G <- Gmean(y1$HPV33)
SD_HPV33_D0_G <- Gsd(y1$HPV33)
GMT_HPV45_Day0_G <- Gmean(y1$HPV45)
SD_HPV45_D0_G <- Gsd(y1$HPV45)
GMT_HPV52_Day0_G <- Gmean(y1$HPV52)
SD_HPV52_D0_G <- Gsd(y1$HPV52)
GMT_HPV58_Day0_G <- Gmean(y1$HPV58)
SD_HPV58_D0_G <- Gsd(y1$HPV58)

GMTtable <- data.frame("HPV_Type" = c("HPV6", "HPV16", "HPV18","HPV31", "HPV33", "HPV45", "HPV52", "HPV58",
                                      "HPV6", "HPV16", "HPV18","HPV31", "HPV33", "HPV45", "HPV52", "HPV58",
                                      "HPV6", "HPV16", "HPV18","HPV31", "HPV33", "HPV45", "HPV52", "HPV58",
                                      "HPV6", "HPV16", "HPV18","HPV31", "HPV33", "HPV45", "HPV52", "HPV58"),
                       "Day" = c("Day 0","Day 0","Day 0","Day 0",
                                 "Day 0","Day 0","Day 0","Day 0",
                                 "Day 0","Day 0","Day 0","Day 0",
                                 "Day 0","Day 0","Day 0","Day 0",
                                 "Day 187","Day 187","Day 187","Day 187",
                                 "Day 187","Day 187","Day 187","Day 187",
                                 "Day 187","Day 187","Day 187","Day 187",
                                 "Day 187","Day 187","Day 187","Day 187"),
                       "Vaccine" = c("Cervarix","Cervarix","Cervarix","Cervarix",
                                     "Cervarix","Cervarix","Cervarix","Cervarix",
                                     "Gardasil", "Gardasil","Gardasil","Gardasil",
                                     "Gardasil","Gardasil","Gardasil","Gardasil",
                                     "Cervarix","Cervarix","Cervarix","Cervarix",
                                     "Cervarix","Cervarix","Cervarix","Cervarix",
                                     "Gardasil","Gardasil","Gardasil","Gardasil",
                                     "Gardasil","Gardasil","Gardasil","Gardasil"),
                       "GMT" = c(GMT_HPV6_Day0_C, GMT_HPV16_Day0_C, GMT_HPV18_Day0_C, GMT_HPV31_Day0_C, 
                                 GMT_HPV33_Day0_C, GMT_HPV45_Day0_C, GMT_HPV52_Day0_C, GMT_HPV58_Day0_C,
                                 GMT_HPV6_Day0_G, GMT_HPV16_Day0_G, GMT_HPV18_Day0_G, GMT_HPV31_Day0_G, 
                                 GMT_HPV33_Day0_G, GMT_HPV45_Day0_G, GMT_HPV52_Day0_G, GMT_HPV58_Day0_G,
                                 GMT_HPV6_Day187_C, GMT_HPV16_Day187_C, GMT_HPV18_Day187_C, GMT_HPV31_Day187_C, 
                                 GMT_HPV33_Day187_C, GMT_HPV45_Day187_C, GMT_HPV52_Day187_C, GMT_HPV58_Day187_C,
                                 GMT_HPV6_Day187_G, GMT_HPV16_Day187_G, GMT_HPV18_Day187_G, GMT_HPV31_Day187_G, 
                                 GMT_HPV33_Day187_G, GMT_HPV45_Day187_G, GMT_HPV52_Day187_G, GMT_HPV58_Day187_G),
                       "sd" = c(SD_HPV6_D0_C, SD_HPV16_D0_C, SD_HPV18_D0_C, SD_HPV31_D0_C, 
                                SD_HPV33_D0_C, SD_HPV45_D0_C, SD_HPV52_D0_C, SD_HPV58_D0_C,
                                SD_HPV6_D0_G, SD_HPV16_D0_G, SD_HPV18_D0_G, SD_HPV31_D0_G, 
                                SD_HPV33_D0_G, SD_HPV45_D0_G, SD_HPV52_D0_G, SD_HPV58_D0_G,
                                SD_HPV6_D187_C, SD_HPV16_D187_C, SD_HPV18_D187_C, SD_HPV31_D187_C, 
                                SD_HPV33_D187_C, SD_HPV45_D187_C, SD_HPV52_D187_C, SD_HPV58_D187_C,
                                SD_HPV6_D187_G, SD_HPV16_D187_G, SD_HPV18_D187_G, SD_HPV31_D187_G, 
                                SD_HPV33_D187_G, SD_HPV45_D187_G, SD_HPV52_D187_G, SD_HPV58_D187_G))

GMTtable$Day <- as.factor(GMTtable$Day)
GMTtable$HPV_Type <- as.factor(GMTtable$HPV_Type)
GMTtable$Vaccine <- as.factor(GMTtable$Vaccine)

Pz_combined <- GMTtable |> 
  ggplot(aes(Day, GMT, color = HPV_Type, group = HPV_Type)) + 
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = GMT - sd, ymax = GMT + sd), width = 0.2,
                position = position_dodge(0.05)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "", color = "", y = "log(GMT)", x = "") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  facet_wrap(~Vaccine, scales = "free_y") + # Separate facets for each vaccine
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    axis.text.x = element_text(size = 12, colour = "black"),  # X-axis text color
    axis.text.y = element_text(size = 10, colour = "black"),  # Y-axis text color
    axis.title.x = element_text(size = 12, colour = "black"),  # X-axis title color
    axis.title.y = element_text(size = 12, colour = "black"),  # Y-axis title color
    legend.text = element_text(size = 12, colour = "black"),  # Legend text color
    strip.text = element_text(size = 12, face = "bold", colour = "black") # Facet label color
  )

Pz_combined
ggsave(filename = "GMT Neutralization.pdf", plot = Pz_combined, device = "pdf", width = 25, height = 15, units = "cm", dpi = "retina")

#Produce descriptive statistics by group
#n – The number of observations for each treatment.
#mean – The mean value for each treatment.
#sd – The standard deviation of each treatment.
#stderr – The standard error of each treatment.  That is the standard deviation / sqrt (n).
#LCL, UCL – The upper and lower confidence intervals of the mean.  That is to say, you can be 95% certain that the true mean falls between the lower and upper values specified for each treatment group assuming a normal distribution. 
#median – The median value for each treatment.
#min, max – The minimum and maximum value for each treatment.
#IQR – The inner quartile range of each treatment. That is the 75th percentile –  25th percentile.
#LCLmed, UCLmed – The 95% confidence interval for the median.

Summary_HPV6 <- pbna |> select(Visit, Vax, HPV6) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV6, na.rm = TRUE), 
            sd = sd(HPV6, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV6, na.rm = TRUE),
            min = min(HPV6, na.rm = TRUE), 
            max = max(HPV6, na.rm = TRUE),
            IQR = IQR(HPV6, na.rm = TRUE),
            LCLmed = MedianCI(HPV6, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV6, na.rm=TRUE)[3])

Summary_HPV16 <- pbna |> select(Visit, Vax, HPV16) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV16, na.rm = TRUE), 
            sd = sd(HPV16, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV16, na.rm = TRUE),
            min = min(HPV16, na.rm = TRUE), 
            max = max(HPV16, na.rm = TRUE),
            IQR = IQR(HPV16, na.rm = TRUE),
            LCLmed = MedianCI(HPV16, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV16, na.rm=TRUE)[3])

Summary_HPV18 <- pbna |> select(Visit, Vax, HPV18) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV18, na.rm = TRUE), 
            sd = sd(HPV18, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV18, na.rm = TRUE),
            min = min(HPV18, na.rm = TRUE), 
            max = max(HPV18, na.rm = TRUE),
            IQR = IQR(HPV18, na.rm = TRUE),
            LCLmed = MedianCI(HPV18, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV18, na.rm=TRUE)[3])

Summary_HPV31 <- pbna |> select(Visit, Vax, HPV31) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV31, na.rm = TRUE), 
            sd = sd(HPV31, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV31, na.rm = TRUE),
            min = min(HPV31, na.rm = TRUE), 
            max = max(HPV31, na.rm = TRUE),
            IQR = IQR(HPV31, na.rm = TRUE),
            LCLmed = MedianCI(HPV31, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV31, na.rm=TRUE)[3])

Summary_HPV33 <- pbna |> select(Visit, Vax, HPV33) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV33, na.rm = TRUE), 
            sd = sd(HPV33, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV33, na.rm = TRUE),
            min = min(HPV33, na.rm = TRUE), 
            max = max(HPV33, na.rm = TRUE),
            IQR = IQR(HPV33, na.rm = TRUE),
            LCLmed = MedianCI(HPV33, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV33, na.rm=TRUE)[3])

Summary_HPV45 <- pbna |> select(Visit, Vax, HPV45) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV45, na.rm = TRUE), 
            sd = sd(HPV45, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV45, na.rm = TRUE),
            min = min(HPV45, na.rm = TRUE), 
            max = max(HPV45, na.rm = TRUE),
            IQR = IQR(HPV45, na.rm = TRUE),
            LCLmed = MedianCI(HPV45, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV45, na.rm=TRUE)[3])

Summary_HPV52 <- pbna |> select(Visit, Vax, HPV52) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV52, na.rm = TRUE), 
            sd = sd(HPV52, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV52, na.rm = TRUE),
            min = min(HPV52, na.rm = TRUE), 
            max = max(HPV52, na.rm = TRUE),
            IQR = IQR(HPV52, na.rm = TRUE),
            LCLmed = MedianCI(HPV52, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV52, na.rm=TRUE)[3])

Summary_HPV58 <- pbna |> select(Visit, Vax, HPV58) |>  group_by(Vax, Visit) |> 
  summarise(n = n(), 
            mean = mean(HPV58, na.rm = TRUE), 
            sd = sd(HPV58, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(HPV58, na.rm = TRUE),
            min = min(HPV58, na.rm = TRUE), 
            max = max(HPV58, na.rm = TRUE),
            IQR = IQR(HPV58, na.rm = TRUE),
            LCLmed = MedianCI(HPV58, na.rm=TRUE)[2],
            UCLmed = MedianCI(HPV58, na.rm=TRUE)[3])

Summary_HPV6$Type <- "HPV6"
Summary_HPV16$Type <- "HPV16"
Summary_HPV18$Type <- "HPV18"
Summary_HPV31$Type <- "HPV31"
Summary_HPV33$Type <- "HPV33"
Summary_HPV45$Type <- "HPV45"
Summary_HPV52$Type <- "HPV52"
Summary_HPV58$Type <- "HPV58"

Summary <- rbind(Summary_HPV16, Summary_HPV18, Summary_HPV31, Summary_HPV33,
                 Summary_HPV45, Summary_HPV52, Summary_HPV58, Summary_HPV6)

write.csv(Summary, "Summary Statistics Neutralizing Antibodies.csv", row.names=FALSE)

#Perform the Mann-Whitney U test
#W – This value represents the Wilcoxon test statistic.  
#The Wilcoxon test statistic is the sum of the ranks in sample 1 minus n1*(n1+1)/2. 
#n1 is the number of observations in sample 1.
#p-value – The p-value corresponding to the two-sided test based on the standard normal (Z) distribution. 
#95% confidence interval – The 95% confidence interval on the difference between the number of bugs that survived under the effects of spray C vs spray D.
#difference in location – This value corresponds to the Hodges-Lehmann Estimate of the location parameter differences between sprays C and D.

x <- pbna |>
  filter(Visit == "Day_187")

HPV6_MWU <- wilcox.test(HPV6 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV16_MWU <- wilcox.test(HPV16 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV18_MWU <- wilcox.test(HPV18 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV31_MWU <- wilcox.test(HPV31 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV33_MWU <- wilcox.test(HPV33 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV45_MWU <- wilcox.test(HPV45 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV52_MWU <- wilcox.test(HPV52 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)
HPV58_MWU <- wilcox.test(HPV58 ~ Vax, data=x, na.rm=TRUE, exact=FALSE, conf.int=TRUE)

MWU <- rbind(HPV6_MWU, HPV16_MWU, HPV18_MWU, HPV31_MWU, 
             HPV33_MWU, HPV45_MWU, HPV52_MWU, HPV58_MWU)

write.csv(MWU, "MannWhitneyU Test Neutralizing Antibodies Day 187.csv", row.names=TRUE)

#All HPV types pbna boxplots in one figure per vaccine per visit
x$logHPV6 <- log10(x$HPV6)
x$logHPV16 <- log10(x$HPV16)
x$logHPV18 <- log10(x$HPV18)
x$logHPV31 <- log10(x$HPV31)
x$logHPV33 <- log10(x$HPV33)
x$logHPV45 <- log10(x$HPV45)
x$logHPV52 <- log10(x$HPV52)
x$logHPV58 <- log10(x$HPV58)

colnames(x)[11:18] <- c("6", "16", "18", "31", "33","45", "52","58")
colnames(x)[19:26] <- c("HPV6", "HPV16", "HPV18", "HPV31", "HPV33","HPV45", "HPV52","HPV58")

y <- cbind(stack(x[, 19:26]), Vaccine = x$Vaccine)
y$day <- "Day 187"
y

ggplot(y, aes(x = day, y = values, fill = Vaccine, color = Vaccine)) +
  geom_boxplot(alpha = 0.5) +
  stat_boxplot(geom = "errorbar") + 
  facet_wrap(~ ind, scales = "free") +
  labs(
    y = "log10(anti-HPV)",
    x = "",
    title = "",
    fill = "Vaccine") +
  scale_fill_manual(values = c("Cervarix" = "#E4AD2E", "Gardasil" = "#007A73")) +
  theme_classic()

y %>% 
  tidyplot(x = day, y = values, color = Vaccine) %>%
  adjust_colors(new_colors = c("Cervarix" = "#E4AD2E", "Gardasil" = "#007A73")) %>%
  add_data_points_beeswarm(white_border = T) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "wilcoxon", p.adjust.method = "BH", label = "p.adj.signif", hide_info = T) %>%
  add_test_pvalue(method = "wilcoxon", p.adjust.method = "BH",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Log10(anti-HPV)", family = "helvetica", fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0,7)) %>%
  adjust_legend_title("Vaccine", family = "helvetica", fontsize = 8, color = "black", face = "bold") %>%
  adjust_font(family = "Helvetica", fontsize = 10, color = "black", face = "bold") %>%
  adjust_size(width = 100, height = 100, unit = "mm") %>%
  split_plot(by = ind,width = 40, height = 40, unit = "mm")

ggsave(filename = "GMT Neutralization per type per vaccine.pdf", width = 25, height = 20, units = "cm", dpi = "retina", device = "pdf")

#Euclidian Distance
pbna <- titers[(titers$Assay == "PBNA"),]
pbna$HPV11 <- NULL
pbna$HPV6 <- NULL
pbna <- pbna[(pbna$Visit == "Day_187"),]
pbna[,11:17] <- log(pbna[,11:17])

euclidean <- function(a, b) sqrt(sum((a - b)^2))

a <- pbna[pbna$SubjectID == '3', 11:17]
b <- pbna[pbna$SubjectID == '4', 11:17]
c <- pbna[pbna$SubjectID == '5', 11:17]
d <- pbna[pbna$SubjectID == '6', 11:17]
e <- pbna[pbna$SubjectID == '7', 11:17]
f <- pbna[pbna$SubjectID == '8', 11:17]
g <- pbna[pbna$SubjectID == '9', 11:17]
h <- pbna[pbna$SubjectID == '10', 11:17]
i <- pbna[pbna$SubjectID == '11', 11:17]
j <- pbna[pbna$SubjectID == '12', 11:17]
k <- pbna[pbna$SubjectID == '13', 11:17]
l <- pbna[pbna$SubjectID == '14', 11:17]

euclidiandistances <- data.frame("1" = euclidean(a, b),
                                 "2" = euclidean(c, d),
                                 "3" = euclidean(e, f),
                                 "4" = euclidean(g, h),
                                 "5" = euclidean(i,j),
                                 "6" = euclidean(k, l))

colnames(euclidiandistances) <- c("Twin pair 1", "Twin pair 2", "Twin pair 3", "Twin pair 4", "Twin pair 5", "Twin pair 6")
rownames(euclidiandistances) <- "EuclidianDistance"
write.csv(euclidiandistances, "Euclidian Distances between Twins based on HPV titres.csv")
euclidiandistances <- t(euclidiandistances)
euclidiandistances <- as.data.frame(euclidiandistances)
euclidiandistances$Twinpair <- c("1", "2", "3", "4", "5", "6")

ggplot(data=euclidiandistances, aes(x=Twinpair, y=EuclidianDistance)) +
  geom_bar(stat="identity",
           color = c("#ed9d82", "#89baf0","#a1cfa1","#f2af44","#edbeeb", "#f0ec7d"),
           fill = c("#ed9d82", "#89baf0","#a1cfa1","#f2af44","#edbeeb", "#f0ec7d")) +
  geom_text(aes(label=EuclidianDistance), vjust = -0.2) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 12),  # Base font settings
    axis.text.x = element_text(size = 12, colour = "black"),  # X-axis text color
    axis.text.y = element_text(size = 12, colour = "black"),  # Y-axis text color
    axis.title.x = element_text(size = 12, colour = "black"),  # X-axis title color
    axis.title.y = element_text(size = 12, colour = "black"),  # Y-axis title color
    legend.text = element_text(size = 12, colour = "black"),  # Legend text color
    strip.text = element_text(size = 12, face = "bold", colour = "black") # Facet label color
  )
