library(tidyplots)

cyto <- Valentino_plates_MSD_01182024

cyto <- na.omit(cyto)
cytoD187 <- cyto[cyto$Day == 'Day 187',]

cytoD187$Concentration <- as.numeric(cytoD187$Concentration)
cytoD187$Concentration_log10 <- log10(cytoD187$Concentration)
cytoD187$Vaccine1 <- as.factor(cytoD187$Vaccine1)


cytoD187$Cytokine <- gsub("TNF-Œ±", "TNF-α", cytoD187$Cytokine)
cytoD187$Cytokine <- gsub("IFN-Œ±2a", "IFN-α2a", cytoD187$Cytokine)
cytoD187$Cytokine <- gsub("IFN-Œ≤", "IFN-β", cytoD187$Cytokine)
cytoD187$Cytokine <- gsub("IFN-Œ≥", "IFN-γ", cytoD187$Cytokine)

vaccine_colors <- c("Cervarix" = "#E4AD2E", "Gardasil" = "#007A73")

cytoD187 %>% 
  tidyplot(x = Vaccine1, y = Concentration, color = Vaccine1) %>%
  adjust_colors(new_colors = vaccine_colors) %>%
  add_data_points_beeswarm(white_border = T) %>%
  #add_data_labels_repel(label = twin) %>%
  add_boxplot(alpha = 0.5) %>%
  #add_mean_dash() %>%
  #add_sem_errorbar() %>%
  #add_test_asterisks(method = "fisher_test", p.adjust.method = "bonferroni", label = "p.signif", hide_info = T) %>%
  add_test_pvalue(method = "wilcox_test", p.adjust.method = "BH",  hide_info = T, hide.ns = F) %>%
  adjust_x_axis_title("") %>%
  adjust_x_axis(labels = c("", ""), rotate_labels = 45) %>%
  adjust_y_axis_title("Concentration", family = "helvetica", fontsize = 8, color = "black") %>%
  adjust_y_axis(limits = c(0, NA), padding = c(0.1,0.15)) %>%
  adjust_legend_title("Vaccine", family = "helvetica", fontsize = 10, color = "black") %>%
  adjust_font(family = "Helvetica", fontsize = 8, color = "black", face = "bold") %>%
  #adjust_size(width = 100, height = 100, unit = "mm") %>%
  split_plot(by = Cytokine)