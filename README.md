# CFH_rotation
---
title: "Candice_pvalue"
output: html_notebook
---


```{r}

## Thickness bar graph
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)

data <- readxl::read_excel('/Users/candice/Desktop/Iowa_24Fall/MullinsLab/Choroidal_thickness.xlsx')

data <- na.omit(data)
data
# calculate Total, OS, OD average
summary_data <- data %>%
  group_by(Genotype) %>%
  summarize(
    Total = mean(Average_µm, na.rm = TRUE),
    OS = mean(Average_µm[Eye == "os"], na.rm = TRUE),
    OD = mean(Average_µm[Eye == "od"], na.rm = TRUE),
    .groups = "drop"
  )

# sample size
sample_sizes <- data %>%
  group_by(Genotype, Eye) %>%
  summarize(SampleSize = n(), .groups = "drop") %>%
  spread(Eye, SampleSize, fill = 0)

summary_data <- summary_data %>%
  mutate(
    Total_n = rowSums(sample_sizes[,-1]),  
    OS_n = sample_sizes$os,
    OD_n = sample_sizes$od
  )
summary_data
# long data format
summary_data_long <- melt(
  summary_data,
  id.vars = c("Genotype"),
  measure.vars = c("Total", "OS", "OD"),
  variable.name = "Group",
  value.name = "Average_µm"
)


sample_sizes_long <- melt(
  summary_data,
  id.vars = c("Genotype"),
  measure.vars = c("Total_n", "OS_n", "OD_n"),
  variable.name = "Group",
  value.name = "SampleSize"
)
sample_sizes_long$Group <- gsub("_n", "", sample_sizes_long$Group)  # get rid of _n
sample_sizes_long

plot_data <- merge(summary_data_long, sample_sizes_long, by = c("Genotype", "Group"))
plot_data

total_data <- subset(plot_data, Group == "Total")

# bar graph 
pic <- ggplot(plot_data, aes(x = Genotype, y = Average_µm, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # border
  geom_text(aes(label = paste0("n=", SampleSize)), 
            position = position_dodge(width = 0.8), vjust = -0.5) +  
  labs(
    title = "Average (µm) by Genotype and Group",
    x = "Genotype",
    y = "Average (µm)"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 250)) +
  scale_fill_brewer(palette = "Set2", name = "Group") +
  theme(legend.position = "top") 

ggsave(
  filename = "/Users/candice/Desktop/Iowa_24Fall/MullinsLab/Average_thickness_by_Genotype.pdf",  
  plot = pic,  
  device = "pdf",  
  width = 8,  
  height = 6,  
  units = "in"  
)

```

```{r}
## P-value calculation
# Create total data by averaging OS and OD values for each DonorID, Genotype, and MeasureType
total_data <- reshaped_data %>%
  filter(Eye %in% c("os", "od")) %>%  # Select OS and OD data
  group_by(DonorID, Genotype, MeasureType) %>%  # Group by DonorID, Genotype, and MeasureType
  summarise(
    Thickness = mean(Thickness, na.rm = TRUE),  # Average Thickness for OS and OD
    .groups = "drop"
  ) %>%
  mutate(Eye = "total")  # Add "total" label for Eye

# Add total data back to the original dataset
reshaped_data_with_total <- bind_rows(reshaped_data, total_data)

# Filter data for Total, OS, and OD
total_data <- reshaped_data_with_total %>% filter(Eye == "total")
os_data <- reshaped_data_with_total %>% filter(Eye == "os")
od_data <- reshaped_data_with_total %>% filter(Eye == "od")

# Define ANOVA analysis function
perform_anova <- function(data) {
  anova_result <- aov(Thickness ~ Genotype, data = data)  # Perform ANOVA
  summary_result <- summary(anova_result)                # ANOVA summary
  tukey_result <- TukeyHSD(anova_result)                 # Tukey HSD post-hoc test
  return(list("ANOVA" = summary_result, "TukeyHSD" = tukey_result))
}

# Perform ANOVA for Total, OS, and OD
total_results <- perform_anova(total_data)
os_results <- perform_anova(os_data)
od_results <- perform_anova(od_data)

# Perform unpaired t-test for OS vs OD (independent samples)
# Perform unpaired t-test for OS vs OD (independent samples)
os_vs_od_results <- reshaped_data %>%
  filter(Eye %in% c("os", "od")) %>%  # Select OS and OD data
  group_by(Genotype, MeasureType) %>%  # Group by Genotype and MeasureType
  summarise(
    p_value = t.test(
      Thickness[Eye == "os"],  # Subset OS data
      Thickness[Eye == "od"],  # Subset OD data
      paired = FALSE,          # Independent samples t-test
      na.rm = TRUE
    )$p.value,
    .groups = "drop"
  ) %>%
filter(MeasureType == "Average_µm") 


# Print results
print("Total Results (ANOVA)")
print(total_results$ANOVA)
print("Tukey HSD for Total")
print(total_results$TukeyHSD)

print("OS Results (ANOVA)")
print(os_results$ANOVA)
print("Tukey HSD for OS")
print(os_results$TukeyHSD)

print("OD Results (ANOVA)")
print(od_results$ANOVA)
print("Tukey HSD for OD")
print(od_results$TukeyHSD)

print("OS vs OD Results (Unpaired t-test)")
print(os_vs_od_results)

```

```{r}

## p-value table generation
library(dplyr)
library(tidyr)
library(gt)
library(readxl)

# Function to extract p-values from Tukey HSD results
extract_tukey_results <- function(tukey_result) {
  results <- as.data.frame(tukey_result$Genotype)  # Extract Tukey results for Genotype
  results <- tibble::rownames_to_column(results, var = "Comparison")  # Add row names as a column
  results <- results %>% 
    select(Comparison, `p adj`) %>%  # Select only Comparison and adjusted p-value
    rename(p_value = `p adj`)        # Rename column for clarity
  return(results)
}

# Extract Tukey HSD results for Total, OS, and OD
total_tukey_results <- extract_tukey_results(total_results$TukeyHSD)
os_tukey_results <- extract_tukey_results(os_results$TukeyHSD)
od_tukey_results <- extract_tukey_results(od_results$TukeyHSD)

# Add OS vs OD t-test results
os_vs_od_table <- os_vs_od_results %>%
  mutate(Comparison = paste(Genotype, "(OS vs OD)")) %>%  # Create comparison labels
  select(Comparison, p_value)

# Combine all results into a single table
final_table <- bind_rows(
  total_tukey_results %>% mutate(Category = "Total"),   # Add Category for Total
  os_tukey_results %>% mutate(Category = "OS"),        # Add Category for OS
  od_tukey_results %>% mutate(Category = "OD"),        # Add Category for OD
  os_vs_od_table %>% mutate(Category = "OS vs OD")     # Add Category for OS vs OD
)

# Add bold borders for specified rows
bold_border_rows <- c(3, 6, 9)

# Create the final table using gt
library(gt)
gt_table <- final_table %>%
  gt() %>%
  tab_header(title = "choroid thickness p-value in genotypes") %>%
  fmt_number(
    columns = vars(p_value),  # Format p-value to 4 decimal places
    decimals = 4
  ) %>%
  cols_label(
    Category = "Category",
    Comparison = "Comparison",
    p_value = "P-value"
  ) %>%
  tab_style(
    style = list(
      cell_borders(
        sides = "bottom",  # Add bold borders at the bottom of specified rows
        color = "black",
        weight = px(1)
      )
    ),
    locations = cells_body(rows = bold_border_rows)
  )

# Print the table
print(gt_table)

# Save the table as an HTML file (optional)
gtsave(gt_table, "/Users/candice/Desktop/Iowa_24Fall/MullinsLab/P-value_Results.html")



```

```{r}
##Age bar graph
# Summarize Age data by Genotype (calculate Total, OS, and OD averages and sample sizes)
age_summary_data <- data %>%
  na.omit() %>%
  group_by(Genotype) %>%
  summarize(
    Total_Age = mean(Age, na.rm = TRUE),              # Average Age for Total
    Total_n = sum(!is.na(Age)),                       # Sample size for Total
    OS_Age = mean(Age[Eye == "os"], na.rm = TRUE),    # Average Age for OS
    OS_n = sum(Eye == "os" & !is.na(Age)),            # Sample size for OS
    OD_Age = mean(Age[Eye == "od"], na.rm = TRUE),    # Average Age for OD
    OD_n = sum(Eye == "od" & !is.na(Age)),            # Sample size for OD
    .groups = "drop"
  )

# Convert Age averages to long format
age_summary_data_long <- melt(
  age_summary_data,
  id.vars = c("Genotype"),
  measure.vars = c("Total_Age", "OS_Age", "OD_Age"),
  variable.name = "Group",
  value.name = "Average_Age"
)

# Convert sample sizes to long format
age_sample_sizes_long <- melt(
  age_summary_data,
  id.vars = c("Genotype"),
  measure.vars = c("Total_n", "OS_n", "OD_n"),
  variable.name = "Group",
  value.name = "SampleSize"
)

# Remove suffixes to unify Group names
age_summary_data_long$Group <- gsub("_Age", "", age_summary_data_long$Group)
age_sample_sizes_long$Group <- gsub("_n", "", age_sample_sizes_long$Group)

# Merge Age averages and sample sizes into a single dataset
age_plot_data <- merge(age_summary_data_long, age_sample_sizes_long, by = c("Genotype", "Group"))

# Create bar plot for Average Age by Genotype and Group
age_plot <- ggplot(age_plot_data, aes(x = Genotype, y = Average_Age, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Add borders to bars
  geom_text(aes(label = paste0("n=", SampleSize)),                                        # Add sample size labels
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(
    title = "Average Age by Genotype and Group",
    x = "Genotype",
    y = "Average Age"
  ) +
  theme_minimal() +  # Use minimal theme
  scale_y_continuous(limits = c(0, max(age_plot_data$Average_Age, na.rm = TRUE) + 5)) +
  scale_fill_brewer(palette = "Set2", name = "Group") +  # Use color palette for groups
  theme(legend.position = "top")  # Place legend at the top

# Save the plot as a PDF
ggsave(
  filename = "/Users/candice/Desktop/Iowa_24Fall/MullinsLab/Average_Age_by_Genotype.pdf",
  plot = age_plot,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in"
)



```
```{r}

## Age Bar Graph
# Filter data for Age
age_data <- reshaped_data %>% filter(MeasureType == "Age")

# Create total data for Age
age_total_data <- age_data %>%
  filter(Eye %in% c("os", "od")) %>%
  group_by(DonorID, Genotype) %>%
  summarise(
    Thickness = mean(Thickness, na.rm = TRUE),  # Average Age for OS and OD
    .groups = "drop"
  ) %>%
  mutate(Eye = "total")  # Add "total" label for Eye

# Add total data back to Age dataset
age_data_with_total <- bind_rows(age_data, age_total_data)

# Filter Total, OS, and OD data for Age
age_total_data <- age_data_with_total %>% filter(Eye == "total")
age_os_data <- age_data_with_total %>% filter(Eye == "os")
age_od_data <- age_data_with_total %>% filter(Eye == "od")

# Perform ANOVA for Age (Total, OS, OD)
age_total_results <- perform_anova(age_total_data)
age_os_results <- perform_anova(age_os_data)
age_od_results <- perform_anova(age_od_data)

# Perform unpaired t-test for OS vs OD for Age
age_os_vs_od_results <- age_data %>%
  filter(Eye %in% c("os", "od")) %>%
  group_by(Genotype) %>%
  summarise(
    p_value = t.test(
      Thickness[Eye == "os"],  # Age for OS
      Thickness[Eye == "od"],  # Age for OD
      paired = FALSE,          # Independent samples t-test
      na.rm = TRUE
    )$p.value,
    .groups = "drop"
  )

# Print results for Age
print("Total Results (ANOVA) for Age")
print(age_total_results$ANOVA)
print("Tukey HSD for Total (Age)")
print(age_total_results$TukeyHSD)

print("OS Results (ANOVA) for Age")
print(age_os_results$ANOVA)
print("Tukey HSD for OS (Age)")
print(age_os_results$TukeyHSD)

print("OD Results (ANOVA) for Age")
print(age_od_results$ANOVA)
print("Tukey HSD for OD (Age)")
print(age_od_results$TukeyHSD)

print("OS vs OD Results (Unpaired t-test) for Age")
print(age_os_vs_od_results)



```
```{r}
##Age p-value generation
# Function to extract p-values from Tukey HSD results for Age
extract_tukey_results <- function(tukey_result) {
  results <- as.data.frame(tukey_result$Genotype)  # Extract Tukey results for Genotype
  results <- tibble::rownames_to_column(results, var = "Comparison")  # Add row names as a column
  results <- results %>% 
    select(Comparison, `p adj`) %>%  # Select only Comparison and adjusted p-value
    rename(p_value = `p adj`)        # Rename column for clarity
  return(results)
}

# Extract Tukey HSD results for Total, OS, and OD for Age
age_total_tukey_results <- extract_tukey_results(age_total_results$TukeyHSD)
age_os_tukey_results <- extract_tukey_results(age_os_results$TukeyHSD)
age_od_tukey_results <- extract_tukey_results(age_od_results$TukeyHSD)

# Add OS vs OD t-test results for Age
age_os_vs_od_table <- age_os_vs_od_results %>%
  mutate(Comparison = paste(Genotype, "(OS vs OD)")) %>%  # Create comparison labels
  select(Comparison, p_value)

# Combine all Age results into a single table
age_final_table <- bind_rows(
  age_total_tukey_results %>% mutate(Category = "Total"),   # Add Category for Total
  age_os_tukey_results %>% mutate(Category = "OS"),        # Add Category for OS
  age_od_tukey_results %>% mutate(Category = "OD"),        # Add Category for OD
  age_os_vs_od_table %>% mutate(Category = "OS vs OD")     # Add Category for OS vs OD
)

# Add bold borders for specified rows
bold_border_rows <- c(3, 6, 9)

# Create the final table for Age using gt
library(gt)
age_gt_table <- age_final_table %>%
  gt() %>%
  tab_header(title = "Age p-value in genotypes") %>%
  fmt_number(
    columns = vars(p_value),  # Format p-value to 4 decimal places
    decimals = 4
  ) %>%
  cols_label(
    Category = "Category",
    Comparison = "Comparison",
    p_value = "P-value"
  ) %>%
  tab_style(
    style = list(
      cell_borders(
        sides = "bottom",  # Add bold borders at the bottom of specified rows
        color = "black",
        weight = px(1)
      )
    ),
    locations = cells_body(rows = bold_border_rows)
  )

# Print the table for Age
print(age_gt_table)

# Save the table as an HTML file (optional)
gtsave(age_gt_table, "/Users/candice/Desktop/Iowa_24Fall/MullinsLab/Age_P-value_Results.html")


```


