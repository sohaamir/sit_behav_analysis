# Load required libraries
library(tidyverse)

# Read the CSV file
data <- read_csv("data/raw/otree/otree_second_pilot_1200_5_2_25_d6lvy1ol/submission_2025-02-05.csv")

library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape data to long format and remove NA values
data_long <- data %>%
  select(starts_with("player.")) %>%
  select(player.task_understanding, player.task_difficulty, 
         player.engagement, player.influence, 
         player.real_players, player.attention_focus) %>%
  drop_na() %>%
  pivot_longer(cols = everything(),
               names_to = "measure",
               values_to = "value") %>%
  mutate(measure = gsub("player.", "", measure))

# Calculate means for labels
means_data <- data_long %>%
  group_by(measure) %>%
  summarize(mean_value = mean(value))

# Create the violin plot
ggplot(data_long, aes(x = measure, y = value)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.1), 
             size = 1, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.3, color = "black") +
  geom_text(data = means_data, 
            aes(x = measure, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 32, face = "bold"),
    axis.title.y = element_text(size = 32, face = "bold")
  ) +
  labs(x = "Measure", y = "Rating")









# Load required libraries
library(tidyverse)

# Read the CSV file
data <- read_csv("data/raw/otree/otree_second_pilot_1200_5_2_25_d6lvy1ol/main_task_2025-02-05.csv")

library(ggplot2)
library(tidyr)
library(dplyr)

# Create vector of IDs to remove
ids_to_remove <- c(9, 10, 13, 14, 15, 16, 17, 18, 19, 20)

# Subset data excluding those IDs 
subset_data <- data[!data$participant.id_in_session %in% ids_to_remove,]

# Save as CSV
write.csv(subset_data, "rewards_data.csv", row.names = FALSE)

subset_data <- subset_data[, c("participant.id_in_session", "group.trial_number",
                               "player.choice1_accuracy", "player.choice2_accuracy")]

write.csv(subset_data, "rewards.csv", row.names = FALSE)

# Calculate totals for each participant
totals <- subset_data %>%
  group_by(participant.id_in_session) %>%
  summarize(
    choice1_total = sum(player.choice1_accuracy),
    choice2_total = sum(player.choice2_accuracy)
  )

# Convert to long format for plotting
totals_long <- totals %>%
  pivot_longer(
    cols = c(choice1_total, choice2_total),
    names_to = "choice_type",
    values_to = "total"
  ) %>%
  mutate(choice_type = factor(choice_type, 
                              labels = c("Choice 1", "Choice 2")))

# Calculate means for labels
means_data <- totals_long %>%
  group_by(choice_type) %>%
  summarize(mean_value = mean(total))

# Create the plot
ggplot(totals_long, aes(x = choice_type, y = total)) +
  geom_violin(fill = "#3182BD", alpha = 0.6) +
  geom_point(size = 2, alpha = 0.7, position = position_jitter(width = 0.1)) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.3, color = "black") +
  geom_text(data = means_data, 
            aes(x = choice_type, y = mean_value, label = round(mean_value, 2)),
            vjust = -0.5, size = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 20)),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  ) +
  labs(x = "Choice Type", 
       y = "Trials with 70% option chosen (out of 60") +
  scale_y_continuous(breaks = seq(0, max(totals_long$total), by = 5))




