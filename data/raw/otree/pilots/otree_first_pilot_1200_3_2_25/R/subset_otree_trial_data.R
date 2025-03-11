# Load required libraries
library(tidyverse)

# Read the CSV file
data <- read_csv("data/raw/otree/otree_first_pilot_1200_3_2_25/all_apps_wide_2025-02-03.csv")

# First, create a vector of the IDs
participant_ids <- c("66b8830e5812d8df82589ef3",
                     "6274f9115634950cf5970f1c",
                     "5f305c9cf7fabe4adaf0de5d",
                     "5dd9a7ff98981d926df9ea18",
                     "63654bbd9fba7764c698ee95",
                     "645d03acd1302079adbc39f9",
                     "66b239b6a6912115ce1d9683",
                     "676edbd71b2d6c518eae7b45",
                     "5f24a57c030da70af5ffd2bf",
                     "5bcedec641f0710001c74dd5")

# Filter for specific participants and remove unwanted columns
cleaned_data <- data %>%
  # Filter for specific participant labels
  filter(participant.label %in% participant_ids) %>%
  # Remove columns starting with 'session'
  select(-starts_with("session")) %>%
  # Remove columns that are all NA
  select_if(~!all(is.na(.))) %>%
  # Remove columns that are all zeros
  select_if(~!all(. == 0))

base_patterns <- c(
  "main_task\\.[0-9]+\\.player\\.player1_choice_one",
  "main_task\\.[0-9]+\\.player\\.player2_choice_one",
  "main_task\\.[0-9]+\\.player\\.player3_choice_one",
  "main_task\\.[0-9]+\\.player\\.player4_choice_one",
  "main_task\\.[0-9]+\\.player\\.player1_choice_two",
  "main_task\\.[0-9]+\\.player\\.player2_choice_two",
  "main_task\\.[0-9]+\\.player\\.player3_choice_two",
  "main_task\\.[0-9]+\\.player\\.player4_choice_two",
  "main_task\\.[0-9]+\\.player\\.player2_choice1_accuracy",
  "main_task\\.[0-9]+\\.player\\.player3_choice1_accuracy",
  "main_task\\.[0-9]+\\.player\\.player2_choice2_accuracy",
  "main_task\\.[0-9]+\\.player\\.player3_choice2_accuracy",
  "main_task\\.[0-9]+\\.player\\.player1_loss_or_gain",
  "main_task\\.[0-9]+\\.player\\.player2_loss_or_gain",
  "main_task\\.[0-9]+\\.player\\.player3_loss_or_gain",
  "main_task\\.[0-9]+\\.player\\.player4_loss_or_gain",
  "main_task\\.[0-9]+\\.player\\.player2_choice1_computer",
  "main_task\\.[0-9]+\\.player\\.player3_choice1_computer",
  "main_task\\.[0-9]+\\.player\\.player2_choice2_computer",
  "main_task\\.[0-9]+\\.player\\.player3_choice2_computer",
  "main_task\\.[0-9]+\\.player\\.player4_choice2_computer",
  "main_task\\.[0-9]+\\.player\\.manual_second_choice",
  "main_task\\.[0-9]+\\.player\\.disconnection_streak",
  "main_task\\.[0-9]+\\.player\\.last_connection_time",
  "main_task\\.[0-9]+\\.player\\.last_check_time",
  "main_task\\.[0-9]+\\.player\\.consecutive_missed_checks",
  "main_task\\.[0-9]+\\.group\\.id_in_subsession",
  "main_task\\.[0-9]+\\.group\\.current_round",
  "main_task\\.[0-9]+\\.group\\.trial_number",
  "main_task\\.[0-9]+\\.group\\.round_reward_A",
  "main_task\\.[0-9]+\\.group\\.intertrial_interval",
  "main_task\\.[0-9]+\\.group\\.second_bet_timer_ended_executed",
  "main_task\\.[0-9]+\\.group\\.next_round_transition_time",
  "main_task\\.[0-9]+\\.group\\.reward_probability_A",
  "main_task\\.[0-9]+\\.group\\.reward_probability_B",
  "main_task\\.[0-9]+\\.group\\.seventy_percent_image",
  "main_task\\.[0-9]+\\.group\\.thirty_percent_image",
  "main_task\\.[0-9]+\\.group\\.bet_container_displayed",
  "main_task\\.[0-9]+\\.group\\.remaining_images_displayed",
  "main_task\\.[0-9]+\\.group\\.round_reward_set",
  "main_task\\.[0-9]+\\.group\\.disconnection_streaks"
)

# Remove these columns and save
cleaned_data <- cleaned_data %>%
  select(-matches(paste(base_patterns, collapse="|")))

# Add the new patterns
additional_patterns <- c(
  "main_task\\.[0-9]+\\.player\\.left_image",
  "main_task\\.[0-9]+\\.player\\.right_image",
  "main_task\\.[0-9]+\\.player\\.trial_reward",
  "main_task\\.[0-9]+\\.player\\.chosen_image_one",
  "main_task\\.[0-9]+\\.player\\.chosen_image_one_binary",
  "main_task\\.[0-9]+\\.player\\.chosen_image_two",
  "main_task\\.[0-9]+\\.player\\.chosen_image_two_binary",
  "main_task\\.[0-9]+\\.player\\.choice1_with",
  "main_task\\.[0-9]+\\.player\\.choice1_against",
  "main_task\\.[0-9]+\\.player\\.choice2_with",
  "main_task\\.[0-9]+\\.player\\.choice2_against",
  "main_task\\.[0-9]+\\.player\\.choice1_accuracy",
  "main_task\\.[0-9]+\\.player\\.choice2_accuracy",
  "main_task\\.[0-9]+\\.player\\.choice1_earnings",
  "main_task\\.[0-9]+\\.player\\.choice2_earnings",
  "main_task\\.[0-9]+\\.player\\.choice1_sum_earnings",
  "main_task\\.[0-9]+\\.player\\.choice2_sum_earnings",
  "main_task\\.[0-9]+\\.player\\.bonus_payment_score",
  "main_task\\.[0-9]+\\.player\\.base_payoff",
  "main_task\\.[0-9]+\\.player\\.loss_or_gain",
  "main_task\\.[0-9]+\\.player\\.player1_choice1_accuracy",
  "main_task\\.[0-9]+\\.player\\.player4_choice1_accuracy",
  "main_task\\.[0-9]+\\.player\\.player1_choice1_computer",
  "main_task\\.[0-9]+\\.player\\.player4_choice1_computer",
  "main_task\\.[0-9]+\\.group\\.my_page_load_time",
  "main_task\\.[0-9]+\\.group\\.round_reward_B"
)

# Combine all patterns
all_patterns <- c(base_patterns, additional_patterns)

# Remove these columns and save
cleaned_data <- cleaned_data %>%
  select(-matches(paste(all_patterns, collapse="|")))

# Save the cleaned data
write_csv(cleaned_data, "cleaned_data.csv")








extract_trial_data <- function(data) {
  patterns <- c(
    "player\\.id_in_group",
    "player\\.choice1",
    "player\\.choice2",
    "player\\.computer_choice1",
    "player\\.computer_choice2",
    "player\\.bet1",
    "player\\.bet2",
    "player\\.chosen_image_computer",
    "player\\.chosen_image_computer_two",
    "player\\.switch_vs_stay",
    "player\\.my_page_load_time",
    "player\\.individual_page_load_time",
    "player\\.initial_choice_time",
    "player\\.initial_bet_time",
    "player\\.second_choice_time",
    "player\\.second_bet_time",
    "player\\.computer_choice_one",
    "player\\.computer_bet_one",
    "player\\.computer_choice_two",
    "player\\.computer_bet_two",
    "player\\.player1_choice2_accuracy",
    "subsession\\.round_number"
  )
  
  long_data <- data %>%
    select(participant.label, 
           matches("^main_task\\.[0-9]+\\.")) %>%
    select(participant.label,
           matches(paste0("main_task\\.[0-9]+\\.", paste(patterns, collapse="|")))) %>%
    # Convert all columns to character type
    mutate(across(everything(), as.character)) %>%
    pivot_longer(
      cols = -participant.label,
      names_to = c("task", "trial", "component", "measure"),
      names_pattern = "(main_task)\\.([0-9]+)\\.(player|subsession)\\.(.*)",
      values_to = "value"
    ) %>%
    mutate(trial = as.numeric(trial)) %>%
    pivot_wider(
      id_cols = c(participant.label, trial),
      names_from = measure,
      values_from = value
    ) %>%
    arrange(participant.label, trial)
  
  return(long_data)
}

# Apply the transformation and save
trial_data <- extract_trial_data(cleaned_data)

write_csv(trial_data, "trial_by_trial_data.csv")
