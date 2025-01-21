library(readr)
merged_test_data <- read_csv("data/preprocessed/merged_test_data_10.csv")

data <- merged_test_data

# Get key parameters
n_trials <- length(unique(data$group.trial_number))
n_total_participants <- 300

# Create expanded dataframe structure
new_data <- data.frame()

for(trial in 1:n_trials) {
  # Get template data for this trial
  trial_template <- data[data$group.trial_number == trial, ][1:10, ]
  
  # Create data for all 300 participants for this trial
  trial_data <- trial_template[rep(1:10, length.out = n_total_participants), ]
  
  # Update participant IDs and labels
  trial_data$participant.id_in_session <- 1:n_total_participants
  trial_data$participant.label <- paste0("P", 1:n_total_participants)
  
  # Generate new unique codes
  trial_data$participant.code <- replicate(n_total_participants, 
                                           paste0(sample(LETTERS, 8, replace=TRUE), collapse=""))
  
  # Set choice_switch_across_trials to NA for trial 1
  if(trial == 1) {
    trial_data$choice_switch_across_trials <- NA
  }
  
  # Append to main dataframe
  new_data <- rbind(new_data, trial_data)
}

# Verify structure
head(new_data[, c("participant.id_in_session", "group.trial_number")], 20)

# Save the expanded dataset
write.csv(new_data, "merged_test_data.csv", row.names = FALSE)