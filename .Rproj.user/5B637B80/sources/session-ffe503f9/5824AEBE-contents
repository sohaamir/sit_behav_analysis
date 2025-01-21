############################################################################################
############################################################################################
############################################################################################
# subsetting the data - only keeping columns that we are interested in

# Load required libraries
library(tidyverse)

# Read the CSV file
data <- read_csv("data/raw/otree/test_data_two_groups.csv")

# Remove unneeded columns
data <- data[, !(names(data) %in% c(
  "participant._is_bot", "participant._index_in_pages", "player.last_check_time",
  "participant._max_page_index", "participant._current_app_name", 
  "participant._current_page_name", "participant.time_started_utc",
  "participant.visited", "participant.mturk_worker_id",
  "participant.mturk_assignment_id", "participant.payoff",
  "player.role", "player.payoff", "player.consecutive_missed_checks",
  "player.individual_page_load_time", "player.base_payoff",
  "player.bonus_payoff", "player.total_payoff",
  "player.computer_choice_one", "player.computer_bet_one",
  "player.computer_choice_two", "player.computer_bet_two",
  "player.player_1_computer_choice_one", "player.player_2_computer_choice_one",
  "player.player_3_computer_choice_one", "player.player_4_computer_choice_one",
  "player.player_1_computer_choice_two", "player.player_2_computer_choice_two",
  "player.player_3_computer_choice_two", "player.player_4_computer_choice_two",
  "player.manual_second_choice", "player.disconnection_streak",
  "player.is_bot", "player.last_connection_time",
  "group.current_round", "group.second_bet_timer_ended_executed",
  "group.next_round_transition_time", "group.reversal_rounds",
  "group.bet_container_displayed", "group.remaining_images_displayed",
  "group.round_reward_set", "group.all_players_loaded",
  "group.players_loaded_count", "group.disconnected_players",
  "group.bot_players", "group.active_bots",
  "group.disconnection_streaks", "session.code",
  "session.label", "session.mturk_HITId",
  "session.mturk_HITGroupId", "session.comment",
  "session.is_demo"
))]


# Convert rewards to -1 and 1 instead of 0 and 1
data$player.trial_reward[data$player.trial_reward == 0] <- -1
data$player.loss_or_gain[data$player.loss_or_gain == 0] <- -1
data$player.loss_or_gain_player1[data$player.loss_or_gain_player1 == 0] <- -1
data$player.loss_or_gain_player2[data$player.loss_or_gain_player2 == 0] <- -1
data$player.loss_or_gain_player3[data$player.loss_or_gain_player3 == 0] <- -1
data$player.loss_or_gain_player4[data$player.loss_or_gain_player4 == 0] <- -1



# Convert choice with/against to decimals instead of raw numbers

# Create a function to apply the division
divide_by_4 <- function(x) {
  ifelse(x == 0, 0, x/4)
}

# Apply the transformation to each column
data$player.choice1_with <- divide_by_4(data$player.choice1_with)
data$player.choice1_against <- divide_by_4(data$player.choice1_against)
data$player.choice2_with <- divide_by_4(data$player.choice2_with)
data$player.choice2_against <- divide_by_4(data$player.choice2_against)


# Save as CSV
write.csv(data, "data/preprocessed/otree/subsetted_test_data.csv", row.names = FALSE)


############################################################################################
############################################################################################
############################################################################################
# creating the dataList


# Create new list
testdataList <- list()

# Set nSubjects and nTrials
testdataList$nSubjects <- max(data$participant.id_in_session)
testdataList$nTrials <- as.integer(max(data$group.trial_number))

# Add 1D arrays for group assignments
testdataList$group <- array(NA, dim = testdataList$nSubjects)
testdataList$groupId <- array(NA, dim = testdataList$nSubjects)

# Create 2D arrays
testdataList$choice1 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice2 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$bet1 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$bet2 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice1Acc <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice2Acc <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$reward <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$winprob <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice1with <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice1against <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice2with <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$choice2against <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))
testdataList$chswitch <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials))

# Create all arrays
testdataList$otherChoice1 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials, 4))
testdataList$otherChoice2 <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials, 4))
testdataList$otherChoice1Acc <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials, 4))
testdataList$otherChoice2Acc <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials, 4))
testdataList$otherReward <- array(NA, dim = c(testdataList$nSubjects, testdataList$nTrials, 4))

# Fill the arrays
for(subject in 1:testdataList$nSubjects) {
  for(trial in 1:testdataList$nTrials) {
    subj_trial_idx <- which(data$participant.id_in_session == subject & 
                              data$group.trial_number == trial)
    
    if(length(subj_trial_idx) > 0) {
      # Only fill group assignments on first trial
      if(trial == 1) {
        testdataList$group[subject] <- data$group.id_in_subsession[subj_trial_idx]
        testdataList$groupId[subject] <- data$player.id_in_group[subj_trial_idx]
      }
      
      # 3D arrays
      testdataList$otherChoice1[subject, trial, ] <- c(
        data$player.player_1_choice_one[subj_trial_idx],
        data$player.player_2_choice_one[subj_trial_idx],
        data$player.player_3_choice_one[subj_trial_idx],
        data$player.player_4_choice_one[subj_trial_idx]
      )
      
      testdataList$otherChoice2[subject, trial, ] <- c(
        data$player.player_1_choice_two[subj_trial_idx],
        data$player.player_2_choice_two[subj_trial_idx],
        data$player.player_3_choice_two[subj_trial_idx],
        data$player.player_4_choice_two[subj_trial_idx]
      )
      
      testdataList$otherChoice1Acc[subject, trial, ] <- c(
        data$player.player1_choice1_accuracy[subj_trial_idx],
        data$player.player2_choice1_accuracy[subj_trial_idx],
        data$player.player3_choice1_accuracy[subj_trial_idx],
        data$player.player4_choice1_accuracy[subj_trial_idx]
      )
      
      testdataList$otherChoice2Acc[subject, trial, ] <- c(
        data$player.player1_choice2_accuracy[subj_trial_idx],
        data$player.player2_choice2_accuracy[subj_trial_idx],
        data$player.player3_choice2_accuracy[subj_trial_idx],
        data$player.player4_choice2_accuracy[subj_trial_idx]
      )
      
      testdataList$otherReward[subject, trial, ] <- c(
        data$player.loss_or_gain_player1[subj_trial_idx],
        data$player.loss_or_gain_player2[subj_trial_idx],
        data$player.loss_or_gain_player3[subj_trial_idx],
        data$player.loss_or_gain_player4[subj_trial_idx]
      )
      
      # 2D arrays
      testdataList$choice1[subject, trial] <- data$player.chosen_image_one_binary[subj_trial_idx]
      testdataList$choice2[subject, trial] <- data$player.chosen_image_two_binary[subj_trial_idx]
      testdataList$bet1[subject, trial] <- data$player.bet1[subj_trial_idx]
      testdataList$bet2[subject, trial] <- data$player.bet2[subj_trial_idx]
      testdataList$choice1Acc[subject, trial] <- data$player.choice1_accuracy[subj_trial_idx]
      testdataList$choice2Acc[subject, trial] <- data$player.choice2_accuracy[subj_trial_idx]
      testdataList$reward[subject, trial] <- data$player.trial_reward[subj_trial_idx]
      testdataList$winprob[subject, trial] <- data$group.reward_probability_A[subj_trial_idx]
      testdataList$choice1with[subject, trial] <- data$player.choice1_with[subj_trial_idx]
      testdataList$choice1against[subject, trial] <- data$player.choice1_against[subj_trial_idx]
      testdataList$choice2with[subject, trial] <- data$player.choice2_with[subj_trial_idx]
      testdataList$choice2against[subject, trial] <- data$player.choice2_against[subj_trial_idx]
      testdataList$chswitch[subject, trial] <- data$player.switch_vs_stay[subj_trial_idx]
    }
  }
}

# Save the testdataList object
save(testdataList, file = "data/rdata/testdata.RData")
