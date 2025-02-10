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
