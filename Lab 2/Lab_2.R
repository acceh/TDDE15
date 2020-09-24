library(HMM)

##### TASK 1 #####

# 10 Different places
states <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
symbol_states <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Equal probability of it moving and staying
transition_probs <- t(matrix(
                   c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5,
                     0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5), 
                     nrow=10, ncol=10))


# Equal probability of it being i-2,i-1,i,i+1 and i+2
emission_probs <-  matrix(
                    c(0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2,
                      0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2,
                      0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0,
                      0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0,
                      0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0,
                      0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0,
                      0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0,
                      0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2,
                      0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2,
                      0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2), 
                      nrow=10, ncol=10)


hmm <-
  initHMM(States = states,
          Symbols = symbol_states,
          startProbs = rep(0.1, 10),
          transProbs = transition_probs,
          emissionProbs = emission_probs
          )



##### TASK 2 #####

set.seed(12345)
simulated_steps <- simHMM(hmm, 100)

##### TASK 3 #####

filter_function <- function(alphas) {
  filtered <- matrix(NA,
                     nrow = dim(alphas)[1],
                     ncol = dim(alphas)[2])
  
  for (t in 1:dim(alphas)[2]) {
    filtered[,t] <- alphas[, t] / sum(alphas[, t])
  }
  
  return(filtered)
}


smooth_function <- function(alphas, betas) {
  smoothed <- matrix(NA,
                     nrow = dim(alphas)[1],
                     ncol = dim(alphas)[2])
  
  for (t in 1:dim(alphas)[2]) {
    smoothed[, t] <-
      (alphas[, t] * betas[, t]) / (sum(alphas[, t] * betas[, t]))
  }
  
  return (smoothed)
}


hmm_forward_alpha <-
  exp(forward(hmm = hmm,
              observation = simulated_steps$observation))

hmm_backward_beta <-
  exp(backward(hmm = hmm,
               observation = simulated_steps$observation))

filtered <- filter_function(hmm_forward_alpha)

smoothed <-
  smooth_function(hmm_forward_alpha, hmm_backward_beta)

viterbi <- viterbi(hmm, simulated_steps$observation)


##### TASK 4 #####

accuracy_function <- function(prediction, true) {
  confusion_matrix <- table(prediction, true)
  accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
  return (accuracy)
}

filtered_prediction <- apply(t(filtered), MARGIN = 1, which.max)
smoothed_prediction <- apply(t(smoothed), MARGIN = 1, which.max)


accuracy_function(filtered_prediction, simulated_steps$states)
accuracy_function(smoothed_prediction, simulated_steps$states)
accuracy_function(viterbi, simulated_steps$states)


##### TASK 5 #####
#Different seed
set.seed(67890)
simulated_steps_seeded <- simHMM(hmm, 100)

hmm_forward_alpha_seeded <-
  exp(forward(hmm = hmm,
              observation = simulated_steps_seeded$observation))

hmm_backward_beta_seeded <-
  exp(backward(hmm = hmm,
               observation = simulated_steps_seeded$observation))

filtered_seeded <- filter_function(hmm_forward_alpha_seeded)

smoothed_seeded <-
  smooth_function(hmm_forward_alpha_seeded, hmm_backward_beta_seeded)

viterbi_seeded <- viterbi(hmm, simulated_steps$observation)

filtered_prediction_seeded <-
  apply(t(filtered_seeded), MARGIN = 1, which.max)
smoothed_prediction_seeded <-
  apply(t(smoothed_seeded), MARGIN = 1, which.max)


accuracy_function(filtered_prediction_seeded, simulated_steps_seeded$states)
accuracy_function(smoothed_prediction_seeded, simulated_steps_seeded$states)
accuracy_function(viterbi_seeded, simulated_steps_seeded$states)


##### TASK 6 #####
library(entropy)


simulated_steps_300 <- simHMM(hmm, 300)

hmm_forward_alpha_300 <-
  exp(forward(hmm = hmm,
              observation = simulated_steps_300$observation))

hmm_backward_beta_300 <-
  exp(backward(hmm = hmm,
               observation = simulated_steps_300$observation))

filtered_300 <- filter_function(hmm_forward_alpha_300)

smoothed_300 <-
  smooth_function(hmm_forward_alpha_300, hmm_backward_beta_300)

viterbi_300 <- viterbi(hmm, simulated_steps_300$observation)

smoothed_prediction_300 <-
  apply(t(smoothed_300), MARGIN = 1, which.max)
filtered_prediction_300 <-
  apply(t(filtered_300), MARGIN = 1, which.max)

entropy_filtered_100 <- entropy.empirical(filtered_prediction)
entropy_filtered_300 <- entropy.empirical(filtered_prediction_300)
entropy_smoothed_100 <- entropy.empirical(smoothed_prediction)
entropy_smoothed_300 <- entropy.empirical(smoothed_prediction_300)





print("Entropy for 100 simulated samples for filtered distribution:")
entropy_filtered_100
print("Entropy for 300 simulated samples for filtered distribution:")
entropy_filtered_300
print("Entropy for 100 simulated samples for smoothed distribution:")
entropy_smoothed_100
print("Entropy for 300 simulated samples for smoothed distribution:")
entropy_smoothed_300




print("Accuracy for 100 simulated samples for filtered distribution:")
accuracy_function(filtered_prediction, simulated_steps$states)
print("Accuracy for 300 simulated samples for filtered distribution:")
accuracy_function(filtered_prediction_300, simulated_steps_300$states)
print("Accuracy for 100 simulated samples for smoothed distribution:")
accuracy_function(smoothed_prediction, simulated_steps$states)
print("Accuracy for 300 simulated samples for smoothed distribution:")
accuracy_function(smoothed_prediction_300, simulated_steps_300$states)



##### TASK 7 #####

step_101 <- transition_probs %*% t(filtered)[100,]

print(step_101)

