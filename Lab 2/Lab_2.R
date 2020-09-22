library(HMM)



##### TASK 1 #####

# 10 Different places
states <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
symbol_states <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Equal probability of it moving and staying
transition_probs <- matrix(
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
                     nrow=10, ncol=10)


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
simulated_steps <- simHMM(hmm,100)

simulated_steps

##### TASK 3 #####




filter_function <- function(alphas) {
  filtered <- matrix(NA,
                     nrow = dim(alphas)[1],
                     ncol = dim(alphas)[2])
  for (t in 1:dim(alphas)[1]) {
    # ??(zt)/???zt ??(zt)
    filtered[t,] <- alphas[t,] / sum(alphas[t,])
  }
  
  return(filtered)
}


smooth_function <- function(alphas, betas) {
  smoothed <- matrix(NA,
                     nrow = dim(alphas)[1],
                     ncol = dim(alphas)[2])
  
  for (t in 1:dim(alphas)[1]) {
    # ??(zt)*??(zt) /???(zt ??(zt)?? (zt))
    smoothed[t, ] <-
      (alphas[t, ] * betas[t, ]) / (sum(alphas[t, ] * betas[t, ]))
  }
  
  return (smoothed)
}


hmm_forward_alpha <-
  exp(forward(
    hmm = hmm,
    observation = simulated_steps$observation
  ))

hmm_forward_beta <-
 exp( backward(
    hmm = hmm,
    observation = simulated_steps$observation
  ))

filtered <- filter_function(hmm_forward_alpha)

smoothed <-
  smooth_function(hmm_forward_alpha, hmm_forward_beta)

viterbi <- viterbi(hmm, simulated_steps$observation)


##### TASK 4 #####

accuracy_function <- function(prediction, true) {
  confusion_matrix <- table(prediction, true)
  accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
  return (accuracy)
}

filtered_prediction <- apply(t(filtered), MARGIN = 1, which.max)
smoothed_prediction <- apply(smoothed, MARGIN = 1, which.max)


accuracy_function(filtered_prediction)
