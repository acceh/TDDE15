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

simulated_steps$states <- c()

##### TASK 3 #####


#Function for alphas for the forward part of the backward forward algorithm
create_alpha <- function(poss_states, hmm, obs, alpha_prev) {
  alpha <- NULL
  
  for (state in poss_states) {
    emission <- hmm$emissionProbs[state, obs]
    
    
    weighted_transition <- sum(sapply(poss_states, function(z) {
      alpha_prev[z] * hmm$transProbs[z, state]
    }))
    
    alpha[state] <- emission * weighted_transition
  }
  return (alpha)
}

create_beta <- function(poss_states, hmm, obs_next, beta_next) {
  beta <- NULL
  
  
  for (state in poss_states) {
    beta[state] <- sum(sapply(poss_states, function(z) {
      beta_next[z] * hmm$emissionProbs[z, obs_next] * hmm$transProbs[state, z]
    }))
  }
  return (beta)
}



forward_backward_function <-
  function(poss_states,
           hmm,
           obs_vars) {
    #Gets the length of the list with observed vars from the simulated hmm steps.
    T <- length(obs_vars)
    # Initalize alpha and beta for the algorithm
    alpha <- matrix(NA, nrow = T, ncol = length(poss_states))
    beta <- matrix(NA, nrow = T, ncol = length(poss_states))
    
    
    first_obs <- obs_vars[1]
    initial <- hmm$startProbs
    alpha[1,] <- hmm$emissionProbs[, first_obs] * initial
    
    #Forward
    for (t in 2:T) {
      # Getting alpha for time t
      alpha[t,] <- create_alpha(
        poss_states = poss_states,
        hmm = hmm,
        obs = obs_vars[t],
        alpha_prev = alpha[t - 1,]
      )
    }
    
    
    beta[T,] <- 1
    
    #Backward
    for (t in (T - 1):1) {
      # Getting beta for time step t
      beta[t,] <- create_beta(
        poss_states = poss_states,
        hmm = hmm,
        obs = obs_vars[t + 1],
        beta_next = beta[t + 1,]
      )
    }
    return (list(alpha = alpha, beta = beta))
    
    
  }


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


hmm_forward_backward <-
  forward_backward_function(
    poss_states = 1:10,
    hmm = hmm,
    obs_vars = simulated_steps$observation
  )

filtered <- filter_function(hmm_forward_backward$alpha)

smoothed <-
  smooth_function(hmm_forward_backward$alpha, hmm_forward_backward$beta)




