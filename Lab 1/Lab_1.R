library(bnlearn)
library(Rgraphviz)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")
BiocManager::install("RBGL")
library(gRain)

data("asia")

data <- asia



##### TASK 1 #####
# Show that multiple runs of the hill-climbing algorithm can
# return non-equivalent Bayesian
# network (BN) structures. Explain why this happens

# Runs hill climbing algorithm with 5 restarts.
set.seed(12345)
bn_hc <- hc(x = data, restart = 5)
bn_hc <- cpdag(bn_hc)
plot(bn_hc)

set.seed(12345)
bn_hc_2 <- hc(x = data,
              score = "aic")
bn_hc_2 <- cpdag(bn_hc_2)
bn_hc <- cpdag(bn_hc_2)

plot(bn_hc_2)

# Check equality
print(all.equal(bn_hc, bn_hc_2))



##### TASK 2 #####
# Learn a BN from 80% of the data set. Use BN to classify the rest of the data.
# Use exact or approximate inference with the help of the bnlearn and gRain
# packages. You are not allowed to use functions such as predict. Report the
# confusion matrix, i.e. true/false positives/negatives. Compare your results
# with those of the true Asia BN.


set.seed(12345)
n <- dim(data)[1]
id <- sample(1:n, floor(n * 0.8))
train <- data[id,]
test <- data[-id,]

bn_structure <- hc(train, restart = 5)
plot(bn_structure)

bn_fit <- bn.fit(bn_structure, train, method = "bayes")

# True Asia BN
bn_true <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
bn_true_fit <- bn.fit(bn_structure, train, method = "bayes")
plot(bn_true)


#Converts the fit from above into a gRain-object
bn_grain <- as.grain(bn_fit)
bn_true_grain <- as.grain(bn_true_fit)

#Compiles the junction
j_tree <- compile(bn_grain)
j_tree_true <- compile(bn_true_grain)



prediction <- function(tree, data, obs, pred) {
  predictions <- c()
  
  for (i in 1:dim(data)[1]) {
    x_vals <- c()
    
    for (j in obs) {
      x_vals[j] <- if (data[i, j] == "yes")
        "yes"
      else
        "no"
    }
    
    
    evidence <- setEvidence(tree, obs, x_vals)
    probability <- querygrain(evidence, pred)$S["yes"]
    predictions[i] <- if (probability >= 0.5)
      "yes"
    else
      "no"
    
  }
  return(predictions)
}

bn_pred <-
  prediction(j_tree, test, c("A", "T", "L", "B", "E", "X", "D"), c("S"))
bn_pred_true <-
  prediction(j_tree_true, test, c("A", "T", "L", "B", "E", "X", "D"), c("S"))


# Confusion matrices by comparing with original data
confusion_matrix_fit <- table(bn_pred, test$S)
print(confusion_matrix_fit)
confusion_matrix_true <- table(bn_pred_true, test$S)
print(confusion_matrix_true)

# As one can see the confusion matrix is the exact same as the ones above.
# The reason for this is that the node is independent from

##### TASK 3 #####
# In the previous exercise, you classified the variable S given observations for all the
# rest of the variables. Now, you are asked to classify S given observations only for the
# so-called Markov blanket of S, i.e. its parents plus its children plus the parents of its
# children minus S itself. Report again the confusion matrix

bn_pred_mb <-
  prediction(j_tree, test, mb(bn_fit, "S"), c("S"))
bn_pred_mb_true <-
  prediction(j_tree_true, test, mb(bn_true_fit, "S"), c("S"))

confusion_matrix_fit_mb <- table(bn_pred_mb, test$S)
print(confusion_matrix_fit_mb)
confusion_matrix_true_mb <- table(bn_pred_mb_true, test$S)
print(confusion_matrix_true_mb)


##### TASK 4 #####
# Repeat the exercise (2) using a naive Bayes classifier, i.e. the predictive variables are
# independent given the class variable. See p. 380 in Bishop's book or Wikipedia for
# more information on the naive Bayes classifier. Model the naive Bayes classifier as a
# BN. You have to create the BN by hand, i.e. you are not allowed to use the function
# naive.bayes from the bnlearn package.

bn_naive <- empty.graph(c("A", "S", "T", "L", "B", "E", "X", "D"))
arc.set <- matrix(
  c("S", "A",
    "S", "T",
    "S", "L",
    "S", "B",
    "S", "E",
    "S", "X",
    "S", "D"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to"))
)
# Above is the BN set. In a naive bayes network all nodes are dependent on the independent classifier
# node. If S is observed, all other nodes are independent

arcs(bn_naive) <- arc.set

bn_naive_fit <- bn.fit(bn_naive, data=train)

bn_naive_grain <- as.grain(bn_naive_fit)

bn_naive_comp <- compile(bn_naive_grain)

bn_naive_pred <- prediction(bn_naive_comp, test, c("A", "T", "L", "B", "E", "X", "D"), c("S"))

confusion_matrix_fit_naive <- table(bn_naive_pred, test$S)
confusion_matrix_fit_naive





































