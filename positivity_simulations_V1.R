# This script is a mordified version of the R script provided in Kang et al.
# (2016). This conducts simulation studies to assess the effect of arbitrary
# cutoffs of propensity scores.

# Load required libraries.
suppressMessages(library(survey))
suppressMessages(library(gbm))
suppressMessages(library(tree))
suppressMessages(library(rpart))
suppressMessages(library(MatchIt))
suppressMessages(library(MASS))
suppressMessages(library(cem))
suppressMessages(library(randomForest))

# Set up dimension parameters.
N.OBS <- 1000
N.COVARIATES <- c(4)
N.NPOSITIVITY <- c(1)
N.SAMPLES <- 200
N.VALIDATIONS <- 10
N.STRATA <- 40
N.CUTS <- 4
N.SIMULATIONS <- length(N.COVARIATES)
outcome.formula <- "Y.observed ~ exposure"

# Define general functions.
expit <- function(a) {
  exp(a) / (1 + exp(a))
}

ComputeExposureFormula <- function(n.covariates) {
  exposure.formula <- as.formula(
    paste("exposure ~ 0 +",
          paste(paste("X", 1:n.covariates, sep = ""), collapse = "+"),
          sep = "")
  )
}

ComputeStabilizedWeights <- function(exposure, propensity.scores) {
  first.term <- exposure * mean(exposure) / propensity.scores
  second.term <- (1 - exposure) * mean(1 - exposure) / (1 - propensity.scores)
  weights <- first.term + second.term
  return(weights)
}

GenerateData <- function(n.obs, n.covariates, n.npositivity) {
  # Generate confounders.
  X <- matrix(NA, n.obs, n.covariates)
  for (j in 1:n.covariates) {
    X[, j] <- sample(x = seq(0.1, 1, 0.1), size = n.obs, replace = TRUE)
  }
  colnames(X) <- paste("X", 1:n.covariates, sep = "")

  # Define the "non-positivity" subjects by one categorical confounder, X1, and
  # the sum of two categorical confounders, X12 = X1 + X2.
  if (n.npositivity == 1) {
    idx.np.low <- (X[, 1] == 0.1)
    idx.np.upp <- (X[, 1] == 1.0)
    idx.p <- (X[, 1] > 0.1 & X[, 1] < 1.0)
  } else if (n.npositivity == 2) {
    x12 <- apply(X[, 1:2], 1, sum)
    idx.np.low <- (x12 <= quantile(x12, probs = 0.1))
    idx.np.upp <- (x12 <= quantile(x12, probs = 0.9))
    idx.p <- (x12 > quantile(x12, probs = 0.1) & x12 < quantile(x12, probs = 0.9))
  }

  #beta <- sample(x = c (1, -1), size = n.covariates, replace = TRUE)
  beta <- c(rep(1, (n.covariates / 2)), rep(-1, (n.covariates / 2)))
  eta <- X[idx.p, ] %*% beta

  # Generate the true propensity scores: P(T = 1 | x) = expit(X * beta).
  propensity.scores <- rep(NA, n.obs)
  propensity.scores[idx.p] <- exp(eta) / (1 + exp(eta))

  # Structurally zero, no treatment group will be found, and hence this group
  # needs to be removed.
  propensity.scores[idx.np.low] <- 0
  propensity.scores[idx.np.upp] <- 1

  # Generate the binary exposure variable.
  exposure <- rep(NA, n.obs)
  exposure <- 1 * (propensity.scores > runif(n.obs))

  # Propensity structure for the simulation study.
  Y1 <- rep(NA, n.obs)
  Y0 <- rep(NA, n.obs)

  # For e(x) = 0.
  Y1[idx.np.low] <-  rbinom(n = sum(idx.np.low), size = 1, prob = 0.1)
  Y0[idx.np.low] <-  rbinom(n = sum(idx.np.low), size = 1, prob = 0.9)

  # For e(x) = 1.
  Y1[idx.np.upp] <-  rbinom(n = sum(idx.np.upp), size = 1, prob = 0.1)
  Y0[idx.np.upp] <-  rbinom(n = sum(idx.np.upp), size = 1, prob = 0.9)

  # For 0 < e(x) < 1.
  Y1[idx.p] <-  rbinom(n = sum(idx.p), size = 1, prob = 0.6)
  Y0[idx.p] <-  rbinom(n = sum(idx.p), size = 1, prob = 0.4)

  Y.observed <- exposure * Y1 + (1 - exposure) * Y0
  true.weights <- ComputeStabilizedWeights(exposure, propensity.scores)
  true.weights[is.na(true.weights)] <- 9999
  simulated.data <- data.frame(X, exposure, Y1, Y0, Y.observed, true.weights,
                               propensity.scores, idx.p)
  true.sdsn <- subset(
    svydesign(id = ~ 1, weights = ~ true.weights, data = simulated.data),
    propensity.scores > 0 & propensity.scores < 1)
  true.sglm <- svyglm(Y.observed ~ exposure, design = true.sdsn,
                      family = quasibinomial())
  true.parameters <- coef(true.sglm)[2]
  return(list(simulated.data = simulated.data, true.parameters = true.parameters))
}

# Bias is the difference between sample estimate with estimated propensity score
# and estimate with true propensity score.
col.names <- c("tree_0", "tree_0.025", "tree_0.05", "tree_0.1",
               # Matching method.
               "M.near.logit", "M.cem.logit", "M.full.logit",
               "M.near.rpart", "M.cem.rpart", "M.full.rpart",
               # Generalized linear model.
               "IPW.glm", "IPW.glm_0.025", "IPW.glm_0.05", "IPW.glm_0.1",
               # Generalized boosted model.
               "IPW.gbm", "IPW.gbm_0.025", "IPW.gbm_0.05", "IPW.gbm_0.1",
               # Random forest.
               "IPW.rf", "IPW.rf_0.025", "IPW.rf_0.05", "IPW.rf_0.1")
bias.list <- vector("list", length = N.SIMULATIONS)
bias <- lapply(bias.list, function(bias.list) {
  matrix(NA, N.SAMPLES, length(col.names),
         dimnames = list(1:N.SAMPLES, col.names))
  }
)

for (s in 1:N.SIMULATIONS) {
  for (b in 1:N.SAMPLES) {
    set.seed(1234)
    current.data <- GenerateData(n.obs = N.OBS, n.covariates = N.COVARIATES[s],
                                 n.npositivity = N.NPOSITIVITY[s])
    data <- current.data$simulated.data
    parameters <- current.data$true.parameters
    exposure.formula <- ComputeExposureFormula(N.COVARIATES[s])

    # JCART: The propensity score model is estimated using CART. The prediction
    # model is based on Jackknife subsamples from the data set.
    n.strata <- N.STRATA
    n.cuts <- N.CUTS
    random.data <- sample(x = 1:N.OBS, replace = FALSE, size = N.OBS)
    idx.strata <- rep(1:n.strata, each = as.integer(N.OBS / n.strata))

    estimates <- matrix(NA, n.strata, n.cuts)
    for(j in 1:n.strata) {
      data.current.stratum <- data[random.data[idx.strata != j], ]

      # Perform recurive partitioning and regression tree for the current
      # stratum.
      fitted.rpart <- rpart(formula = exposure.formula,
                            data = data.current.stratum,
                            control = rpart.control(xval = N.VALIDATIONS))
      current.p.scores <- predict(fitted.rpart, newdata = data.current.stratum)
      data.current.stratum$where <- fitted.rpart$where

      for(j2 in 1:n.cuts) {
        if (j2 == 1) {
          trimmed.rule <- (current.p.scores > 0 & current.p.scores < 1)
        } else if (j2 == 2) {
          trimmed.rule <- (current.p.scores >= 0.025 & current.p.scores <= 0.975)
        } else if (j2 == 3) {
          trimmed.rule <- (current.p.scores >= 0.05 & current.p.scores <= 0.95)
        } else if (j2 == 4) {
          trimmed.rule <- (current.p.scores >= 0.1 & current.p.scores <= 0.9)
        }
        trimmed.data.current.stratum <- data.current.stratum[trimmed.rule, ]

        getY01Samples <- function(data){
          fitted.glm <- glm(outcome.formula, family = binomial, data = data)
          betas <- coef(fitted.glm)
          Y1 <- expit(sum(betas))
          Y0 <- expit(betas[1])
          Y01 <- c(Y0, Y1)
          names(Y01) <- c("Y0", "Y1")
          return(Y01)
        }

        # Causal effect estimates for each of the strata.
        Y01.list <- lapply(
          split(x = trimmed.data.current.stratum,
                f = as.factor(trimmed.data.current.stratum$where),
                drop = TRUE),
          getY01Samples)
        n.current.obs <- nrow(trimmed.data.current.stratum)
        where.prop <- table(trimmed.data.current.stratum$where) / n.current.obs
        n.where.prop <- length(where.prop)
        probabilties <- apply(
          do.call(cbind, Y01.list) *
            matrix(where.prop, 2, n.where.prop, byrow = TRUE), 1, sum)
        Y0.probs <- probabilties["Y0"]
        Y1.probs <- probabilties["Y1"]
        estimates[j, j2] <-
          (log(Y1.probs / (1 - Y1.probs) / (Y0.probs / (1 - Y0.probs))))
      }
    }
    bias[[s]][b, "tree_0"] <- mean(estimates[, 1]) - parameters
    bias[[s]][b, "tree_0.025"] <- mean(estimates[, 2]) - parameters
    bias[[s]][b, "tree_0.05"] <- mean(estimates[, 3]) - parameters
    bias[[s]][b, "tree_0.1"] <- mean(estimates[, 4]) - parameters

    # Matching methods: The propensity scores are estimated using matching
    # methods with distance criterion determined by "nearest-neighbor",
    # "coarsened exact matching", and "full matching". The logistic regression
    # model ("logit") and random forest ("rpart") are used in the propensity
    # score model.
    match.data1 <- matchit(exposure.formula, data = data, distance = "logit",
                           method = "nearest", discard = "both")
    fitted.glm1 <- glm(outcome.formula, data = match.data(match.data1),
                       family = binomial)
    bias[[s]][b, "M.near.logit"] <- coef(fitted.glm1)[2] - parameters

    match.data2 <- matchit(exposure.formula, data = data, distance = "rpart",
                           method = "nearest", discard = "both")
    fitted.glm2 <- glm(outcome.formula, data = match.data(match.data2),
                       family = binomial)
    bias[[s]][b, "M.near.rpart"] <- coef(fitted.glm2)[2] - parameters

    is.Error <- tryCatch(
      match.data3 <- matchit(exposure.formula, data = data, distance = "logit",
                             method = "cem", discard = "both"),
      error = function(e) e
    )
    if(!inherits(is.Error, "error")) {
      fitted.glm3 <- glm(outcome.formula, data = match.data(match.data3),
                         family = binomial)
      bias[[s]][b, "M.cem.logit"] <- coef(fitted.glm3)[2] - parameters
    }

    is.Error <- tryCatch(
      match.data4 <- matchit(exposure.formula, data = data, distance = "rpart",
                             method = "cem", discard = "both"),
      error = function(e) e
    )
    if(!inherits(is.Error, "error")) {
      fitted.glm4 <- glm(outcome.formula, data = match.data(match.data4),
                         family = binomial)
      bias[[s]][b, "M.cem.rpart"] <- coef(fitted.glm4)[2] - parameters
    }

    match.data5 <- matchit(exposure.formula, data = data, distance = "logit",
                           method = "full", discard = "both")
    dsn.match5 <- svydesign(id = ~ 1, weights = ~ weights,
                            data = match.data(match.data5))
    fitted.glm5 <- svyglm(outcome.formula, design = dsn.match5,
                          family = quasibinomial())
    bias[[s]][b, "M.full.logit"] <- coef(fitted.glm5)[2] - parameters

    match.data6 <- matchit(exposure.formula, data = data, distance = "rpart",
                           method = "full", discard = "both")
    dsn.match6 <- svydesign(id = ~ 1, weights = ~ weights,
                            data = match.data(match.data6))
    fitted.glm6 <- svyglm(outcome.formula, design = dsn.match6,
                          family = quasibinomial())
    bias[[s]][b, "M.full.rpart"] <- coef(fitted.glm6)[2] - parameters
    #
    # # The inverse of propensity weighting method. Estimations without removing
    # # the extreme group results in bias.
    #
    # # Propensity scores are estimated using logistic regression. The estimation
    # # model uses the IPW method to estimate causal effects.
    # data$glm.p.scores <- glm(exposure.formula, data = data)$fitted
    # data$glm.weights <- ComputeStabilizedWeights(data$exposure, data$glm.p.scores)
    # data$glm.weights[is.na(data$glm.weights) | data$glm.weights < 0] <- 9999
    #
    # dsn.glm <- subset(svydesign(id = ~1,weights = ~glm.weights, data = data),
    #                   data$glm.weights != 9999)
    # sglm.glm <- svyglm(outcome.formula, design = dsn.glm, family = quasibinomial())
    # bias[[s]][b, "IPW.glm"] <- coef(sglm.glm)[2]- parameters
    #
    # dsn.glm1 <- subset(svydesign(id = ~1,weights = ~glm.weights, data = data),
    #                    data$glm.weights != 9999 & data$glm.weights >= 0.025 & data$glm.weights <= 0.975)
    # sglm.glm1 <- svyglm(outcome.formula, design = dsn.glm1, family = quasibinomial())
    # bias[[s]][b, "IPW.glm_0.025"] <- coef(sglm.glm1)[2]- parameters
    #
    # dsn.glm2 <- subset(svydesign(id = ~1, weights = ~glm.weights, data = data),
    #                    data$glm.weights != 9999 & data$glm.weights >= 0.05 & data$glm.weights <= 0.95)
    # sglm.glm2 <- svyglm(outcome.formula, design = dsn.glm2, family = quasibinomial())
    # bias[[s]][b, "IPW.glm_0.05"] <- coef(sglm.glm2)[2]- parameters
    #
    # dsn.glm3 <- subset(svydesign(id = ~1, weights = ~glm.weights, data = data),
    #                    data$glm.weights != 9999 & data$glm.weights >= 0.1 & data$glm.weights <= 0.9)
    # sglm.glm3 <- svyglm(outcome.formula, design = dsn.glm3, family = quasibinomial())
    # bias[[s]][b, "IPW.glm_0.1"] <- coef(sglm.glm3)[2]- parameters
    #
    # # Propensity scores are estimated using GBM. The estimaton model uses the
    # # IPW method to estimate causal effects.
    # fitted.gbm <- gbm(exposure.formula, data = data, interaction.depth = 2,
    #                   n.trees = 3000, distribution = "bernoulli", verbose = FALSE)
    # data$gbm.p.scores <- predict(fitted.gbm, newdata = data, type = "response",
    #                              n.trees = 3000)
    # data$gbm.weights <- ComputeStabilizedWeights(data$exposure, data$gbm.p.scores)
    # data$gbm.weights[is.na(data$gbm.weights) | data$gbm.weights < 0] <- 9999
    #
    # dsn.gbm <- subset(svydesign(id = ~1,weights = ~gbm.weights, data = data),
    #                   data$gbm.weights != 9999)
    # sglm.gbm <- svyglm(outcome.formula, design = dsn.gbm, family = quasibinomial())
    # bias[[s]][b, "IPW.gbm"] <- coef(sglm.gbm)[2]- parameters
    #
    # dsn.gbm1 <- subset(svydesign(id = ~1,weights = ~gbm.weights, data = data),
    #                    data$gbm.weights != 9999 & data$gbm.weights >= 0.025 & data$gbm.weights <= 0.975)
    # sglm.gbm1 <- svyglm(outcome.formula, design = dsn.gbm1, family = quasibinomial())
    # bias[[s]][b, "IPW.gbm_0.025"] <- coef(sglm.gbm1)[2]- parameters
    #
    # dsn.gbm2 <- subset(svydesign(id = ~1,weights = ~gbm.weights, data = data),
    #                    data$gbm.weights != 9999 & data$gbm.weights >= 0.05 & data$gbm.weights <= 0.95)
    # sglm.gbm2 <- svyglm(outcome.formula, design = dsn.gbm2, family = quasibinomial())
    # bias[[s]][b, "IPW.gbm_0.05"] <- coef(sglm.gbm2)[2]- parameters
    #
    # dsn.gbm3 <- subset(svydesign(id = ~1,weights = ~gbm.weights, data = data),
    #                    data$gbm.weights != 9999 & data$gbm.weights >= 0.1 & data$gbm.weights <= 0.9)
    # sglm.gbm3 <- svyglm(outcome.formula, design = dsn.gbm3, family = quasibinomial())
    # bias[[s]][b, "IPW.gbm_0.1"] <- coef(sglm.gbm3)[2]- parameters
    #
    # # Propensity scores are estimated using random forest. The estimaton model
    # # uses the IPW method to estimate causal effects.
    # data$factored.exposure <- as.factor(data$exposure)
    # fexposure.formula <- as.formula(
    #   paste("factored.exposure ~ 0 +",
    #         paste(paste("X", 1:N.COVARIATES[s], sep = ""), collapse = "+"),
    #         sep = ""))
    # fitted.rf <- randomForest(fexposure.formula, data = data, ntree = 3000,
    #                           importance = TRUE, proximity = TRUE)
    # data$rf.p.scores <- predict(fitted.rf, newdata = data, type = "prob")[, "1"]
    # data$rf.weights <- ComputeStabilizedWeights(data$exposure, data$rf.p.scores)
    # data$rf.weights[is.na(data$rf.weights) | data$rf.weights < 0] <- 9999
    #
    # dsn.rf <- subset(svydesign(id = ~1, weights = ~rf.weights, data = data),
    #                  data$rf.weights != 9999)
    # sglm.rf<- svyglm(outcome.formula, design = dsn.rf, family = quasibinomial())
    # bias[[s]][b,"IPW.rf"]  <- coef(sglm.rf)[2] - parameters
  }
}

results <- vector("list", length = N.SIMULATIONS)
for (i in 1:N.SIMULATIONS) {
  results[[i]] <- data.frame(
    bias.mean = round(apply(bias[[i]], 2, mean), 2),
    bias.median = round(apply(bias[[i]], 2, median), 2),
    bias.rmse = round(apply(bias[[i]], 2, function(col) sqrt(mean(col ^ 2))), 2)
  )
}
