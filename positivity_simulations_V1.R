# This script is a mordified version of the R script provided in Kang et al.
# (2016). This conducts simulation studies to assess the effect of arbitrary
# cutoffs of propensity scores.

# Load required libraries.
suppressMessages(library(cem))
suppressMessages(library(gbm))
suppressMessages(library(MatchIt))
suppressMessages(library(MASS))
suppressMessages(library(randomForest))
suppressMessages(library(rpart))
suppressMessages(library(survey))
suppressMessages(library(tree))

is.runJCART <- TRUE
is.runMatching <- TRUE
is.runIPW <- FALSE

# Set up dimension effect.
N.OBS <- 1000
N.COVARIATES <- c(4, 4, 40, 40)
N.NPOSITIVITY <- c(1, 1, 2, 2)
N.SIMULATIONS <- length(N.COVARIATES)
N.SAMPLES <- 200
N.VALIDS <- 40 # Number of cross-validations for JCART.
N.STRATA <- 4 # Number of strata/groups for performing JCART.
N.CUTS <- 4 # Number of propensity score trimming rules for JCART.

# Set up response-exposure formula.
outcome.formula <- "Y.observed ~ exposure"

# Set up exposure-confounder formula.
ComputeExposureFormula <- function(n.covariates) {
  exposure.formula <- as.formula(
    paste("exposure ~ 0 +",
          paste(paste("X", 1:n.covariates, sep = ""), collapse = "+"),
          sep = "")
  )
}

# Define general functions.
expit <- function(a) {
  exp(a) / (1 + exp(a))
}

ComputePropensityWeights <- function(exposure, p.scores) {
  first.term <- exposure * mean(exposure) / p.scores
  second.term <- (1 - exposure) * mean(1 - exposure) / (1 - p.scores)
  weights <- first.term + second.term
  return(weights)
}

CreateSurveyDesign <- function(weights, data, subset.rule) {
  subset(svydesign(id = ~ 1, weights = ~ weights, data = data),
         subset.rule)
}

GenerateData <- function(n.obs, n.covariates, n.npositivity) {
  # Generates data for four simulation study designs (Sim1 to Sim4).

  # Generate categorical confounders, which were randomly sampled with equal
  # probability from the interval [0.1, 1].
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
    idx.np.upp <- (x12 >= quantile(x12, probs = 0.9))
    idx.p <- (x12 > quantile(x12, probs = 0.1) &
                x12 < quantile(x12, probs = 0.9))
  }

  # Generate the true propensity scores: e(X) = P(Z = 1 | X) = expit(X * beta).
  beta <- c(rep(1, (n.covariates / 2)), rep(-1, (n.covariates / 2)))
  eta <- X[idx.p, ] %*% beta
  true.p.scores <- rep(NA, n.obs)
  true.p.scores[idx.p] <- exp(eta) / (1 + exp(eta))

  # Set propensity scores to 0 or 1 for those non-positivity subjects.
  true.p.scores[idx.np.low] <- 0
  true.p.scores[idx.np.upp] <- 1

  # Generate the binary exposure variable.
  exposure <- rep(NA, n.obs)
  exposure <- 1 * (true.p.scores > runif(n.obs))

  # Propensity structure for the simulation study.
  Y1 <- rep(NA, n.obs)
  Y0 <- rep(NA, n.obs)

  # For e(x) = 0. If the covariate X1 had values 0.1 or 1, the potential outcome
  # Y(1) for these subjects would be sampled with fixed probability 0.1 whereas
  # if the covariate X1 had values between 0.1 and 1, the potential outcome Y(1)
  # would be sam,pled with probabillity 0.6.
  Y1[idx.np.low] <-  rbinom(n = sum(idx.np.low), size = 1, prob = 0.1)
  Y0[idx.np.low] <-  rbinom(n = sum(idx.np.low), size = 1, prob = 0.9)

  # For e(x) = 1.
  Y1[idx.np.upp] <-  rbinom(n = sum(idx.np.upp), size = 1, prob = 0.1)
  Y0[idx.np.upp] <-  rbinom(n = sum(idx.np.upp), size = 1, prob = 0.9)

  # For 0 < e(x) < 1.
  Y1[idx.p] <-  rbinom(n = sum(idx.p), size = 1, prob = 0.6)
  Y0[idx.p] <-  rbinom(n = sum(idx.p), size = 1, prob = 0.4)

  #The observed outcome Y was defined as Y = Z * Y(1) + (1 − Z) * Y(0).
  Y.observed <- exposure * Y1 + (1 - exposure) * Y0
  true.weights <- ComputePropensityWeights(exposure, true.p.scores)
  true.weights[is.na(true.weights)] <- 9999

  simulated.data <- data.frame(idx.p, Y.observed, exposure, X, Y1, Y0,
                               true.p.scores, true.weights)

  # Note that it is impossible to estimate E(Y(1)) for subjects with e(x) = 0
  # and so they should be removed from the study. Similarly, subjects with e(x)
  # = 1 should also be removed.
  true.sdsn <- CreateSurveyDesign(
    weights = true.weights,
    data = simulated.data,
    subset.rule = true.p.scores > 0 & true.p.scores < 1)
  true.sglm <- svyglm(outcome.formula, design = true.sdsn,
                      family = quasibinomial())
  true.effect <- coef(true.sglm)[2]
  return(list(simulated.data = simulated.data, true.effect = true.effect))
}

# Initalize an empty list for the bias matrices for each simulation. Bias is
# defined as the difference between sample estimate with estimated propensity
# score and estimate with true propensity score.
col.names <- c("tree_0", "tree_0.025", "tree_0.05", "tree_0.1",
               "M.near.logit", "M.cem.logit", "M.full.logit",
               "M.near.rpart", "M.cem.rpart", "M.full.rpart",
               "IPW.glm", "IPW.glm_0.025", "IPW.glm_0.05", "IPW.glm_0.1",
               "IPW.gbm", "IPW.gbm_0.025", "IPW.gbm_0.05", "IPW.gbm_0.1",
               "IPW.rf", "IPW.rf_0.025", "IPW.rf_0.05", "IPW.rf_0.1")
bias.list <- vector("list", length = N.SIMULATIONS)
bias <- lapply(bias.list, function(bias.list) {
  matrix(NA, N.SAMPLES, length(col.names),
         dimnames = list(1:N.SAMPLES, col.names))
  }
)

for (s in 1:N.SIMULATIONS) {
  set.seed(1234)
  exposure.formula <- ComputeExposureFormula(N.COVARIATES[s])

  for (b in 1:N.SAMPLES) {
    current.data <- GenerateData(n.obs = N.OBS, n.covariates = N.COVARIATES[s],
                                 n.npositivity = N.NPOSITIVITY[s])
    data <- current.data$simulated.data
    effect <- current.data$true.effect

    if (isTRUE(is.runJCART)) {
      # JCART: The propensity score model is estimated using CART. The
      # prediction model is based on Jackknife subsamples from the data set.
      n.strata <- N.STRATA
      n.cuts <- N.CUTS

      # Divide the data into "n.strata" groups randomly. Each group has (1 /
      # n.strata) * 100% of the original data sample.
      random.data <- sample(x = 1:N.OBS, replace = FALSE, size = N.OBS)
      idx.strata <- rep(1:n.strata, each = as.integer(N.OBS / n.strata))

      estimates <- matrix(NA, n.strata, n.cuts)
      for(j in 1:n.strata) {
        # Delete group j out of "n.strata" randomly divided groups.
        data.current.stratum <- data[random.data[idx.strata != j], ]

        # Use the rest of the "n.strata - 1" groups for growing the CART with
        # multiple folds cross-validation.
        fitted.rpart <- rpart(formula = exposure.formula,
                              data = data.current.stratum,
                              control = rpart.control(xval = N.VALIDS))
        current.p.scores <- predict(fitted.rpart,
                                    newdata = data.current.stratum)

        # "where" is an integer vector (N.Obs) that contains the row number of
        # data frame corresponding to the leaf node that each observations falls
        # into.
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

          # The causal estimator of L j is given by:
          #              logit(E(Ysj (1))) − logit(E(Ysj (0))),
          # and E(Ysj (1)) is estimated by Sum_sj E(Ysj | zsj = 1) P(sj), where
          # E(Ysj | zsj ) can be simply modeled via glm.
          Y01.list <- lapply(
            split(x = trimmed.data.current.stratum,
                  f = as.factor(trimmed.data.current.stratum$where),
                  drop = TRUE),
            getY01Samples)
          n.curr.obs <- nrow(trimmed.data.current.stratum)
          where.prop <- table(trimmed.data.current.stratum$where) / n.curr.obs
          where.prop.matrix <- matrix(where.prop, 2, length(where.prop),
                                      byrow = TRUE)
          probabilties <- apply(do.call(cbind, Y01.list) *  where.prop.matrix,
                                1, sum)
          Y0.probs <- probabilties["Y0"]
          Y1.probs <- probabilties["Y1"]
          estimates[j, j2] <-
            (log(Y1.probs / (1 - Y1.probs) / (Y0.probs / (1 - Y0.probs))))
        }
      }
      bias[[s]][b, "tree_0"] <- mean(estimates[, 1]) - effect
      bias[[s]][b, "tree_0.025"] <- mean(estimates[, 2]) - effect
      bias[[s]][b, "tree_0.05"] <- mean(estimates[, 3]) - effect
      bias[[s]][b, "tree_0.1"] <- mean(estimates[, 4]) - effect
    }

    if (isTRUE(is.runMatching)) {
      # Matching methods: The propensity scores are estimated using matching
      # methods with distance criterion determined by "nearest-neighbor",
      # "coarsened exact matching", and "full matching". The logistic regression
      # model ("logit") and random forest ("rpart") are used in the propensity
      # score model.
      PropensityViaMatching <- function(method, distance) {
        is.Error <- tryCatch(
          match.data <- matchit(exposure.formula,
                                data = data,
                                distance = distance,
                                method = method,
                                discard = "both"),
          error = function(e) e
        )

        if(!inherits(is.Error, "error")) {
          fitted.glm <- glm(outcome.formula, data = match.data(match.data),
                            family = binomial)
          bias.estimate <- coef(fitted.glm)[2] - effect
        }
        return(bias.estimate)
      }

      # Note: "rpart" stands for random forest, need to set a random seed for
      # reproducibility.
      bias[[s]][b, "M.near.logit"] <- PropensityViaMatching("nearest", "logit")
      bias[[s]][b, "M.near.rpart"] <- PropensityViaMatching("nearest", "rpart")
      bias[[s]][b, "M.cem.logit"] <- PropensityViaMatching("cem", "logit")
      bias[[s]][b, "M.cem.rpart"] <- PropensityViaMatching("cem", "rpart")
      bias[[s]][b, "M.full.logit"] <- PropensityViaMatching("full", "logit")
      bias[[s]][b, "M.full.rpart"] <- PropensityViaMatching("full", "rpart")
    }

    if (isTRUE(is.runIPW)) {
      # Propensity scores are estimated using logistic regression. The
      # estimation model uses the IPW method to estimate causal effects.
      data$glm.p.scores <- glm(exposure.formula, data = data)$fitted
      data$glm.weights <- ComputePropensityWeights(data$exposure,
                                                   data$glm.p.scores)
      data$glm.weights[is.na(data$glm.weights) | data$glm.weights < 0] <- 9999

      PropensityViaIPWglm <- function(low.cut, high.cut) {
        dsn.glm <- CreateSurveyDesign(
          weights = data$glm.weights, data = data,
          subset.rule = (data$glm.weights != 9999 &
                           data$glm.weights >= low.cut &
                           data$glm.weights <= high.cut))
        sglm.glm <- svyglm(outcome.formula, design = dsn.glm,
                           family = quasibinomial())
        bias.estimate <- coef(sglm.glm)[2]- effect
        return(bias.estimate)
      }

      bias[[s]][b, "IPW.glm"] <-
        PropensityViaIPWglm(low.cut = -99, high.cut = 99)
      bias[[s]][b, "IPW.glm_0.025"] <-
        PropensityViaIPWglm(low.cut = 0.025, high.cut = 0.975)
      bias[[s]][b, "IPW.glm_0.05"] <-
        PropensityViaIPWglm(low.cut = 0.05, high.cut = 0.95)
      bias[[s]][b, "IPW.glm_0.1"] <-
        PropensityViaIPWglm(low.cut = 0.1, high.cut = 0.9)

      # Propensity scores are estimated using generalized boosted regression.
      # The estimaton model uses the IPW method to estimate causal effects.
      fitted.gbm <- gbm(exposure.formula, data = data, interaction.depth = 2,
                        n.trees = 10 ^ 3, distribution = "bernoulli",
                        verbose = FALSE)
      data$gbm.p.scores <- predict(fitted.gbm, newdata = data,
                                   type = "response", n.trees = 10 ^ 3)
      data$gbm.weights <- ComputePropensityWeights(data$exposure,
                                                   data$gbm.p.scores)
      data$gbm.weights[is.na(data$gbm.weights) | data$gbm.weights < 0] <- 9999

      PropensityViaIPWgbm <- function(low.cut, high.cut) {
        dsn.gbm <- CreateSurveyDesign(
          weights = data$gbm.weights, data = data,
          subset.rule = (data$gbm.weights != 9999 &
                           data$gbm.weights >= low.cut &
                           data$gbm.weights <= high.cut))
        sglm.gbm <- svyglm(outcome.formula, design = dsn.gbm,
                           family = quasibinomial())
        bias.estimate <- coef(sglm.gbm)[2]- effect
        return(bias.estimate)
      }

      bias[[s]][b, "IPW.gbm"] <-
        PropensityViaIPWgbm(low.cut = -99, high.cut = 99)
      bias[[s]][b, "IPW.gbm_0.025"] <-
        PropensityViaIPWgbm(low.cut = 0.025, high.cut = 0.975)
      bias[[s]][b, "IPW.gbm_0.05"] <-
        PropensityViaIPWgbm(low.cut = 0.05, high.cut = 0.95)
      bias[[s]][b, "IPW.gbm_0.1"] <-
        PropensityViaIPWgbm(low.cut = 0.1, high.cut = 0.9)

      # Propensity scores are estimated using random forest. The estimaton model
      # uses the IPW method to estimate causal effects.
      data$factored.exposure <- as.factor(data$exposure)
      fexposure.formula <- as.formula(
        paste("factored.exposure ~ 0 +",
              paste(paste("X", 1:N.COVARIATES[s], sep = ""), collapse = "+"),
              sep = ""))
      fitted.rf <- randomForest(fexposure.formula, data = data, ntree = 3000,
                                importance = TRUE, proximity = TRUE)
      data$rf.p.scores <- predict(fitted.rf, newdata = data, type = "prob")[, "1"]
      data$rf.weights <- ComputePropensityWeights(data$exposure,
                                                  data$rf.p.scores)
      data$rf.weights[is.na(data$rf.weights) | data$rf.weights < 0] <- 9999

      PropensityViaIPWrf <- function(low.cut, high.cut) {
        dsn.rf <- CreateSurveyDesign(
          weights = data$rf.weights, data = data,
          subset.rule = (data$rf.weights != 9999 &
                           data$rf.weights >= low.cut &
                           data$rf.weights <= high.cut))
        sglm.rf <- svyglm(outcome.formula, design = dsn.rf,
                          family = quasibinomial())
        bias.estimate <- coef(sglm.rf)[2]- effect
        return(bias.estimate)
      }
    }
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
results
