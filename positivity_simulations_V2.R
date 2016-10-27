# This script is a modified version of the R script provided in Petersen et al.
# (2012). This conducts simulation studies to investigate the performance of the
# parametric bootstrap as a diagnostic tool for assessing the positivity
# assumption.

suppressMessages(library(assertthat))

GenerateDataForSim1 <- function(n.obs) {
# Simulation 1: Generate a sample of size n.
# X = X1, X2, bivariate normal, mu1 = 0.5, mu2 = 1, sigma = (2, 1, 1, 1).
# Z = PHI(0.3* (0.5 + 0.25 * X1 + 0.75 * X2)) (probit model).
# Y = 1 + Z + X1 + 2 * X2 + epsilon, epsilon ~ N(0,1).
# true.effect = 1.
	sigma <- matrix(c(2, 1, 1, 1), ncol = 2)
	covariates <- matrix(rnorm(n.obs * 2), ncol = nrow(sigma)) %*% chol(sigma)
  covariates <- (covariates +
    matrix(rep(c(0.5, 1), each = n.obs), byrow = FALSE, ncol = 2))
	colnames(covariates) <- c("X1", "X2")
	exposure <- (0.2 * (0.5 + 0.25 * covariates[, "X1"]
	                    + 0.75 * covariates[, "X2"]) + rnorm(n.obs) > 0)
  exposure <- as.integer(exposure)
	outcome <- (1 + exposure + covariates[,"X1"] + 2 * covariates[,"X2"] +
	  rnorm(n.obs))
	return(data.frame(outcome, exposure, covariates))
}

bound <- function(x, bounds){
  # Truncates values at (min, max) bounds.
	x[x < min(bounds)] <- min(bounds)
	x[x > max(bounds)] <- max(bounds)
	return(x)
}

tmle <- function(outcome, exposure, Q, p.scores, alpha){
  # Targeted maximum likelihood estimation of the additive treatment effect
  # based on external estimates of Q and g1W. Targeting through a logistic
  # fluctuation, with predicted probabilities bounded away from (0,1) at (alpha,
  # 1-alpha).
	outcome.range <- range(outcome)
	Q <- bound(Q, outcome.range)
	Q <- qlogis(bound((Q - outcome.range[1]) / diff(outcome.range),
	                  c(alpha, 1 - alpha)))
	new.outcome <- (outcome - outcome.range[1]) / diff(outcome.range)
	h <- cbind((exposure / p.scores - (1 - exposure) / (1 - p.scores)),
	           1 / p.scores, -1 / (1 - p.scores))
	suppressWarnings(
	  eps <- coef(glm(new.outcome ~ -1 + offset(Q[, 1]) + h[, 1],
	                  family = "binomial")))
	new.Q <- plogis(Q + eps * h) * diff(outcome.range) + outcome.range[1]
	effect <- mean(new.Q[,2]) - mean(new.Q[,3])
	return(effect)
}

ipwt <- function(outcome, exposure, Qfamily, weights){
  # Inverse probability of treatment weighting estimation of the additive
  # treatment effect using wts based on inverse of g1W.
	if (identical(Qfamily, "binomial") | identical(Qfamily, binomial)) {
		EY1 <- exposure * weights * outcome
		EY0 <- (1 - exposure)* weights * outcome
		Qmodel <- EY1 * exposure + EY0 * (1 - exposure)
		effect <- mean(EY1) - mean(EY0)
	} else {
		fitted.glm <- glm(outcome ~ exposure, weights = weights, family = Qfamily)
		Qmod <- predict(fitted.glm, type = "response")
		Qmod.with.exposure <- predict(fitted.glm, newdata = data.frame(exposure = 1))
		Qmod.no.exposure <- predict(fitted.glm, newdata = data.frame(exposure = 0))
		effect <- Qmod.with.exposure - Qmod.no.exposure
	}
	return(effect)
}

bias.pboot <- function(data, n.samples, gbounds = c(0, 1), Qmethod, gmethod,
                       Qformula, gformula, Qfamily = "gaussian",
                       gfamily = "binomial" (link = "probit"),
                       is.pboot = TRUE) {

	# Check argument inputs.
  assert_that(Qmethod %in% c("glm", "step1", "step2"))
  assert_that(gmethod %in% c("glm", "step"))
  assert_that(min(gbounds) >= 0 & max(gbounds) <= 1)
	if (length(gbounds) == 1){
	  gbounds <- c(gbounds, 1 - gbounds)
	}

  # Set outcome, exposure and covariate columns.
  outcome.col <- colnames(data) %in% "outcome"
  exposure.col <- colnames(data) %in% "exposure"
  covariate.cols <- !colnames(data) %in% c("outcome", "exposure")
  covariate.names <- paste(colnames(data)[covariate.cols], collapse = "+")

  EstimatePropensityModel <- function(g.formula) {
    if(gmethod == "glm"){
      fitted.gglm <- glm(g.formula, data = data, family = gfamily)
    } else {
      fitted.gglm <- glm("exposure ~ 1", data = data, family = gfamily)
      upper.scope <- paste("~", covariate.names)
      glm.scope <- list(upper = upper.scope, lower = ~ 1)
      fitted.gglm <- step(fitted.gglm, scope = glm.scope, trace = 0,
                          direction = "forward", data = data[, -outcome.col])
    }

    p.scores <- bound(predict(fitted.gglm, type = "response"), gbounds)
    weights <- 1 / p.scores
    weights[data$exposure == 0] <- 1 / (1 - p.scores[data$exposure == 0])

    return(list(fitted.gglm = fitted.gglm, p.scores = p.scores,
                weights = weights))
  }

  EstimateOutcomeModel <- function(Q.formula) {
    if (Qmethod == "glm") {
      fitted.Qglm <- glm(Q.formula, data = data, family = Qfamily)
    } else {
      Qformula <- "outcome ~ exposure"
      lower.scope <- "~ exposure"
      if (Qmethod == "step2" & length(unique(p.scores)) > 1) {
        Qformula <- "outcome ~ exposure + p.scores"
        lower.scope <- "~ exposure + p.scores"
      }
      fitted.Qglm <- glm(as.formula(Qformula),
                         data = data.frame(data, p.scores), family = Qfamily)
      upper.scope <- paste(lower.scope, "+", covariate.names)
      glm.scope <- list(upper = upper.scope, lower = lower.scope)
      fitted.Qglm <- step(fitted.Qglm, scope = glm.scope, trace = 0,
                          direction = "forward", data.frame(data, p.scores))
    }

    Qmod <- predict(fitted.Qglm)
    data.with.exposure <- data.frame(data[, -exposure.col], exposure = 1)
    Qmod.with.exposure <- predict(fitted.Qglm, newdata = data.with.exposure)
    data.no.exposure <- data.frame(data[, -exposure.col], exposure = 0)
    Qmod.no.exposure <- predict(fitted.Qglm, newdata = data.no.exposure)

    return(list(fitted.Qglm = fitted.Qglm, Qmod = Qmod,
                Qmod.with.exposure = Qmod.with.exposure,
                Qmod.no.exposure = Qmod.no.exposure))
  }

  EstimateViaParaBootstrap <- function(Qmodel, gmodel, covariates) {
    GenerateData <- function() {
      covariate.samples <- covariates[idx.boot.samples, ]
      probabilties <- predict(gmodel,
                              newdata = data.frame(covariate.samples),
                              type = "response")
      exposure <- rbinom(n.obs, 1, probabilties)
      test.data <- data.frame(exposure, covariate.samples, probabilties)
      outcome <- predict(Qmodel, newdata = test.data, type = "response") +
        rnorm(n.obs)
      data <- data.frame(outcome, exposure, covariate.samples)
      return(data)
    }

    n.obs <- nrow(covariates)
    idx.boot.samples <- sample(1:n.obs, n.obs, replace = TRUE)
    data <- GenerateData()

    # Estimate P(Z = 1 | X) on bootstrap samples.
    p.scores <- EstimatePropensityModel(gmodel$formula)$p.scores
    weights <- EstimatePropensityModel(gmodel$formula)$weights

    # Estimate P(Y = y | Z, X) based on bootstrap samples.
    Qmod <- EstimateOutcomeModel(Qmodel$formula)$Qmod
    Qmod.with.exposure <- EstimateOutcomeModel(Qmodel$formula)$Qmod.with.exposure
    Qmod.no.exposure <- EstimateOutcomeModel(Qmodel$formula)$Qmod.no.exposure

    # Get estimates on bootstrap samples.
    IPWT.effect <- ipwt(data$outcome, data$exposure, Qfamily, weights)
    TMLE.effect <- tmle(data$outcome, data$exposure,
                        Q = cbind(Qmod, Qmod.with.exposure, Qmod.no.exposure),
                        p.scores = p.scores, alpha = 0.99999)
    Gcomp.effect <- mean(Qmod.with.exposure) - mean(Qmod.no.exposure)
    AIPWT.effect <- mean((data$exposure - (1 - data$exposure)) * weights *
                           (data$outcome - Qmod) + Qmod.with.exposure -
                           Qmod.no.exposure)
    return(c(Gcomp = Gcomp.effect, IPWT = IPWT.effect, AIPWT = AIPWT.effect,
             TMLE = TMLE.effect))
  }

  # Estimate P(Z = 1 | X) on original sample.
  fitted.gglm <- EstimatePropensityModel(gformula)$fitted.gglm
  p.scores <- EstimatePropensityModel(gformula)$p.scores
  weights <- EstimatePropensityModel(gformula)$weights

  # Estimate P(Y = y | Z, X) based on original sample.
  fitted.Qglm <- EstimateOutcomeModel(Qformula)$fitted.Qglm
  Qmod <- EstimateOutcomeModel(Qformula)$Qmod
  Qmod.with.exposure <- EstimateOutcomeModel(Qformula)$Qmod.with.exposure
  Qmod.no.exposure <- EstimateOutcomeModel(Qformula)$Qmod.no.exposure

  # Get estimates on original sample.
  IPWT.effect <- ipwt(data$outcome, data$exposure, Qfamily, weights)
  TMLE.effect <- tmle(data$outcome, data$exposure,
                      Q = cbind(Qmod, Qmod.with.exposure, Qmod.no.exposure),
                      p.scores = p.scores,
                      alpha = 0.99999)
  Gcomp.effect <- mean(Qmod.with.exposure) - mean(Qmod.no.exposure)
  AIPWT.effect <- mean((data$exposure - (1 - data$exposure)) * weights *
                         (data$outcome - Qmod) + Qmod.with.exposure -
                         Qmod.no.exposure)

  # Call pboot function to get estimates bootstrap samples.
  if (isTRUE(is.pboot)) {
    bootstrap.samples <- replicate(n.samples,
      EstimateViaParaBootstrap(Qmodel = fitted.Qglm, gmodel = fitted.gglm,
                              covariates = data[, covariate.cols, drop = FALSE])
    )
    mean.bootstrap.samples <- rowMeans(bootstrap.samples)
    ETA.bias <- rowMeans(bootstrap.samples - Gcomp.effect)
  } else {
    bootstrap.samples <- mean.bootstrap.samples <- ETA.bias <- NULL
  }
  return(list(estimated.effects= c(Gcomp = Gcomp.effect, IPWT = IPWT.effect,
                                   AIPWT = AIPWT.effect, TMLE = TMLE.effect),
              bootstrap.effects = mean.bootstrap.samples,
              ETA.bias = ETA.bias,
              bootstrap.samples = bootstrap.samples))
}
