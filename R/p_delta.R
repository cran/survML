#' Wrapper for various p_delta algorithms
#'
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return An fitted binary regression for (complement of)
#' probability of censoring
#'
#' @noRd
p_delta <- function(event,
                    X,
                    learner = "SuperLearner",
                    SL_control,
                    xgb_control){
  if (learner == "SuperLearner"){
    fit <- p_delta_SuperLearner(event = event,
                                X = X,
                                SL_control = SL_control)
  }

  return(fit)
}
