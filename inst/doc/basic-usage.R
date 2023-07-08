## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(survML)
library(ggplot2)

## ----stackG_example-----------------------------------------------------------
# This is a small simulation example
set.seed(123)
n <- 500
X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))

S0 <- function(t, x){
  pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), lower.tail = FALSE)
}
T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))

G0 <- function(t, x) {
  as.numeric(t < 15) *.9*pexp(t,
                              rate = exp(-2 -.5*x[,1]-.25*x[,2]+.5*x[,1]*x[,2]),
                              lower.tail=FALSE)
}
C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
C[C > 15] <- 15

time <- pmin(T, C)
event <- as.numeric(T <= C)

# note that this a very small library, just for demonstration
SL.library <- c("SL.mean", "SL.glm", "SL.gam")

fit <- stackG(time = time,
              event = event,
              X = X,
              newX = X,
              newtimes = seq(0, 15, .1),
              direction = "prospective",
              bin_size = 0.02,
              time_basis = "continuous",
              time_grid_approx = sort(unique(time)),
              surv_form = "exp",
              SL_control = list(SL.library = SL.library,
                                V = 5))

## ----plot_stackG_example------------------------------------------------------
plot_dat <- data.frame(fitted = fit$S_T_preds[1,], 
                       true = S0(t =  seq(0, 15, .1), X[1,]))

p <- ggplot(data = plot_dat, mapping = aes(x = true, y = fitted)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  theme_bw() + 
  ylab("fitted") +
  xlab("true") + 
  ggtitle("Global survival stacking example (event time distribution)")

p

## ----plot_stackG_example_cens-------------------------------------------------
plot_dat <- data.frame(fitted = fit$S_C_preds[1,], 
                       true = G0(t =  seq(0, 15, .1), X[1,]))

p <- ggplot(data = plot_dat, mapping = aes(x = true, y = fitted)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  theme_bw() + 
  ylab("fitted") +
  xlab("true") + 
  ggtitle("Global survival stacking example (censoring time distribution)")

p

## ----stackL_example-----------------------------------------------------------
fit <- stackL(time = time,
              event = event,
              X = X,
              newX = X,
              newtimes = seq(0, 15, .1),
              direction = "prospective",
              bin_size = 0.02,
              time_basis = "continuous",
              SL_control = list(SL.library = SL.library,
                                V = 5))

## ----plot_stackL_example------------------------------------------------------
plot_dat <- data.frame(fitted = fit$S_T_preds[1,], 
                       true = S0(t =  seq(0, 15, .1), X[1,]))

p <- ggplot(data = plot_dat, mapping = aes(x = true, y = fitted)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  theme_bw() + 
  ylab("fitted") +
  xlab("true") + 
  ggtitle("Local survival stacking example")

p

