## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(survML)
library(ggplot2)
library(dplyr)
set.seed(102524)

## -----------------------------------------------------------------------------
# Simulate some current status data
n <- 250
x <- cbind(2*rbinom(n, size = 1, prob = 0.5)-1,
           2*rbinom(n, size = 1, prob = 0.5)-1)
t <- rweibull(n,
              shape = 0.75,
              scale = exp(0.8*x[,1] - 0.4*x[,2]))
y <- rweibull(n,
              shape = 0.75,
              scale = exp(0.8*x[,1] - 0.4*x[,2]))

# Round y to nearest quantile of y, just so there aren't so many unique values
# This will speed computation in this example analysis
quants <- quantile(y, probs = seq(0, 1, by = 0.025), type = 1)
for (i in 1:length(y)){
  y[i] <- quants[which.min(abs(y[i] - quants))]
}
delta <- as.numeric(t <= y)

dat <- data.frame(y = y, delta = delta, x1 = x[,1], x2 = x[,2])

dat$delta[dat$y > 1.65] <- NA
dat$y[dat$y > 1.65] <- NA

## -----------------------------------------------------------------------------
eval_region <- c(0.02, 1.5)
res <- currstatCIR(time = dat$y,
                   event = dat$delta,
                   X = dat[,3:4],
                   SL_control = list(SL.library = c("SL.mean", "SL.glm"),
                                     V = 2),
                   HAL_control = list(n_bins = c(5),
                                      grid_type = c("equal_mass", "equal_range"),
                                      V = 2),
                   eval_region = eval_region,
                   n_eval_pts = 1000)


## -----------------------------------------------------------------------------
# use Monte Carlo to approximate the true survival function
n_test <- 5e5
x_test <- cbind(2*rbinom(n_test, size = 1, prob = 0.5)-1,
                2*rbinom(n_test, size = 1, prob = 0.5)-1)
t_test <- rweibull(n_test,
                   shape = 0.75,
                   scale = exp(0.8*x_test[,1] - 0.4*x_test[,2]))

S0 <- function(x){
  return(mean(t_test > x))
}

other_data <- data.frame(t = seq(min(res$t), max(res$t), length.out = 1000))
other_data$y <- apply(as.matrix(other_data$t), MARGIN = 1, FUN = S0)

# plot the results
p1 <- ggplot(data = res, aes(x = t)) +
  geom_step(aes(y = S_hat_est)) +
  geom_step(aes(y = S_hat_cil), linetype = "dashed") +
  geom_step(aes(y = S_hat_ciu), linetype = "dashed") +
  geom_smooth(data = other_data, aes(x = t, y = y), color = "red") +
  theme_bw() +
  ylab("Estimated survival probability") +
  xlab("Time") + 
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle("Covariate-adjusted survival curve")
p1

