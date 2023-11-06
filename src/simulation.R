library(deSolve)
library(tidyverse)
library(rootSolve)

theme_set(theme_minimal())

seirmod <- function(t, y, parms) {
  # Pull state variables from y vector
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]

  # Pull parameter values from parms vector
  beta <- parms["beta"]
  sigma <- parms["sigma"]
  mu <- parms["mu"]
  gamma <- parms["gamma"]
  N <- parms["N"]

  # Define equations
  dS <- mu * (N - S) - beta * S * I / N
  dE <- beta * S * I / N - sigma * E
  dI <- sigma * E - (mu + gamma) * I
  dR <- gamma * I - mu * R
  res <- c(dS, dE, dI, dR)

  # Return list of gradients
  list(res)
}


times <- seq(0, 26, by = 1 / 10)
parms <- c(mu = 0, N = 1, beta = 2, sigma = 1, gamma = 1 / 2)
start <- c(S = 0.999, E = 0.0, I = 0.001, R = 0)

out <- ode(y = start, times = times, func = seirmod, parms = parms)
out_df <- as_tibble(out) %>%
  pivot_longer(cols = -time, names_to = "state", values_to = "number") %>%
  mutate(
    time = as.numeric(time),
    number = as.numeric(number),
    state = factor(state, levels = c("S", "E", "I", "R")),
    number = round(number, 6)
  )

ggplot(out_df, aes(x = time, y = number, color = state)) +
  geom_line(linewidth = 2) +
  labs(x = "Time", y = "Number", color = "State")


# Candidate values for R0 and beta
R0 <- seq(0.1, 5, length = 50)
betas <- R0 * 1 / 2

# Calculate proportion infected for each value of R0
# map2_dfr is a {purrr} function that applies a function to two vectors i.e., it is a vectorized version of a for loop, and returns a data frame
final_size_df <- map2_dfr(
  .x = betas,
  .y = R0,
  .f = function(.x, .y) {
    equil <- runsteady(
      y = c(S = 1 - 1E-5, E = 0.0, I = 1E-5, R = 0),
      times = c(0, 1E5),
      func = seirmod,
      parms = c(mu = 0, N = 1, beta = .x, sigma = 1, gamma = 1 / 2)
    )

    tibble(
      R0 = .y,
      final_size = equil$y["R"]
    )
  }
)

ggplot(final_size_df, aes(x = R0, y = final_size)) +
  geom_line(linewidth = 2) +
  labs(x = "R0", y = "Final size")

