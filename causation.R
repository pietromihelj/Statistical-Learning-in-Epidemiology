library(tidyverse)
library(ggplot2)

##
data <- data.frame(
  id = 1:10,
  y_chocolate = c(4, 4, 6, 5, 6, 5, 6, 7, 5, 6),
  y_vanilla = c(1, 3, 4, 5, 5, 6, 8, 6, 3, 5)
)

data <- data |>
  mutate(causal_effect = y_chocolate - y_vanilla)

data

##
data |>
  summarize(
    avg_chocolate = mean(y_chocolate),
    avg_vanilla = mean(y_vanilla),
    avg_causal_effect = mean(causal_effect)
  )

##
data_observed <- data |>
  mutate(
    exposure = case_when(
      # people who like chocolate more chose that
      y_chocolate > y_vanilla ~ "chocolate",
      # people who like vanilla more chose that
      y_vanilla >= y_chocolate ~ "vanilla"
    ),
    observed_outcome = case_when(
      exposure == "chocolate" ~ y_chocolate,
      exposure == "vanilla" ~ y_vanilla
    )
  ) |>
  # we can only observe the exposure and one potential outcome!
  select(id, exposure, observed_outcome)
data_observed

##
data_observed |>
  group_by(exposure) |>
  summarise(avg_outcome = mean(observed_outcome))

##
set.seed(11)
data_observed <- data |>
  mutate(
    # change the exposure to randomized, generate from a binomial distribution
    # with a probability 0.5 for being in either group
    exposure = case_when(
      rbinom(10, 1, 0.5) == 1 ~ "chocolate",
      TRUE ~ "vanilla"
    ),
    observed_outcome = case_when(
      exposure == "chocolate" ~ y_chocolate,
      exposure == "vanilla" ~ y_vanilla
    )
  ) |>
  # as before, we can only observe the exposure and one potential outcome
  select(id, exposure, observed_outcome)
data_observed |>
  group_by(exposure) |>
  summarise(avg_outcome = mean(observed_outcome) )
