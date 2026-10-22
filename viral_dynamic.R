#Load data
df <- read.csv("path/ATCG315.csv")

#install package 
library(nlme)
library(dplyr)
library(ggplot2)

##############  EDA ############################# 
head(df)

# Check column names
names(df)

# Check data types and structure
str(df)

# Check the number of rows and columns
dim(df)

# Summary statistics
summary(df)

# Detection limit on the log10 scale
LLOQ_LOG10 <- 2.0

# Analyze observations through day 196
MAX_DAY <- 196

# How to handle values below the detection limit:
# "half_LLOQ" = replace with half the detection limit
# "exclude"   = remove these observations
CENSOR_METHOD <- "half_LLOQ"

# Keep usable observations and remove observations after day 196
dat <- df %>%
  filter(
    !is.na(Patid),
    !is.na(Day),
    !is.na(log10.RNA.),
    Day <= MAX_DAY
  ) %>%
  arrange(Patid, Day)

# Get each patient's first available CD4 measurement
baseline_cd4 <- dat %>%
  group_by(Patid) %>%
  summarise(
    baseline_cd4 = if (all(is.na(CD4))) {
      NA_real_
    } else {
      first(CD4[!is.na(CD4)])
    },
    .groups = "drop"
  )

# Add baseline CD4 and center it around the population median
dat <- dat %>%
  left_join(baseline_cd4, by = "Patid") %>%
  mutate(
    CD4B = baseline_cd4 / 100
  )

median_cd4 <- median(dat$CD4B, na.rm = TRUE)

dat <- dat %>%
  mutate(
    CD4B = CD4B - median_cd4,
    cens = log10.RNA. <= LLOQ_LOG10
  )

# Report how much data is below the detection limit
cat("Number of patients:", n_distinct(dat$Patid), "\n")
cat("Number of observations:", nrow(dat), "\n")
cat(
  "Percentage at or below LLOQ:",
  round(mean(dat$cens) * 100, 1), "%\n"
)

# Handle below-limit values
if (CENSOR_METHOD == "exclude") {
  dat <- dat %>% filter(!cens)
} else if (CENSOR_METHOD == "half_LLOQ") {
  dat <- dat %>%
    mutate(
      log10.RNA. = ifelse(
        cens,
        log10((10^LLOQ_LOG10) / 2),
        log10.RNA.
      )
    )
} else {
  stop("CENSOR_METHOD must be 'half_LLOQ' or 'exclude'.")
}

# Keep the columns used by the model
dat <- dat %>%
  filter(!is.na(CD4B)) %>%
  transmute(
    patid = factor(Patid),
    day = Day,
    DV = log10.RNA.,     #DV = Dependent variable
    CD4B = CD4B,
    cens = cens
  ) %>%
  as.data.frame()

# Check the prepared dataset
head(dat)
summary(dat)

p_raw <- ggplot(
  dat,
  aes(x = day, y = DV, group = patid)
) +
  geom_line(alpha = 0.4) +
  labs(
    title = "ACTG 315: Observed viral-load profiles",
    x = "Days since treatment start",
    y = "HIV-1 RNA (log10 copies/mL)"
  ) +
  theme_classic()

print(p_raw)


##############  Fit model ############################# 

# Define the predicted log10 viral load
# p1 is the log of the starting amount of phase-1 virus.
# exp(ll1) is the decay rate. 
# Estimating it on the log scale (ll1) forces the rate to be positive.
# p2 is the log of the starting amount of phase-2 virus.
# l2 is the decay rate, not log-transformed here, so it could go negative.
model_fun <- function(day, p1, ll1, p2, l2) {
  log10(
    exp(p1 - exp(ll1) * day) +
      exp(p2 - l2 * day)
  )
}

# Starting values are initial guesses, not final estimates
start_values <- c(
  p1 = 11,
  ll1 = log(0.3),
  p2 = 6,
  l2_intercept = 0.01,
  l2_CD4B = 0
)

# Try different random-effect structures
random_options <- list(
  "p1 + ll1 + p2 + l2" =
    pdDiag(p1 + ll1 + p2 + l2 ~ 1),
  
  "p1 + ll1 + l2" =
    pdDiag(p1 + ll1 + l2 ~ 1),
  
  "p1 + ll1" =
    pdDiag(p1 + ll1 ~ 1)
)

fit <- NULL
random_effects_used <- NULL

for (option_name in names(random_options)) {
  
  cat("\nTrying random-effects model:", option_name, "\n")
  
  fit <- tryCatch(
    nlme(
      DV ~ log10(
        exp(p1 - exp(ll1) * day) +
          exp(p2 - l2 * day)
      ),
      data = dat,
      
      # Population-level parameters
      fixed = list(
        p1 ~ 1,
        ll1 ~ 1,
        p2 ~ 1,
        l2 ~ CD4B
      ),
      
      # Patient-level variation
      random = random_options[[option_name]],
      groups = ~ patid,
      
      # Initial guesses for the five fixed-effect coefficients
      start = unname(start_values),
      
      control = nlmeControl(
        maxIter = 200,
        pnlsMaxIter = 20,
        msMaxIter = 200,
        tolerance = 1e-5
      )
    ),
    error = function(e) {
      message("This model failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(fit)) {
    random_effects_used <- option_name
    break
  }
}

if (is.null(fit)) {
  stop(
    "None of the models converged. Check the data and starting values."
  )
}

cat("\nRandom-effects structure selected:",
    random_effects_used, "\n")

# Review the model fit
summary(fit)

# Extract population-level parameter estimates
fixed_effects <- fixef(fit)

cat("\nPopulation fixed-effect estimates:\n")
print(fixed_effects)

# Calculate approximate 95% confidence intervals
fixed_intervals <- tryCatch(
  intervals(fit, which = "fixed")$fixed,
  error = function(e) {
    message("Could not calculate confidence intervals: ",
            conditionMessage(e))
    NULL
  }
)

if (!is.null(fixed_intervals)) {
  print(fixed_intervals)
}

# Derived quantity: first-phase half-life
first_phase_rate <- exp(fixed_effects["ll1"])

first_phase_half_life <- log(2) / first_phase_rate

cat(
  "First-phase decay rate:",
  round(first_phase_rate, 4),
  "per day\n"
)

cat(
  "First-phase half-life:",
  round(first_phase_half_life, 2),
  "days\n"
)

# Individual predictions include patient-level random effects
diagnostics <- dat %>%
  mutate(
    IPRED = as.numeric(fitted(fit, level = 1)), #individual predictions (include random effects)
    PRED = as.numeric(fitted(fit, level = 0)), # population predictions (fixed effects only)
    RES = as.numeric(resid(fit, type = "pearson")) #standardized residuals
  )

# Observed versus individual predictions
p_obs_ipred <- ggplot(
  diagnostics,
  aes(x = IPRED, y = DV)
) +
  geom_point(alpha = 0.6) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Observed vs Individual Predictions",
    x = "Individual predicted viral load",
    y = "Observed viral load"
  ) +
  theme_classic()

print(p_obs_ipred)

# Residuals versus time
p_res_time <- ggplot(
  diagnostics,
  aes(x = day, y = RES)
) +
  geom_point(alpha = 0.6) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Residuals vs Time",
    x = "Days since treatment start",
    y = "Pearson residual"
  ) +
  theme_classic()

print(p_res_time)


# Extra Plot 1: Residuals over time with a smooth trend line
p_res_smooth <- ggplot(
  diagnostics,
  aes(x = day, y = RES)
) +
  geom_point(alpha = 0.5) +
  
  # The horizontal line means zero error
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  # The curved line shows the overall trend
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    color = "blue"
  ) +
  
  labs(
    title = "Residuals vs Time with Smooth Trend",
    x = "Days since treatment start",
    y = "Pearson residual"
  ) +
  theme_classic()

print(p_res_smooth)
#Result --> Good sign
# The blue curve stays approximately horizontal around zero. 
# there is no strong systematic pattern in the residuals over time.

# Extra Plot 2: Check whether the size of errors changes over time

p_res_spread <- ggplot(
  diagnostics,
  aes(x = day, y = abs(RES))
) +
  geom_point(alpha = 0.5) +
  
  # Smooth line shows how error size changes over time
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    color = "red"
  ) +
  
  labs(
    title = "Absolute Residuals vs Time",
    x = "Days since treatment start",
    y = "Absolute Pearson residual"
  ) +
  theme_classic()

print(p_res_spread)

#result --> 
# Red line approximately flat: No evidence that residual magnitude changes systematically over time.
# Red line rises: Errors tend to become larger later in the study.
# Red line falls: Errors tend to become smaller over time.

p_individual <- ggplot(
  diagnostics,
  aes(x = day)
) +
  geom_point(
    aes(y = DV),
    alpha = 0.6,
    size = 1
  ) +
  geom_line(
    aes(y = IPRED),
    color = "firebrick"
  ) +
  geom_line(
    aes(y = PRED),
    color = "steelblue",
    linetype = "dashed"
  ) +
  facet_wrap(~ patid) +
  labs(
    title = "Observed and Predicted Viral-Load Profiles",
    x = "Days since treatment start",
    y = "HIV-1 RNA (log10 copies/mL)"
  ) +
  theme_classic()

print(p_individual)

summary(fit)

##############  VPC #############################

# ---- USER-EDITABLE SETTINGS ----
N_SIM <- 1000; N_BINS <- 8; SEED <- 1
# --------------------------------
set.seed(SEED)
# the model’s uncertainty is split into three parts: 
#1.fixed effects, 2. between-patient SDs, and 3.residual noise.
fe  <- fixef(fit)
sdr <- sqrt(diag(pdMatrix(fit$modelStruct$reStruct$patid))) * fit$sigma
names(sdr) <- colnames(ranef(fit))
sig <- fit$sigma
ids <- as.character(dat$patid); uid <- unique(ids)
print(fe); print(sdr)

simulate <- function() {
  eta <- sapply(sdr, function(s) rnorm(length(uid), 0, s))
  rownames(eta) <- uid
  g  <- function(nm) eta[ids, nm]
  l2 <- fe[["l2.(Intercept)"]] + fe[["l2.CD4B"]] * dat$CD4B + g("l2.(Intercept)")
  y  <- log10(exp(fe[["p1"]] + g("p1") - exp(fe[["ll1"]] + g("ll1")) * dat$day) +
                exp(fe[["p2"]] + g("p2") - l2 * dat$day)) +
    rnorm(nrow(dat), 0, sig)
  pmax(y, log10(10^LLOQ_LOG10 / 2))
}

y <- simulate(); stopifnot(!anyNA(y)); summary(y)   # check before the loop

dat$bin <- cut(dat$day,
               unique(quantile(dat$day, seq(0, 1, length.out = N_BINS + 1))),
               include.lowest = TRUE)
pct <- function(x) quantile(x, c(.05, .5, .95))

sim <- do.call(rbind, lapply(seq_len(N_SIM), function(i)
  dplyr::reframe(data.frame(bin = dat$bin, y = simulate()),
                 q = c("P05", "P50", "P95"), v = pct(y), .by = bin)))

vpc <- sim |> dplyr::group_by(bin, q) |>
  dplyr::summarise(lo = quantile(v, .05), hi = quantile(v, .95), .groups = "drop")
obs <- dplyr::reframe(dat, q = c("P05", "P50", "P95"), v = pct(DV), .by = bin)
mid <- dat |> dplyr::group_by(bin) |> dplyr::summarise(day = median(day))
vpc <- merge(merge(vpc, obs, by = c("bin", "q")), mid, by = "bin")

ggplot(vpc, aes(day)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = q), alpha = .3) +
  geom_line(aes(y = v, colour = q), linewidth = .8) +
  labs(x = "Days since treatment start", y = "HIV-1 RNA (log10 copies/mL)",
       title = "VPC: observed percentiles vs 90% simulation intervals",
       fill = "Percentile", colour = "Percentile") +
  theme_classic()

intervals(fit, which = "fixed")           # proper CIs
summary(dat$CD4B)                         # CD4 scale, to fix the units
mean(dat$DV <= log10(50))                 # share of values at the floor

cf  <- coef(fit)                          # individual parameters
cd4 <- tapply(dat$CD4B, dat$patid, unique)
l2i <- cf[, "l2.(Intercept)"] + cf[, "l2.CD4B"] * cd4[rownames(cf)]
summary(l2i); sum(l2i < 0)                # second-phase rates; patients with a rising second phase
