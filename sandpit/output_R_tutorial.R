# load packages
if (!require("haven", character.only = TRUE)) install.packages("haven", dependencies = TRUE)
library(haven)
if (!require("tidyverse", character.only = TRUE)) install.packages("tidyverse", dependencies = TRUE)
library(tidyverse)
library(wtdr)
library(bbmle)
library(data.table)
library(tidyverse)
library(rlang)
library(boot)
library(Deriv)


# load data (only warfarin dispensation filtered)
df <- haven::read_dta("inst/extdata/wtd_tutorial_data.dta")


# make agecat and sex factor
df <- df %>%
  mutate(agecat = as.factor(agecat),
         sex = as.factor(sex))


# 1) fit model without covariates ----
fit1 <- wtdttt(data = df,
               rxdate ~ dlnorm(logitp, mu, lnsigma),
               id = "id",
               start = as.Date('2015-01-01'),
               end = as.Date('2015-12-31'))

# summary of the model
summary(fit1)

# estimation of the duration and its uncertainty
predict(fit1, quantile = 0.9)


# 2) fit model with one covariate on all parameters ----
fit2 <- wtdttt(data = df,
               rxdate ~ dlnorm(logitp, mu, lnsigma),
               id = "id",
               start = as.Date('2015-01-01'),
               end = as.Date('2015-12-31'),
               reverse = T,
               parameters = list(logitp ~ sex, mu ~ sex, lnsigma ~ sex))

# summary of the model
summary(fit2)

# estimation of the duration and its uncertainty
predict(fit2, quantile = 0.9)


# 3) fit model with one covariate on 2 parameters out of 3 ----
fit3 <- wtdttt(data = df,
               rxdate ~ dlnorm(logitp, mu, lnsigma),
               id = "id",
               start = as.Date('2015-01-01'),
               end = as.Date('2015-12-31'),
               reverse = T,
               parameters = list(logitp ~ sex, mu ~ sex, lnsigma ~ 1))

# summary of the model
summary(fit3)

# estimation of the duration and its uncertainty
unique(predict(fit3, quantile = 0.9))


# 4) fit model with random index dates ----
set.seed(489)

fit4 <- ranwtdttt(data = df,
                  rxdate ~ dlnorm(logitp, mu, lnsigma),
                  id = "id",
                  nsamp = 5,
                  start = as.Date('2015-01-01'),
                  end = as.Date('2015-12-31'),
                  reverse = T
)

# summary of the model
summary(fit4)

# estimation of the duration and its uncertainty
predict(fit4, quantile = 0.9)

