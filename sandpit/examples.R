# Quick test

library(wtdr)
library(haven)
library(data.table)
library(tidyverse)
library(rlang)

# load data
# Dataset contains one observation per person, first dispensing for 2014
df <- haven::read_dta(system.file("extdata", "wtddat_dates.dta", package="wtdr"))

df <- cbind(df, sex = sample(c("M","F"), dim(df)[1], replace = T))

# fit waiting time distribution
fit1 <- wtdttt(data = df,
               rx1time ~ dlnorm(logitp, mu, lnsigma),
               id = "pid",
               start = as.Date('2014-01-01'), end = as.Date('2014-12-31'), subset = sex=='F'
               )

# N.B. We get an error if we set reverse = T because this dataset is created having the first dispensing for 2014

summary(fit1)

plot(fit1)

predict(fit1)
predict(fit1, type = "prob")

dfn <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

predict(fit1, type = "prob", prediction.data = dfn, distrx = "last_rxtime")


# repeat with Weibull
# can omit the "id" parameter as the data are pre-processed

fit2 <- wtdttt(data = df,
               rx1time ~ dweib(logitp, lnalpha, lnbeta),
                start = as.Date('2014-01-01'), end = as.Date('2014-12-31'),
)


summary(fit2)
plot(fit2)

predict(fit2, type="prob")

# ... and with Exponential
fit3 <- wtdttt(data = df,
               rx1time ~ dexp(logitp, lnbeta),
               start = as.Date('2014-01-01'), end = as.Date('2014-12-31'),
)


summary(fit3)

plot(fit3)

# Compare results with those obtained from wtdttt_ex.do in Stata ----

## 1) Ordinary WTD analysis with event dates (discrete time) ----
df <- haven::read_dta(system.file("extdata", "wtddat_dates.dta", package="wtdr"))

# fit waiting time distribution
fit1 <- wtdttt(data = df,
               rx1time ~ dlnorm(logitp, mu, lnsigma),
               id = "pid",
               start = as.Date('2014-01-01'), end = as.Date('2014-12-31')
)

summary(fit1)

## 2) Reverse WTD analysis with covariates (continuous time)----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

### fit waiting time distribution (exp)
fit3 <- wtdttt(data = df,
               last_rxtime ~ dexp(logitp, lnbeta),
               start = 0, end = 1, reverse = T, subset = packsize==200
)

summary(fit3)

### fit waiting time distribution (log-normal) ----
fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               # id = "pid",
               start = 0, end = 1, reverse = T
)

summary(fit1)

### fit waiting time distribution (weibull) ----
fit2 <- wtdttt(data = df,
               last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
               start = 0, end = 1, reverse = T
)

summary(fit2)

predict(fit2, type = "dur")

### Use covariate in all three parameter equations ----

# make packsize a factor
df <- df %>%
       mutate(packsize = as.factor(packsize))

# parameters = list(logitp ~ packcat, mu ~ packcat, lnsigma ~ 1)
fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ packsize)
)

summary(fit1)

### Since covariate appears to have little influence on the lnsigma parameter, we estimate a model where number of pills only affect median parameter (mu) and the prevalence (logitp) ----
fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ 1)
)

summary(fit1)

predict(fit1)

# weibull
fit1 <- wtdttt(data = df,
               last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, lnalpha ~ packsize, lnbeta ~ 1)
)

summary(fit1)

predict(fit1)

# exponential
fit1 <- wtdttt(data = df,
               last_rxtime ~ dexp(logitp, lnbeta),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, lnbeta ~ packsize)
)

summary(fit1)

predict(fit1)

### Try with more than one covariate ----
df <- df %>% mutate(sex = sample(c("F","M"), dim(df)[1], replace = T))

df <- df %>%
  mutate(sex = as.factor(sex))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize+sex, mu ~ packsize+sex, lnsigma ~ 1)
)

summary(fit1)

predict(fit1)

## 3) Prediction of treatment probability: A small example showing how treatment probability can be predicted based on the distance between index dates and date of last prescription redemption, while taking covariates (here: pack size) into account. The last fitted WTD is used for this prediction. ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ 1)
)

df <- haven::read_dta(system.file("extdata", "lastRx_index.dta", package="wtdr"))

df <- df %>%
       arrange(packsize, distlast) %>%
       mutate(packsize = as.factor(packsize))

prob_pred <- predict(fit1, type = "prob", prediction.data = df, distrx = "distlast")

# df <- df %>% arrange(packsize, distlast)
# plot(df$distlast, prob_pred, type = "l", ylim = c(0,1))

tmp <- cbind(df, prob_pred)

ggplot(data = tmp, aes(x=distlast, y=prob_pred, group = packsize)) +
  geom_line()

### try if the model without covariates matches stata results ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               # id = "pid",
               start = 0, end = 1, reverse = T
)

df <- haven::read_dta(system.file("extdata", "lastRx_index.dta", package="wtdr"))

df <- df %>% arrange(packsize, distlast)
prob_pred <- predict(fit1, type = "prob", prediction.data = df, distrx = "distlast")

ggplot(data = df, aes(x=distlast, y=prob_pred, group = packsize)) +
  geom_line() +
  geom_point()


## 4) Prediction of prescription duration (based on an estimated WTD). ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>%
  mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ 1)
)

summary(fit1)

predict(fit1, quantile = 0.9)

### with weibull distribution ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>%
  mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, lnalpha ~ 1, lnbeta ~ packsize)
)

summary(fit1)

predict(fit1, quantile = 0.9)

### with exponential distribution ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>%
  mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dexp(logitp, lnbeta),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, lnbeta ~ packsize)
)

summary(fit1)

predict(fit1, quantile = 0.9)


### with a new dataset for prediction ----
df <- haven::read_dta(system.file("extdata", "lastRx_index.dta", package="wtdr"))

df <- df %>%
  mutate(packsize = as.factor(packsize))

predict(fit1, quantile = 0.9, prediction.data = df)

### with more than one covariate (after having implemented model.frame in pred_dur_prob to extract var names) ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(packsize),
                    sex = as.factor(sample(c("F","M"), dim(df)[1], replace = T)))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize+sex, mu ~ packsize+sex, lnsigma ~ 1)
)

predict(fit1, quantile = 0.9)

### with more than one covariate and a new dataset (after having implemented model.frame in pred_dur_prob to extract var names) ----
df <- haven::read_dta(system.file("extdata", "lastRx_index.dta", package="wtdr"))

df <- df %>%
       mutate(packsize = as.factor(packsize),
              sex = as.factor(sample(c("F","M"), dim(df)[1], replace = T)))

predict(fit1, quantile = 0.9, prediction.data = df)

### with different covariates in the estimation and prediction dataset (checking the error message) ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(packsize),
                    sex = as.factor(sample(c("F","M"), dim(df)[1], replace = T)))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize+sex, mu ~ packsize+sex, lnsigma ~ 1)
)

df <- haven::read_dta(system.file("extdata", "lastRx_index.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(packsize),
                    sex = as.factor(sample(c("F","M"), dim(df)[1], replace = T)))

predict(fit1, quantile = 0.9, prediction.data = df)

### with a larger prediction dataset (multiple rows per subject) ----
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ 1)
)

df <- haven::read_dta(system.file("extdata", "ranwtddat_discdates.dta", package="wtdr"))

df <- df %>% mutate(packsize = as.factor(sample(c(100,200), dim(df)[1], replace = T)))

predict(fit1, quantile = 0.9, prediction.data = df)


## try if the model without covariates matches stata results
df <- haven::read_dta(system.file("extdata", "wtddat_covar.dta", package="wtdr"))

# make packsize a factor
df <- df %>%
  mutate(packsize = as.factor(packsize))

fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               # id = "pid",
               start = 0, end = 1, reverse = T
)

predict(fit1, quantile = 0.9)

## 5) Prediction of mean duration ----
fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               start = 0, end = 1, reverse = T, parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ 1)
)

# try if the model without covariates matches stata results
fit1 <- wtdttt(data = df,
               last_rxtime ~ dlnorm(logitp, mu, lnsigma),
               # id = "pid",
               start = 0, end = 1, reverse = T
)

predict(fit1, iadmean = T)


# Compare results with those obtained from ranwtdttt_ex.do in Stata ----

## 4) Random index date
df <- haven::read_dta(system.file("extdata", "ranwtddat_discdates.dta", package="wtdr"))

fit_r <- ranwtdttt(data = df,
                   rxdate ~ dlnorm(logitp, mu, lnsigma),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = T
)

summary(fit_r)


# Try on Jesper fake data (period: 02/01/95 - 31/08/22; drugs: clopidogrel, human insulin, metformin) ----

# load data
df <- haven::read_dta(system.file("extdata", "drugpakud.dta", package="wtdr"))

# filter clopidogrel dispensations
df <- df %>%
       filter(atc=="B01AC04") %>%
       rename(pid = id)

# fit waiting time distribution
fit4 <- wtdttt(data = df,
               rxdate ~ dlnorm(logitp, mu, lnsigma),
               id = "pid",
               start = as.Date('2015-01-01'), end = as.Date('2015-12-31')
)


summary(fit4)

plot(fit4)

predict(fit4)

# Extract a small subset and test with covariates ----
df <- haven::read_dta(system.file("extdata", "drugpakud.dta", package="wtdr"))
set.seed(1)
ids <- sample(unique(df$id), size=1000)
df2 <- df[df$id %in% ids,]

set.seed(1)
fit1 <- ranwtdttt(data = df2,
  rxdate ~ dlnorm(logitp, mu, lnsigma),
  id = "id",
  start = as.Date('2014-01-01'), end = as.Date('2014-12-31'),
  subset = atc=='B01AC04'&!is.na(ddd), reverse=T, nsamp=5,
  parameters=list(mu~log(ddd))
)
summary(fit1)
