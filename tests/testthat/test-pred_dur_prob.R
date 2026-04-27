# Test predict.wtd
# * all distributions DONE
# * date or continuous data DONE
# * forward or reverse DONE
# * different values of delta DONE
# * type=c("dur", "prob") DONE
# * iadmean=c(FALSE, TRUE) DONE
# * quantile=0.8, 0.5 DONE
# * distrx=NULL
# * newdata=NULL
# * newdata and type="prob"
# * complete all cases for newdata with NA
# * se.fit=TRUE
# * linear predictors

testthat::test_that("errors", {

  # stop("Covariates used in the estimation are not in the prediction dataset (new data)")

})

v <- function(x) as.vector(x, mode="numeric")

testthat::test_that("predictions", {

  # durations, simple
  # continuous data

  dt_exp <- readRDS(test_path("fixtures", "dt_exp.rds"))

  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = t ~ dexp(logitp, lnbeta), start=0, end=1),
      "The id variable was not provided"
    ),
    "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1402, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.06040, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.08713, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-exp(x@coef[2])*x@data$t), tolerance=0.001)



  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = t ~ dweib(logitp, lnalpha, lnbeta), start=0, end=1),
      "The id variable was not provided"
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1168, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.03864, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.03087, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-((x@data$t*exp(x@coef[3]))^exp(x@coef[2]))), tolerance=0.001)



  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = t ~ dlnorm(logitp, mu, lnsigma), start=0, end=1),
      "The id variable was not provided"
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1348, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.06367, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.09468, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         pnorm(-(log(x@data$t)-x@coef[2])/exp(x@coef[3])), tolerance=0.001)






  # And for reverse...

  dt_exp$tr <- 1 - dt_exp$t

  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = tr ~ dexp(logitp, lnbeta), start=0, end=1, reverse=TRUE),
      "The id variable was not provided"
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1402, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.06040, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.08713, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-exp(x@coef[2])*x@data$t), tolerance=0.001)


  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = tr ~ dweib(logitp, lnalpha, lnbeta), start=0, end=1, reverse=TRUE),
      "The id variable was not provided"
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1168, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.03864, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.03087, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-((x@data$t*exp(x@coef[3]))^exp(x@coef[2]))), tolerance=0.001)


  testthat::expect_warning(
    testthat::expect_warning(
      x <- wtdttt(dt_exp, form = tr ~ dlnorm(logitp, mu, lnsigma), start=0, end=1, reverse=TRUE),
      "The id variable was not provided"
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(0.1348, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.06367, 62), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.09468, 62), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         pnorm(-(log(x@data$t)-x@coef[2])/exp(x@coef[3])), tolerance=0.001)


  # Date data

  rd <- readRDS(test_path("fixtures", "randat_disc.rds"))

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dexp(logitp, lnbeta),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = FALSE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(68.41, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(29.46, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(42.50, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-exp(x@coef[2])*x@data$rxdate), tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dweib(logitp, lnalpha, lnbeta),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = FALSE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(86.00, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(67.85, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(3.037, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-((x@data$rxdate*exp(x@coef[3]))^exp(x@coef[2]))), tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dlnorm(logitp, mu, lnsigma),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = FALSE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(85.06, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(68.81, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(71.03, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         pnorm(-(log(x@data$rxdate)-x@coef[2])/exp(x@coef[3])), tolerance=0.001)

  # reverse

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dexp(logitp, lnbeta),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = TRUE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(76.19, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(32.81, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(47.34, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-exp(x@coef[2])*x@data$rxdate), tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dweib(logitp, lnalpha, lnbeta),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = TRUE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(88.73, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(70.17, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(3.030, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         exp(-((x@data$rxdate*exp(x@coef[3]))^exp(x@coef[2]))), tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(data = rd,
                   rxdate ~ dlnorm(logitp, mu, lnsigma),
                   id = "pid",
                   start = as.Date('2014-01-01'),
                   end = as.Date('2014-12-31'),
                   reverse = TRUE
    ),
  "Some dates are out of the window"
  )

  testthat::expect_equal(v(predict(x)), rep(87.96, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(69.80, 642), tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(72.49, 642), tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         pnorm(-(log(x@data$rxdate)-x@coef[2])/exp(x@coef[3])), tolerance=0.001)

})

testthat::test_that("linear predictors", {

  dt_coef <- readRDS(test_path("fixtures", "dt_coef.rds"))
  dt_coef$packsize <- factor(dt_coef$packsize) # work around issue #33
  dt_coef$sex <- factor(dt_coef$sex)

  nn<-data.frame(packsize=factor(c("100", "200")), sex=factor(c("M", "F")))
  nn2<-data.frame(packsize=factor(c("100", "200", NA)), sex=factor(c("M", "F", NA)))
  nn3<-data.frame(packsize=factor(c("100", "200", NA)), sex=factor(c("M", "F", NA)),
                  tt=c(0.1,0.15,0.2))

  testthat::expect_warning(
    x <- wtdttt(dt_coef, form = last_rxtime ~ dexp(logitp, lnbeta),
                parameters = list(logitp ~ packsize),
                start=0, end=1, reverse=T),
    "The id variable was not provided"
  )

  testthat::expect_equal(v(predict(x)), rep(0.2192, 1000), tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)), rep(0.09441, 1000),
                         tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)), rep(0.1362, 1000),
                         tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         v(exp(-exp(x@coef[3])*x@data$last_rxtime)),
                         tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn)), rep(0.2192, 2), tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn2)), rep(0.2192, 3), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn2, na.action=na.omit)), rep(0.2192, 2), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt")), c(0.4799, 0.3324, 0.2303), tolerance=0.001)
# BUG
# testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt", na.action=na.omit)), ?????, tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(dt_coef, form = last_rxtime ~ dexp(logitp, lnbeta),
                parameters = list(lnbeta ~ packsize),
                start=0, end=1, reverse=T),
    "The id variable was not provided"
  )

  testthat::expect_equal(v(predict(x)),
                         ifelse(x@data$packsize==200, 0.2479, 0.1822),
                         tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)),
                         ifelse(x@data$packsize==200, 0.1067, 0.07847),
                         tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)),
                         ifelse(x@data$packsize==200, 0.1540, 0.1132),
                         tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         v(exp(-exp(x@coef[2]+ifelse(x@data$packsize==200, x@coef[3], 0))*x@data$last_rxtime)), tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn)), c(0.1822, 0.2479), tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn2)), c(0.1822, 0.2479, NA), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn2, na.action=na.omit)), c(0.1822, 0.2479), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt")), c(0.4134, 0.3776, NA), tolerance=0.001)
# BUG
# testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt", na.action=na.omit)), ?????, tolerance=0.001)

  testthat::expect_warning(
    x <- wtdttt(dt_coef, form = last_rxtime ~ dexp(logitp, lnbeta),
                parameters = list(logitp ~ packsize, lnbeta ~ packsize),
                start=0, end=1, reverse=T),
    "The id variable was not provided"
  )

  testthat::expect_equal(v(predict(x)),
                         ifelse(x@data$packsize==200, 0.2596, 0.1762),
                         tolerance=0.001)
  testthat::expect_equal(v(predict(x, quantile=0.5)),
                         ifelse(x@data$packsize==200, 0.1118, 0.0759),
                         tolerance=0.001)
  testthat::expect_equal(v(predict(x, iadmean=TRUE)),
                         ifelse(x@data$packsize==200, 0.1613, 0.1095),
                         tolerance=0.001)

  testthat::expect_equal(v(predict(x, type="prob")),
                         v(exp(-exp(x@coef[3]+ifelse(x@data$packsize==200, x@coef[4], 0))*x@data$last_rxtime)), tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn)), c(0.1762, 0.2596), tolerance=0.001)

  testthat::expect_equal(v(predict(x, newdata=nn2)), c(0.1762, 0.2596, NA), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn2, na.action=na.omit)), c(0.1762, 0.2596), tolerance=0.001)
  testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt")), c(0.4012, 0.3945, NA), tolerance=0.001)
# BUG
# testthat::expect_equal(v(predict(x, newdata=nn3, type="prob", distrx="tt", na.action=na.omit)), ?????, tolerance=0.001)

#   # repeat for dweib
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
#                 parameters = list(logitp ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
#                 parameters = list(lnalpha ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
#                 parameters = list(lnbeta ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dweib(logitp, lnalpha, lnbeta),
#                 parameters = list(lnalpha ~ packsize, lnbeta ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(logitp ~ packsize, lnsigma ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(logitp ~ packsize, mu ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(mu ~ packsize, lnsigma ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(mu ~ packsize, lnsigma ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(logitp ~ packsize, mu ~ packsize, lnsigma ~ packsize),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX
#
#   testthat::expect_warning(
#     x <- wtdttt(dt_coef, form = last_rxtime ~ dlnorm(logitp, mu, lnsigma),
#                 parameters = list(logitp ~ sex, mu ~ sex + log(packsize)),
#                 start=0, end=1, reverse=T),
#     "The id variable was not provided"
#   )
#
# # XXXX

})

