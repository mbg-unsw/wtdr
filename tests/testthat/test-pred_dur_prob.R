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
