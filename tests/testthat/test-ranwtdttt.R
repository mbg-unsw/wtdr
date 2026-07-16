# test ranwtdttt()
# * all tests from wtdttt()
# - plus to be done for wtdttt():
# - test errors relating to start and end (types, values etc)
# - naming of parameters, including naming of parameter formulae
# - formatted output of print(), summary() [move to wtd-class testing?]
# - subset
# - init

# * nsamp
# * robust

testthat::test_that("errors", {

  testthat::expect_error(
    ranwtdttt(NULL, form = rx1time ~ dlnorm(logitp, mu, lnsigma), start=0, end=0),
    "data must be non-empty"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=double()), form = rx1time ~ dlnorm(logitp, mu, lnsigma), start=0, end=1),
    "data must be non-empty"
  )

  testthat::expect_error(
    wtdttt(data.frame(rx1time=c(1)), form = ~ dlnorm(logitp, mu, lnsigma), start=0, end=1),
    "obstime variable must be specified"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(1)), form = foo ~ dlnorm(logitp, mu, lnsigma), start=0, end=1),
    "is not in data"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(1)), form = rx1time ~ dlnorm(logitp, mu, lnsigma), start=as.Date("2023-01-01"), end=as.Date("2024-01-01")),
    "variables start, end and 'rx1time' must be either all of class Date or all of class numeric"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=as.Date("2024-01-01")), form = rx1time ~ dlnorm(logitp, mu, lnsigma), start=0, end=1),
    "variables start, end and 'rx1time' must be either all of class Date or all of class numeric"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(1)), form = rx1time ~ dlnorm(logitp, mu, lnsigma), start=0, end=1, subset=rx1time!=1),
    "data must be non-empty"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=as.Date("2024-01-01")), form = rx1time ~ dlnorm(logitp, mu, lnsigma),
              start=as.Date("2024-01-01"), end=as.Date("2025-01-01")),
    "The id variable must be provided"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=as.Date("2024-01-01")), form = rx1time ~ dlnorm(logitp, mu, lnsigma),
              start=as.Date("2024-01-01"), end=as.Date("2025-01-01"), id=c(1,2)),
    "The id variable must be provided"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(as.Date("2020-01-01"))),
              form = rx1time ~ dlnorm(logitp, mu, lnsigma),
              start=as.Date("2020-01-01"), end=as.Date("2021-01-01"), id="foo"),
    "is not in data"
  )

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(as.Date("2020-07-01")), id=c(1)),
              form = rx1time ~ dpois(lambda),
              start=as.Date("2020-01-01"), end=as.Date("2021-01-01"), id="id"),
    "model must use one of dlnorm, dweib or dexp"
  )

})

testthat::test_that("preprocessing errors", {

  testthat::expect_error(
    ranwtdttt(data.frame(rx1time=c(as.Date("2020-01-01")), id=c(1)),
              form = rx1time ~ dlnorm(logitp, mu, lnsigma),
              start=as.Date("2017-01-01"), end=as.Date("2018-01-01"),
              id="id"),
    "All dates are out of the window"
  )

  rd <- readRDS(test_path("fixtures", "randat_disc.rds"))
  npid <- length(unique(rd$pid))

  set.seed(127)

  # XXXX repeat all tests a few times to allow random variation in sampling?

  x <- ranwtdttt(data = rd,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = FALSE,
                 robust = FALSE
  )

  testthat::expect_lte(max(x@data$rxshift), 365.5)
  testthat::expect_gte(min(x@data$rxshift), 0.5)
  testthat::expect_lte(length(x@data$pid), npid)
  testthat::expect_equal(length(x@data$pid), length(unique(x@data$pid)))

  # XXXX need to think about which obs are excluded so can test for them

  set.seed(127)
  x <- ranwtdttt(data = rd,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = FALSE
  )

  testthat::expect_lte(max(x@data$rxshift), 365.5)
  testthat::expect_gte(min(x@data$rxshift), 0.5)
  testthat::expect_lte(length(x@data$pid), npid)
  testthat::expect_equal(length(x@data$pid), length(unique(x@data$pid)))

  # now test when nsamp > 1
  # check different scenarios, each one both forward and reverse
  # XXXX

  set.seed(127)

  # XXXX repeat all tests a few times to allow random variation in sampling?

  x <- ranwtdttt(data = rd,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = FALSE,
                 robust = FALSE,
                 nsamp = 5
  )

  testthat::expect_lte(max(x@data$rxshift), 365.5)
  testthat::expect_gte(min(x@data$rxshift), 0.5)
  testthat::expect_lte(length(x@data$pid), 5*npid)
  testthat::expect_lte(length(x@data$pid), 5*length(unique(x@data$pid)))

  x <- ranwtdttt(data = rd,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = FALSE,
                 nsamp = 5
  )

  testthat::expect_lte(max(x@data$rxshift), 365.5)
  testthat::expect_gte(min(x@data$rxshift), 0.5)
  testthat::expect_lte(length(x@data$pid), 5*npid)
  testthat::expect_lte(length(x@data$pid), 5*length(unique(x@data$pid)))

  # Quick check that alphanumeric ids work OK

  rd$apid <- sprintf("A%05d", rd$pid)

  x <- ranwtdttt(data = rd,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "apid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = FALSE
  )

  testthat::expect_lte(max(x@data$rxshift), 365.5)
  testthat::expect_gte(min(x@data$rxshift), 0.5)
  testthat::expect_lte(length(x@data$pid), npid)
  testthat::expect_equal(length(x@data$pid), length(unique(x@data$pid)))
})

# simulate test data
# t <- c(runif(980, -1, 0) + rexp(980, rate=20), runif(30, max=1.5))
# dt_exp <- data.frame(t=t[t>-0.05])

# XXXX using a random sample of the example data instead, converted to rds

testthat::test_that("basics", {

  rd <- readRDS(test_path("fixtures", "randat_disc.rds"))

  set.seed(127)
  x <- ranwtdttt(data = rd,
                     rxdate ~ dlnorm(logitp, mu, lnsigma),
                     id = "pid",
                     start = as.Date('2014-01-01'),
                     end = as.Date('2014-12-31'),
                     reverse = TRUE,
                     robust = FALSE
  )

  testthat::expect_s4_class(x, "mle2")
  testthat::expect_s4_class(x, "wtd")

  testthat::expect_equal(as.vector(x@fullcoef["logitp"]), 1.1489, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["mu"]), 4.312, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnsigma"]), -1.660, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["delta"]), 365)

  # XXXX could check the log likelihood too

  testthat::expect_equal(as.vector(x@vcov), c(0.0133589293, 0.0001877469, 0.0024136993,
                                              0.0001877469, 0.0012496759, -0.004001219,
                                              0.0024136993, -0.0040012191, 0.0283231302), tolerance=0.001)

  set.seed(127)
  x <- ranwtdttt(data = rd,
                 rxdate ~ dweib(logitp, lnalpha, lnbeta),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = FALSE
  )

  testthat::expect_equal(as.vector(x@fullcoef["logitp"]), 1.1425, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnalpha"]), 1.616, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnbeta"]), -4.3886, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["delta"]), 365)

  # XXXX could check the log likelihood too

  testthat::expect_equal(as.vector(x@vcov), c( 0.0130308759, -0.0015652083, -0.0003395173,
                                              -0.0015652083,  0.0315019131, -0.0031927367,
                                              -0.0003395173, -0.0031927367,  0.0009624777), tolerance=0.001)

  set.seed(127)
  x <- ranwtdttt(data = rd,
                 rxdate ~ dexp(logitp, lnbeta),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = FALSE
  )

  testthat::expect_equal(as.vector(x@fullcoef["logitp"]), 1.4518, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnbeta"]), -3.8353, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["delta"]), 365)

  # XXXX could check the log likelihood too

  testthat::expect_equal(as.vector(x@vcov), c( 0.025053120, -0.003949551,
                                              -0.003949551,  0.003587246), tolerance=0.001)

  set.seed(127)
  x <- ranwtdttt(data = rd,
                 rxdate ~ dexp(logitp, lnbeta),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = FALSE,
                 robust = FALSE
  )

  testthat::expect_equal(as.vector(x@fullcoef["logitp"]), 1.2546, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnbeta"]), -3.7065, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["delta"]), 365)

  # XXXX could check the log likelihood too

  testthat::expect_equal(as.vector(x@vcov), c( 0.019456371, -0.003097151,
                                              -0.003097151,  0.003572717), tolerance=0.001)


  # now with robust=T

  set.seed(127)
  x <- ranwtdttt(data = rd,
                 rxdate ~ dexp(logitp, lnbeta),
                 id = "pid",
                 start = as.Date('2014-01-01'),
                 end = as.Date('2014-12-31'),
                 reverse = TRUE,
                 robust = TRUE,
                 nsamp = 5
  )

  testthat::expect_equal(as.vector(x@fullcoef["logitp"]), 1.3134, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["lnbeta"]), -3.8324, tolerance=0.001)
  testthat::expect_equal(as.vector(x@fullcoef["delta"]), 365)

  # XXXX could check the log likelihood too

  testthat::expect_equal(as.vector(x@vcov), c( 0.0115917136, -0.0004325173,
                                              -0.0004325173,  0.0006074480), tolerance=0.001)

})
