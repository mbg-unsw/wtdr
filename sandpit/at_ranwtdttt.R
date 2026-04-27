# at_ranwtdttt - R functions and documentation

#' @include wtd-class.R wtdttt.R
NULL

#' Fit Waiting Time Distribution with antithetic random index times
#'
#' `at_ranwtdttt()` estimates maximum likelihood estimates for parametric
#' Waiting Time Distribution (WTD) based on observed prescription redemptions
#' with adjustment for covariates, using antithetic sampling of index times for
#' each individual. The sampling window \[start, end\] is divided into
#' \code{nint} equal sub-intervals. Within each sub-interval \eqn{(a, b)}, a
#' single uniform variate \eqn{U \sim \mathrm{Uniform}(0, 1)} is drawn per
#' individual, yielding two negatively correlated antithetic index dates:
#' \eqn{a + U(b-a)} and \eqn{a + (1-U)(b-a)}.  This gives \eqn{2 \times
#' \texttt{nint}} index dates per individual in total.  The negative
#' correlation between paired variates reduces Monte Carlo variance in the
#' likelihood relative to independent sampling with the same total number of
#' index dates.
#'
#' @param form an object of class "formula" (or one that can be coerced to that
#'   class): a symbolic description of the model to be fitted. The details of
#'   the model specification are given under 'Details'
#' @param parameters model formulae for distribution parameters
#' @param data an optional data frame, list or environment (or object coercible
#'   by as.data.frame to a data frame) containing the variables in the model.
#'   If not found in data, the variables are taken from environment(formula),
#'   typically the environment from which at_ranwtdttt is called.
#' @param id the name of the variable that identifies distinct individuals
#' @param start start of the sampling window within which random index date(s)
#'   are sampled
#' @param end end of the sampling window within which random index date(s) are
#'   sampled
#' @param reverse logical; Fit the reverse waiting time distribution.
#' @param nint number of equal sub-intervals into which \[start, end\] is
#'   divided. Each sub-interval contributes two antithetic index dates per
#'   individual (one based on \eqn{U}, one on \eqn{1-U}), for a total of
#'   \eqn{2 \times \texttt{nint}} index dates per individual.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param robust logical; compute a robust estimate of variance.
#' @param na.action a function which indicates what should happen when the data
#'   contain NAs. The default is set by the na.action setting of options, and
#'   is na.fail if that is unset. The 'factory-fresh' default is na.omit.
#'   Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param init starting values for the parameters.
#' @param control a list of parameters for controlling the fitting process.
#' @param ... further arguments passed to other methods.
#'
#' @return at_ranwtdttt returns an object of class "wtd" inheriting from "mle".
#' @importFrom data.table data.table setDT := .N .SD as.data.table setnames
#' @importFrom stats runif
#' @export
at_ranwtdttt <- function(data, form, parameters = NULL, start = NA, end = NA,
                         reverse = FALSE, id = NA, nint = 10, subset = NULL,
                         robust = TRUE, na.action = na.omit, init = NULL,
                         control = NULL, ...) {

  if (is.null(data) || (nrow(data) < 1)) {
    stop("data must be non-empty")
  }

  if (!inherits(form, "formula") || attr(terms(form), "response") == 0) {
    stop("obstime variable must be specified in model formula")
  }

  data <- as.data.table(data)

  if (!is.null(substitute(subset))) {

    rows     <- enquo(subset)
    rows_val <- eval_tidy(rows, data)
    data     <- data[rows_val, ]

    if (nrow(data) < 1) {
      stop("data must be non-empty")
    }

  }

  obs.name   <- all.vars(form)[1]
  covar.names <- unique(unlist(lapply(parameters, function(x) all.vars(x)[-1])))

  if (!(obs.name %in% names(data))) {
    stop(paste0("'", obs.name, "'", "is not in data"))
  }

  data <- na.action(data, cols = c(obs.name, covar.names))

  ## Validate date types

  if (!is(data[[obs.name]], "Date") || !is(start, "Date") || !is(end, "Date"))
    stop(paste0("variables start, end and '", obs.name,
                "' must be all of class Date"))

  delta <- as.numeric(end - start)  # total window length in days

  if (is.null(id) || length(id) != 1 || is.na(id)) {
    stop("The id variable must be provided")
  }

  if (!(id %in% names(data))) {
    stop(paste0("'", id, "'", "is not in data"))
  }

  if (!is.numeric(nint) || length(nint) != 1 || nint < 1 ||
      nint != floor(nint)) {
    stop("nint must be a single positive integer")
  }
  nint <- as.integer(nint)

  ## Define 'id' as key so it can be used to assign random offsets

  kc <- c(id, obs.name)
  setkeyv(data, kc)

  .id <- c(id)

  off <- data.table(id = unique(data[[id]]), key = "id")

  ## Antithetic sampler ---------------------------------------------------
  ##
  ## For sub-interval k (k = 1, ..., nint):
  ##   a_num = (k-1) * delta / nint   (numeric offset from start, in days)
  ##   b_num =  k    * delta / nint
  ##   int_len = delta / nint          (sub-interval length in days)
  ##
  ## Draw U ~ Uniform(0,1) independently per individual.
  ## Antithetic index dates (as numeric offsets from start, floored to integer
  ## days):
  ##   offset_1 = floor(a_num + U       * int_len)   -> indda = start + offset_1
  ##   offset_2 = floor(a_num + (1 - U) * int_len)   -> indda = start + offset_2
  ##
  ## Each call to f_at(k) returns the stacked results for both antithetic
  ## variates, contributing 2 index dates per individual per interval.

  int_len <- delta / nint   # may be non-integer; computed once outside f_at

  ## Inner helper: run the data.table look-up for one vector of index dates
  ## (stored in off$indda) and return the shifted dispense times.

  run_one <- function() {

    if (!reverse) {

      data[off, indda := i.indda][
        data[[obs.name]] >= indda & data[[obs.name]] <= (indda + delta),
        .SD[1L], by = .id][,
                           rxshift := .obs - (indda - start), env = list(.obs = obs.name)]

    } else {

      data[off, indda := i.indda][
        data[[obs.name]] <= indda & data[[obs.name]] >= (indda - delta),
        .SD[.N], by = .id][,
                           rxshift := .obs + (end - indda), env = list(.obs = obs.name)]

    }

  }

  ## Outer function: iterate over sub-intervals, drawing antithetic pairs

  f_at <- function(k) {

    a_num <- (k - 1L) * int_len          # left edge of sub-interval (days)

    U <- runif(n = nrow(off))            # one U per individual

    ## First antithetic variate (U)
    off[, indda := start + floor(a_num + U * int_len)]
    res1 <- run_one()

    ## Second antithetic variate (1 - U)
    off[, indda := start + floor(a_num + (1 - U) * int_len)]
    res2 <- run_one()

    rbind(res1, res2)

  }

  tmp <- do.call(rbind, lapply(seq_len(nint), f_at))

  ## Parse distribution family from formula

  disttmp <- attr(terms(form, specials = c("dlnorm", "dweib", "dexp")),
                  "specials")

  dist <- if (isTRUE(disttmp$dlnorm == 2)) "lnorm"
  else if (isTRUE(disttmp$dweib  == 2)) "weib"
  else if (isTRUE(disttmp$dexp   == 2)) "exp"
  else stop("model must use one of dlnorm, dweib or dexp")

  newform <- switch(dist,
                    lnorm = rxshift ~ dlnorm(logitp, mu, lnsigma),
                    weib  = rxshift ~ dweib(logitp, lnalpha, lnbeta),
                    exp   = rxshift ~ dexp(logitp, lnbeta)
  )

  if (nrow(tmp) == 0)
    stop("All dates are out of the window defined by start and end")

  ## Fit the model

  out <- wtdttt(form = newform, parameters = parameters,
                start = start, end = end, reverse = reverse, id = id,
                preprocess = FALSE, init = init, data = tmp)

  if (robust) {
    vcov_s    <- sand_vcov(out)
    out@vcov  <- vcov_s
  }

  return(out)

}
