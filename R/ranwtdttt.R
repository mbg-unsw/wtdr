# ranwtdttt - R functions and documentation

#' @include wtd-class.R wtdttt.R
NULL

#' Fit Waiting Time Distribution with random index times
#'
#' `ranwtdttt()` estimates maximum likelihood estimates for parametric Waiting Time Distribution (WTD)
#' based on observed prescription redemptions with adjustment for covariates
#' using one or more random index times for each individual. It reports estimates
#' of prevalence fraction and specified percentile of inter-arrival density
#' together with regression coefficients.
#'
#' @param form an object of class "formula" (or one that can be coerced to that
#' class): a symbolic description of the model to be fitted. The details of the
#' model specification are given under 'Details'
#' @param parameters model formulae for distribution parameters
#' @param data an optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from environment(formula),
#' typically the environment from which wtdttt is called.
#' @param id the name of the variable that identifies distinct individuals
#' @param start start of the sampling window within which random index value(s) are sampled (date or integer)
#' @param end end of the sampling window within which random index value(s) are sampled (date or integer)
#' @param reverse logical; Fit the reverse waiting time distribution.
#' @param nsamp number of samples to take.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param robust logical; compute a robust estimate of variance.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is set by the na.action setting of options, and is
#' na.fail if that is unset. The 'factory-fresh' default is na.omit. Another
#' possible value is NULL, no action. Value na.exclude can be useful.
#' @param init starting values for the parameters.
#' @param control a list of parameters for controlling the fitting process.
#' @param ... further arguments passed to other methods.
#'
#' @return ranwtdttt returns an object of class "wtd" inheriting from "mle".
#' @importFrom data.table data.table setDT := .N .SD as.data.table setnames
#' @importFrom stats runif
#' @export
ranwtdttt <- function(data, form, parameters=NULL, start=NA, end=NA, reverse=F, id=NA,
                      nsamp=1, subset=NULL, robust=T, na.action=na.omit, init=NULL, control=NULL, ...) {

  if(is.null(data) || (nrow(data)<1)) {
    stop("data must be non-empty")
  }

  if(!inherits(form, "formula") || attr(terms(form), "response")==0) {
    stop("obstime variable must be specified in model formula")
  }

  data <- as.data.table(data)

  if(!is.null(substitute(subset))) {

    rows <- enquo(subset)
    rows_val <- eval_tidy(rows, data)
    data <- data[rows_val ,]

    if(nrow(data)<1) {
      stop("data must be non-empty")
    }

  }


  # 14/04/25
  obs.name <- all.vars(form)[1]
  covar.names <- unique(unlist(lapply(parameters, function(x) all.vars(x)[-1])))

  if(!(obs.name %in% names(data))) {
    stop(paste0("'", obs.name, "'", "is not in data"))
  }

  data <- na.action(data, cols = c(obs.name, covar.names))

  # check data

  if(!(is(data[[obs.name]], "Date") && is(start, "Date") && is(end, "Date"))
   && !(is(data[[obs.name]], "numeric") && is(start, "numeric") && is(end, "numeric")))
  stop(paste0("variables start, end and '", obs.name, "' must be either all of class Date or all of class numeric"))

  if(is(data[[obs.name]], "numeric"))
     warning("numeric data are assumed to be integer-valued")

  # creation of shifted dates

  delta <- as.numeric(end - start)

  if(is.null(id) || length(id)!=1 || is.na(id)) {
    stop("The id variable must be provided")
  }

  if(!(id %in% names(data))) {
    stop(paste0("'", id, "'", "is not in data"))
  }

  # define 'id' as key so it can be used to assign random offsets
  kc <- c(id, obs.name); setkeyv(data, kc) # XXXX consider using indexing

  obs.ind <- which(colnames(data)==obs.name)
  .id <- c(id)

  off <- data.table(id=unique(data[[id]]), key=c("id"))

  f <- function() {

    # all obs with the same id should get the same random offset
    off[, indda := start + floor(runif(n=nrow(off), max=delta+1))]

    if (!reverse) {

      data[off, indda := i.indda][data[[obs.name]] >= indda & data[[obs.name]] <= (indda + delta), .SD[1L], by=.id][, rxshift := .obs - (indda-start), env=list(.obs=obs.name)]

    } else {

      data[off, indda := i.indda][data[[obs.name]] <= indda & data[[obs.name]] >= (indda - delta), .SD[.N], by=.id][, rxshift := .obs + (end - indda), env=list(.obs=obs.name)]

    }

  }

  # bind multiple copies of dataframes
  tmp <- do.call(rbind, replicate(nsamp, f(), simplify = FALSE))

  # XXXX do formula rewrite like in wtdttt() ?
  disttmp <- attr(terms(form, specials=c("dlnorm", "dweib", "dexp")), "specials")

  dist <- if(isTRUE(disttmp$dlnorm==2)) "lnorm" # need isTRUE() as value can be NULL
          else if(isTRUE(disttmp$dweib==2)) "weib"
          else if(isTRUE(disttmp$dexp==2)) "exp"
          else stop("model must use one of dlnorm, dweib or dexp")

  if(dist == "lnorm") {

     newform <- rxshift ~ dlnorm(logitp, mu, lnsigma)

  } else if (dist == "weib") {

     newform <- rxshift ~ dweib(logitp, lnalpha, lnbeta)

  } else if (dist == "exp") {

     newform <- rxshift ~ dexp(logitp, lnbeta)

  }

  if(nrow(tmp)==0)
    stop("All dates are out of the window defined by start and end")


  # XXXX shouldn't apply subset both here and above
  # note, pass 'id' to wtdttt() but set preprocess=F
  # as when nsamp > 1 *all* copies of data should be used

  out <- wtdttt(form = newform, parameters = parameters,
                start = start, end = end, reverse = reverse, id = id,
                preprocess = F, init = init, data = tmp)

  if (!robust) {

    out <- out

  } else {

    vcov_s <- sand_vcov(out)
    out@vcov <- vcov_s

  }


  return(out)

  # parse 'form' to determine the distribution in use and test if it
  # is a supported one, otherwise error

  # parse 'parameters' to test if they match 'form', otherwise error

  # test if start, end are dates, error if a mix of types
  # test if outcome variable is a date
  # if start, end are dates and outcome is not, error
  # test if nsamp is an integer >=1

  # repeat nsamp times:
  #   generate a random start date for each id sample window
  #   extract dispense times within sample window for each id
  #   add to temporary analysis data set

  # call wtdttt using temporary data set

}
