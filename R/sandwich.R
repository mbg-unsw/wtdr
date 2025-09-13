# calculate sandwich vcov

#' Calculate a robust variance-covariance matrix using the sandwich estimator
#'
#' @param fit an object of class "wtd" returned by `ranwtdttt()`
#'
#' @return sand_vcov returns a matrix
#' @importFrom stats model.matrix numericDeriv
sand_vcov <- function(fit) {
  # parse model formula components
  parm_form <- unlist(strsplit(gsub(" ", "", unlist(strsplit(fit@formula, ":", fixed=T))[2]), ",", fixed=T))

#####
# based on code from bbmle::calc_mle2_function
# FIXME: redo without special case code for each density type
# FIXME: redo to handle non-standard parameter naming

  parnames <- grep("delta", names(fit@call$start), value=T, fixed=T, invert=T)

  if (typeof(parm_form)=="character" && length(parm_form)>0) {
    ## linear model specified for some parameters
    vars <- sapply(strsplit(parm_form, "~", fixed=T),"[",1)
    models <-  paste0("~", sapply(strsplit(parm_form, "~", fixed=T),"[",2))
  } else {parm_form <- c(); vars <- c(); models <- c()}

  # FIXME for any missing components, substitute "~1"

    for (i in seq(along=parnames)) {
      if (!(parnames[i] %in% vars)) {
        vars <- c(vars, parnames[i])
        models <- c(models, "~1")
        parm_form <- c(parm_form, paste0(parnames[i], "~1"))
      }
    }

    vpos <- list()
    for (i in seq(along=parm_form)) {
      vname <- vars[i]      ## name of variable
      vpos[[vname]] <- i
    }


    # browser()

    if(fit@dist=="lnorm") {

      mm1 <- model.matrix(formula(models[vpos[["logitp"]]]), data=as.data.frame(fit@data), drop=F)
      mm2 <- model.matrix(formula(models[vpos[["mu"]]]), data=as.data.frame(fit@data), drop=F)
      mm3 <- model.matrix(formula(models[vpos[["lnsigma"]]]), data=as.data.frame(fit@data), drop=F)

      myenv <- new.env()
      myenv$logitp <- mm1 %*% fit@coef[1:ncol(mm1)]
      myenv$mu <- mm2 %*% fit@coef[(ncol(mm1)+1):(ncol(mm1)+ncol(mm2))]
      myenv$lnsigma <- mm3%*% fit@coef[(ncol(mm1)+ncol(mm2)+1):(ncol(mm1)+ncol(mm2)+ncol(mm3))]
      myenv$x <- as.numeric(getElement(fit@data, fit@depvar))

      rm(mm1, mm2, mm3)

      fun_logitp <- function(logitp) {

        dlnorm(myenv$x, logitp, myenv$mu, myenv$lnsigma, delta = fit@delta, log=T)

      }

      fun_mu <- function(mu) {

        dlnorm(myenv$x, myenv$logitp, mu, myenv$lnsigma, delta = fit@delta, log=T)

      }

      fun_lnsigma <- function(lnsigma) {

        dlnorm(myenv$x, myenv$logitp, myenv$mu, lnsigma, delta = fit@delta, log=T)

      }


      sc1 <- numDeriv::grad(fun_logitp, myenv$logitp)
      sc2 <- numDeriv::grad(fun_mu, myenv$mu)
      sc3 <- numDeriv::grad(fun_lnsigma, myenv$lnsigma)

      mm1 <- model.matrix(formula(models[vpos[["logitp"]]]), data=as.data.frame(fit@data), drop=F)
      mm2 <- model.matrix(formula(models[vpos[["mu"]]]), data=as.data.frame(fit@data), drop=F)
      mm3 <- model.matrix(formula(models[vpos[["lnsigma"]]]), data=as.data.frame(fit@data), drop=F)

      score_mat <- cbind(mm1 * sc1, mm2 * sc2, mm3 * sc3)

    } else if(fit@dist=="weib") {

      mm1 <- model.matrix(formula(models[vpos[["logitp"]]]), data=as.data.frame(fit@data), drop=F)
      mm2 <- model.matrix(formula(models[vpos[["lnalpha"]]]), data=as.data.frame(fit@data), drop=F)
      mm3 <- model.matrix(formula(models[vpos[["lnbeta"]]]), data=as.data.frame(fit@data), drop=F)

      myenv <- new.env()
      myenv$logitp <- mm1 %*% fit@coef[1:ncol(mm1)]
      myenv$lnalpha <- mm2 %*% fit@coef[(ncol(mm1)+1):(ncol(mm1)+ncol(mm2))]
      myenv$lnbeta <- mm3%*% fit@coef[(ncol(mm1)+ncol(mm2)+1):(ncol(mm1)+ncol(mm2)+ncol(mm3))]
      myenv$x <- as.numeric(getElement(fit@data, fit@depvar))

      sc1 <- diag(attr(numericDeriv(
        quote(dweib(x, logitp, lnalpha, lnbeta, delta = fit@delta, log=T)), c("logitp"), myenv), "gradient"))
      sc2 <- diag(attr(numericDeriv(
        quote(dweib(x, logitp, lnalpha, lnbeta, delta = fit@delta, log=T)), c("lnalpha"), myenv), "gradient"))
      sc3 <- diag(attr(numericDeriv(
        quote(dweib(x, logitp, lnalpha, lnbeta, delta = fit@delta, log=T)), c("lnbeta"), myenv), "gradient"))

      score_mat <- cbind(mm1 * sc1, mm2 * sc2, mm3 * sc3)

    } else if(fit@dist=="exp") {

      mm1 <- model.matrix(formula(models[vpos[["logitp"]]]), data=as.data.frame(fit@data), drop=F)
      mm2 <- model.matrix(formula(models[vpos[["lnbeta"]]]), data=as.data.frame(fit@data), drop=F)

      myenv <- new.env()
      myenv$logitp <- mm1 %*% fit@coef[1:ncol(mm1)]
      myenv$lnbeta <- mm2 %*% fit@coef[(ncol(mm1)+1):(ncol(mm1)+ncol(mm2))]
      myenv$x <- as.numeric(getElement(fit@data, fit@depvar))

      sc1 <- diag(attr(numericDeriv(
        quote(dexp(x, logitp, lnbeta, delta = fit@delta, log=T)), c("logitp"), myenv), "gradient"))
      sc2 <- diag(attr(numericDeriv(
        quote(dexp(x, logitp, lnbeta, delta = fit@delta, log=T)), c("lnbeta"), myenv), "gradient"))

      score_mat <- cbind(mm1 * sc1, mm2 * sc2)

    }

  colnames(score_mat) <- names(fit@coef)

  bread <- fit@vcov

  applyTapplySum <- function(X,index) apply(X, 2, function(col) tapply(col, index, sum))

  psi <- applyTapplySum(as.matrix(score_mat), getElement(fit@data, fit@idvar)) # sum by cluster

  meat <- crossprod(as.matrix(psi))

  rownames(meat) <- colnames(meat) <- colnames(psi)

  nid <- NROW(unique(getElement(fit@data, fit@idvar)))

  (nid/(nid-1)) * (bread %*% meat %*% bread)
}
