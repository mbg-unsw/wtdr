# =============================================================================
# test_at_ranwtdttt.R
#
# Test script for at_ranwtdttt() – antithetic-sampling WTD estimator.
#
# Steps:
#   1. Load data (drugpakud.dta), filter to ATC = "N06AB06"
#   2. Display rxdate history for the first five individuals
#   3. Reproduce the antithetic index-date sampling for those five individuals
#      and display the resulting index dates (2 * nint per person)
#   4. Fit the model with at_ranwtdttt() (nint = 5)
#   5. Fit the model with ranwtdttt()    (nsamp = 10, matching total index dates)
#   6. Print estimates side-by-side for comparison
# =============================================================================

library(haven)        # read_dta()
library(data.table)   # fast data manipulation
library(wtdr)         # ranwtdttt(), at_ranwtdttt() – load package or source files below
library(bbmle)

## If wtdr is not yet installed, source the individual files instead:
#(HS: I could not make it work without these uncommented?)
source("R/wtd-class.R")
source("R/wtdttt.R")
source("R/ranwtdttt.R")
source("R/sandwich.R")
source("sandpit/at_ranwtdttt.R")

## ---------------------------------------------------------------------------
## Path setup
## Directory layout assumed:
##   wtdr/
##     sandpit/          <- this script lives here
##     inst/extdata/     <- drugpakud.dta lives here
##
## pkg_root resolves to the wtdr package root (one level above sandpit/).
## All file paths are built from pkg_root so the script works regardless of
## the working directory, as long as it is sourced from within the wtdr tree.
## ---------------------------------------------------------------------------

#(HS: perhaps something easier could be done here?)
pkg_root <- normalizePath(file.path(
  dirname(tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) file.path(getwd(), "dummy.R")   # fallback for interactive use
  )),
  "."
))

## Sanity-check: confirm we resolved to the wtdr package root
stopifnot(basename(pkg_root) == "wtdr")

data_file <- file.path(pkg_root, "inst", "extdata", "drugpakud.dta")

set.seed(42)   # reproducibility

# =============================================================================
# 1. Load and prepare data
# =============================================================================

df <- as.data.table(read_dta(data_file))

## Keep only the ATC code of interest and the two variables used in estimation
df <- df[atc == "N06AB06", .(id, rxdate)]

## rxdate is stored as a Stata date (days since 1960-01-01); haven imports it
## as a Date, but confirm:
stopifnot(inherits(df$rxdate, "Date"))

## Restrict to dispensings in 2015–2016.
## 2015 records provide look-back history for index dates placed in 2016;
## 2016 records fall within the sampling window itself.
df <- df[data.table::year(rxdate) %in% 2015:2016]

cat("Records after filtering to N06AB06 and years 2015-2016:", nrow(df), "\n")
cat("Unique individuals:                                     ", uniqueN(df$id), "\n\n")

# =============================================================================
# 2. Display rxdate history for the first five individuals
# =============================================================================

first5_ids <- unique(df$id)[1:5]

cat("----------------------------------------------------------------------\n")
cat("Observed dispensing dates (rxdate) for the first five individuals\n")
cat("----------------------------------------------------------------------\n")

## Display only the dispensings in 2015–2016 for readability
for (pid in first5_ids) {
  dates <- sort(df[id == pid, rxdate])
  cat(sprintf("  id = %s  (%d dispensing(s)): %s\n",
              pid, length(dates), paste(format(dates), collapse = ", ")))
}
cat("\n")

# =============================================================================
# 3. Reproduce antithetic index-date sampling for the first five individuals
#
#    This mirrors exactly what at_ranwtdttt() does internally, using the same
#    seed so the displayed dates correspond to those used in the estimation.
#
#    Parameters:
#      start   = 2016-01-01
#      end     = 2016-12-31
#      nint    = 5
#      reverse = TRUE   (reverse WTD: look back from index date)
# =============================================================================

start   <- as.Date("2016-01-01")
end     <- as.Date("2016-12-31")
nint    <- 5L
reverse <- TRUE

delta   <- as.numeric(end - start)   # 364 days
int_len <- delta / nint               # sub-interval length in days

cat("----------------------------------------------------------------------\n")
cat("Sampling window : ", format(start), "to", format(end), "\n")
cat("delta (days)    : ", delta,   "\n")
cat("nint            : ", nint,    "\n")
cat("int_len (days)  : ", int_len, "\n")
cat("Total index dates per individual: 2 * nint =", 2L * nint, "\n")
cat("----------------------------------------------------------------------\n\n")

## Build the offset table for the first five individuals only (for display)
off5 <- data.table(id = first5_ids, key = "id")

set.seed(42)   # same seed as the estimation call below

## Collect all antithetic index dates for the first five individuals
idx_list <- vector("list", nint)

for (k in seq_len(nint)) {
  a_num <- (k - 1L) * int_len

  U <- runif(n = nrow(off5))

  indda1 <- start + floor(a_num + U         * int_len)
  indda2 <- start + floor(a_num + (1 - U)   * int_len)

  idx_list[[k]] <- data.table(
    id       = rep(first5_ids, 2L),
    interval = k,
    variate  = rep(c("U", "1-U"), each = nrow(off5)),
    indda    = c(indda1, indda2)
  )
}

idx_df <- rbindlist(idx_list)
setorder(idx_df, id, interval, variate)

cat("----------------------------------------------------------------------\n")
cat("Antithetic index dates for the first five individuals\n")
cat("(interval k, both variates U and 1-U)\n")
cat("----------------------------------------------------------------------\n")

for (pid in first5_ids) {
  sub <- idx_df[id == pid]
  cat(sprintf("\n  id = %s\n", pid))
  cat(sprintf("    %-10s %-6s %s\n", "indda", "intv", "variate"))
  for (r in seq_len(nrow(sub))) {
    cat(sprintf("    %-10s  k=%d   %s\n",
                format(sub$indda[r]), sub$interval[r], sub$variate[r]))
  }
}
cat("\n")

# =============================================================================
# 4. Fit with at_ranwtdttt()  (nint = 5, 2*5 = 10 index dates per individual)
# =============================================================================

cat("======================================================================\n")
cat("Fitting model with at_ranwtdttt() – antithetic sampling, nint = 5\n")
cat("======================================================================\n")

set.seed(42)
fit_at <- at_ranwtdttt(
  data    = df,
  form    = rxdate ~ dlnorm(logitp, mu, lnsigma),
  id      = "id",
  start   = start,
  end     = end,
  reverse = reverse,
  nint    = 5L,
  robust  = TRUE
)

cat("\nSummary – at_ranwtdttt():\n")
print(summary(fit_at))

# =============================================================================
# 5. Fit with ranwtdttt()  (nsamp = 10, same total index dates per individual)
# =============================================================================

cat("======================================================================\n")
cat("Fitting model with ranwtdttt() – independent sampling, nsamp = 10\n")
cat("======================================================================\n")

set.seed(42)
fit_ran <- ranwtdttt(
  data    = df,
  form    = rxdate ~ dlnorm(logitp, mu, lnsigma),
  id      = "id",
  start   = start,
  end     = end,
  reverse = reverse,
  nsamp   = 10L,
  robust  = TRUE
)

cat("\nSummary – ranwtdttt():\n")
print(summary(fit_ran))

# =============================================================================
# 6. Side-by-side coefficient comparison
# =============================================================================

cat("\n======================================================================\n")
cat("Side-by-side comparison of estimates\n")
cat("(both methods: reverse WTD, lognormal, same window)\n")
cat("======================================================================\n")

## Extract coefficient tables
coef_at  <- coef(summary(fit_at))
coef_ran <- coef(summary(fit_ran))

## Align rows (parameter names should match)
params <- union(rownames(coef_at), rownames(coef_ran))

compare <- data.frame(
  Parameter   = params,
  AT_Estimate = coef_at [params, "Estimate"],
  AT_SE       = coef_at [params, "Std. Error"],
  RAN_Estimate= coef_ran[params, "Estimate"],
  RAN_SE      = coef_ran[params, "Std. Error"],
  Diff        = coef_at [params, "Estimate"] - coef_ran[params, "Estimate"],
  SE_gain     = (coef_at [params, "Std. Error"] - coef_ran[params, "Std. Error"]) / coef_ran[params, "Std. Error"] *100,
  row.names   = NULL
)

cat("\n")
print(compare, digits = 4)

cat("\n")
cat("Notes:\n")
cat("  AT  = at_ranwtdttt (antithetic sampling, nint = 5, 10 index dates/person)\n")
cat("  RAN = ranwtdttt    (independent sampling, nsamp = 10, 10 index dates/person)\n")
cat("  Diff = AT_Estimate - RAN_Estimate\n")
cat("  SE_gain = % reduction in SE (AT vs RAN)")
cat("  Estimates should be similar; SE differences reflect variance-reduction\n")
cat("  from antithetic sampling.\n")
