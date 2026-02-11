#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr group_by summarise mutate filter select across all_of n n_distinct .data
#' @importFrom tidyr pivot_wider
#' @importFrom stats pchisq qnorm
## usethis namespace: end
NULL

# Suppress R CMD check notes about global variables used in dplyr pipelines
utils::globalVariables(c(
  "strata", "a", "b", "c", "d", "m1", "m0", "n1", "n0",
  "events", "pyears", "stratum", "d0", "d1", "py0", "py1",
  "rate1", "rate0", "rr", "ln_rr", "var_ln_rr", "se_ln_rr",
  "ci_lower", "ci_upper", "or", "or_lower", "or_upper",
  "x", "py", "rate", "rate_per_1000"
))
