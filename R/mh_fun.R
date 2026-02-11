#' Format p-value for display
#'
#' Formats a p-value for display, showing "p<0.0001" when the value rounds to 0.
#'
#' @param p Numeric p-value.
#' @return Character string with formatted p-value.
#' @keywords internal
format_pvalue <- function(p) {
  if (is.na(p)) return(NA_character_)
  rounded <- round(p, 4)
  if (rounded == 0) {
    return("<0.0001")
  } else {
    return(as.character(rounded))
  }
}


#' Calculate Incidence Rates with Confidence Intervals
#'
#' Calculates incidence rates with confidence intervals for cohort or follow-up data.
#' Can compute overall rates or rates stratified by one or more grouping variables.
#' Uses the Stata strate method where SE(ln(rate)) = 1/sqrt(D).
#'
#' @param data A data frame containing the variables for analysis.
#' @param event Character string specifying the name of the event/case count variable.
#'   Should contain numeric values representing the number of events (typically 0 or 1
#'   for individual-level data, or counts for aggregated data).
#' @param time Character string specifying the name of the person-time variable.
#'   Should contain numeric values representing the follow-up time (e.g., person-years).
#' @param strata Character string or vector specifying the stratification variable(s).
#'   If \code{NULL} (the default), calculates the overall rate.
#' @param per Numeric value specifying the multiplier for rate display.
#'   Default is 1000 (rates per 1000 person-time units).
#' @param conf.level Numeric value between 0 and 1 specifying the confidence level
#'   for the confidence interval. Default is 0.95 (95% CI).
#'
#' @return Invisibly returns a data frame with the following columns:
#' \describe{
#'   \item{stratum}{Stratum identifier (if stratified).}
#'   \item{events}{Number of events.}
#'   \item{person_time}{Total person-time at risk.}
#'   \item{rate}{Incidence rate (per specified multiplier).}
#'   \item{ci_lower}{Lower bound of the confidence interval.}
#'   \item{ci_upper}{Upper bound of the confidence interval.}
#' }
#'
#' @details
#' The confidence interval is calculated using the method from Stata's strate command:
#' \itemize{
#'   \item SE(ln(rate)) = 1/sqrt(D) where D is the number of events
#'   \item CI = exp(ln(rate) ± z * SE(ln(rate)))
#' }
#'
#' This method produces confidence intervals that are always positive and have good
#' coverage properties, especially when the number of events is small.
#'
#' @seealso \code{\link{mh_analysis}} for rate ratio comparisons between groups.
#'
#' @examples
#' # Using the included mortality cohort dataset
#' data(mortality_cohort)
#'
#' # Calculate overall mortality rate
#' strate(mortality_cohort, "death", "person_years")
#'
#' # Calculate rates by treatment group
#' strate(mortality_cohort, "death", "person_years", strata = "treatment")
#'
#' # Calculate rates by severity
#' strate(mortality_cohort, "death", "person_years", strata = "severity")
#'
#' # Calculate rates per 100,000 person-years
#' strate(mortality_cohort, "death", "person_years", per = 100000)
#'
#' # Calculate rates stratified by multiple variables
#' strate(mortality_cohort, "death", "person_years", strata = c("treatment", "severity"))
#'
#' @export
strate <- function(data, event, time, strata = NULL, per = 1000, conf.level = 0.95) {

  # Check required variables exist
  required_vars <- c(event, time)
  if (!is.null(strata)) required_vars <- c(required_vars, strata)

  if (!all(required_vars %in% names(data))) {
    stop("One or more variables not found in data")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

  # Create grouping variable

  if (is.null(strata)) {
    data$.stratum <- "Overall"
    group_vars <- ".stratum"
  } else {
    group_vars <- strata
  }

  # Aggregate by strata
  result <- data |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      events = sum(.data[[event]], na.rm = TRUE),
      person_time = sum(.data[[time]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rate = events / person_time,
      se_ln_rate = ifelse(events > 0, 1 / sqrt(events), NA),
      ci_lower = ifelse(events > 0, exp(log(rate) - z * se_ln_rate), NA),
      ci_upper = ifelse(events > 0, exp(log(rate) + z * se_ln_rate), NA),
      # Scale to per multiplier
      rate_scaled = rate * per,
      ci_lower_scaled = ci_lower * per,
      ci_upper_scaled = ci_upper * per
    )

  # Output
  cat("\n")
  cat("Incidence Rates\n")
  cat("===============\n")
  cat("Event:", event, "\n")
  cat("Person-time:", time, "\n")
  if (!is.null(strata)) {
    cat("Stratified by:", paste(strata, collapse = ", "), "\n")
  }
  cat("Confidence level:", conf.level * 100, "%\n")
  cat("Rate per:", format(per, big.mark = ","), "\n")
  cat("\n")

  # Prepare output table
  if (is.null(strata)) {
    output_table <- result |>
      select(events, person_time, rate_scaled, ci_lower_scaled, ci_upper_scaled) |>
      mutate(
        person_time = round(person_time, 1),
        across(c(rate_scaled, ci_lower_scaled, ci_upper_scaled), \(x) round(x, 4))
      )
    names(output_table) <- c("Events", "Person-time",
                             paste0("Rate/", format(per, big.mark = ",")),
                             "CI lower", "CI upper")
  } else {
    output_table <- result |>
      select(all_of(group_vars), events, person_time, rate_scaled, ci_lower_scaled, ci_upper_scaled) |>
      mutate(
        person_time = round(person_time, 1),
        across(c(rate_scaled, ci_lower_scaled, ci_upper_scaled), \(x) round(x, 4))
      )
    rate_col_name <- paste0("Rate/", format(per, big.mark = ","))
    names(output_table) <- c(strata, "Events", "Person-time", rate_col_name, "CI lower", "CI upper")
  }

  print(as.data.frame(output_table), row.names = FALSE)

  # Calculate totals
  total_events <- sum(result$events)
  total_pt <- sum(result$person_time)
  total_rate <- total_events / total_pt
  total_se_ln <- ifelse(total_events > 0, 1 / sqrt(total_events), NA)
  total_ci_lower <- ifelse(total_events > 0, exp(log(total_rate) - z * total_se_ln) * per, NA)
  total_ci_upper <- ifelse(total_events > 0, exp(log(total_rate) + z * total_se_ln) * per, NA)

  if (!is.null(strata)) {
    cat("\n")
    cat("Total:", total_events, "events,", round(total_pt, 1), "person-time\n")
    cat("Overall rate:", round(total_rate * per, 4), "per", format(per, big.mark = ","), "\n")
    cat(conf.level * 100, "% CI: ", round(total_ci_lower, 4), " - ", round(total_ci_upper, 4), "\n", sep = "")
  }

  # Return data frame invisibly
  return_df <- result |>
    select(all_of(c(if (!is.null(strata)) group_vars else NULL,
                    "events", "person_time", "rate", "ci_lower", "ci_upper")))

  invisible(return_df)
}
#'
#' Calculates the Mantel-Haenszel odds ratio for case-control data. For binary
#' exposures, computes crude or stratified odds ratios with optional pooling.
#' For multi-level numeric exposures, performs a score test for trend.
#' Confidence intervals use the Robins-Breslow-Greenland variance estimator.
#'
#' @param data A data frame containing the variables for analysis.
#' @param exposure Character string specifying the name of the exposure variable.
#'   Can be binary (2-levels) or multi-level numeric. For binary exposures, the
#'   higher value is treated as "exposed". For multi-level exposures, a trend
#'   test is performed (exposure must be numeric).
#' @param outcome Character string specifying the name of the outcome variable.
#'   Must be a binary variable (exactly 2 unique values). The higher value is
#'   treated as "case" and the lower value as "control".
#' @param strata_vars Character string or vector of character strings specifying
#'   the stratification variable(s). If \code{NULL} (the default), a crude
#'   (unstratified) odds ratio is calculated. Multiple variables will be combined
#'   using \code{\link[base]{interaction}}. For trend tests, stratification
#'   computes u and v within each stratum and pools across strata.
#' @param conf.level Numeric value between 0 and 1 specifying the confidence level
#'   for the confidence interval. Default is 0.95 (95% CI).
#'
#' @return Invisibly returns a list with the following components:
#' \describe{
#'   \item{OR_MH}{The Mantel-Haenszel odds ratio (pooled, crude, or per-unit for trend).}
#'   \item{CI_lower}{Lower bound of the confidence interval.}
#'   \item{CI_upper}{Upper bound of the confidence interval.}
#'   \item{se_ln_or}{Standard error of the log odds ratio.}
#'   \item{chi2}{Chi-squared test statistic (association or trend).}
#'   \item{p_value}{P-value from the chi-squared test.}
#'   \item{n_total_strata}{Total number of strata formed (NA if unstratified or trend).}
#'   \item{n_informative}{Number of informative strata (NA if unstratified or trend).}
#'   \item{strata_data}{Data frame with stratum-specific cell counts (NULL if unstratified or trend).}
#'   \item{homogeneity_Q}{Q statistic for test of homogeneity (NA if unstratified or trend).}
#'   \item{homogeneity_df}{Degrees of freedom for homogeneity test (NA if unstratified or trend).}
#'   \item{homogeneity_p}{P-value for homogeneity test (NA if unstratified or trend).}
#'   \item{type}{Character string: "trend" if trend test was performed (NULL otherwise).}
#' }
#'
#' @details
#' The function automatically detects the coding of binary variables by sorting
#' unique values. The higher value is treated as the "positive" category (exposed
#' or case), and the lower value as the "negative" category (unexposed or control).
#'
#' For binary exposures with stratification, only informative strata (those
#' where at least one of ad or bc is positive) contribute to the pooled estimate.
#' A test of homogeneity (Breslow-Day type) is also performed to assess whether
#' the stratum-specific odds ratios are consistent.
#'
#' For multi-level numeric exposures, a Mantel extension score test for trend
#' is performed. The odds ratio represents the change in odds per one-unit
#' increase in the exposure variable.
#'
#' The confidence interval uses the Robins-Breslow-Greenland variance formula
#' for binary exposures, which matches output from Stata.
#'
#' @seealso \code{\link{stmh_r}} for stratified rate ratios,
#'   \code{\link{mh_analysis}} for a unified interface.
#'
#' @examples
#' # Using the included lung cancer case-control dataset
#' data(lung_cancer_cc)
#'
#' # Crude (unstratified) odds ratio
#' crude_result <- mh_analysis(lung_cancer_cc, "smoking", "case", measure = "or")
#'
#' # Stratified by age group
#' result <- mh_analysis(lung_cancer_cc, "smoking", "case", strata = "age_group", measure = "or")
#'
#' # Access the odds ratio
#' result$estimate
#'
#' @keywords internal
mh_or <- function(data, exposure, outcome, strata_vars = NULL, conf.level = 0.95) {

  # Check required variables exist
  required_vars <- c(exposure, outcome)
  if (!is.null(strata_vars)) required_vars <- c(required_vars, strata_vars)

  if (!all(required_vars %in% names(data))) {
    stop("One or more variables not found in data")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

  # Detect outcome coding
  out_levels <- sort(unique(data[[outcome]]))
  if (length(out_levels) != 2) stop("Outcome must have exactly 2 levels")

  # Determine if exposure is binary or multi-level
  exp_levels <- sort(unique(data[[exposure]]))
  n_exp_levels <- length(exp_levels)

  if (n_exp_levels < 2) {
    stop("Exposure variable must have at least 2 unique values.")
  }

  is_binary <- (n_exp_levels == 2)

  # For multi-level exposures, check if numeric for trend test
  if (!is_binary) {
    if (!is.numeric(data[[exposure]])) {
      stop("Multi-level exposure must be numeric for trend test.\n",
           "Found ", n_exp_levels, " levels in '", exposure, "': ",
           paste(head(exp_levels, 5), collapse = ", "),
           if (n_exp_levels > 5) "..." else "",
           "\nConvert to numeric scores (e.g., 1, 2, 3) for trend analysis.")
    }
    # Run trend test for odds ratio
    return(mh_or_trend(data, exposure, outcome, strata_vars, conf.level))
  }

  # If no strata, create a dummy stratum
  if (is.null(strata_vars)) {
    data$strata <- "Overall"
  } else {
    data$strata <- interaction(data[, strata_vars], drop = TRUE)
  }

  exp_yes <- exp_levels[2]
  exp_no  <- exp_levels[1]
  out_yes <- out_levels[2]
  out_no  <- out_levels[1]

  cat("Variable coding detected:\n")
  cat("  Exposure:", exposure, "-> Exposed =", exp_yes, ", Unexposed =", exp_no, "\n")
  cat("  Outcome:", outcome, "-> Case =", out_yes, ", Control =", out_no, "\n\n")

  # Get stratum-specific 2x2 cells
  results <- data |>
    group_by(strata) |>
    summarise(
      a = sum(.data[[exposure]] == exp_yes & .data[[outcome]] == out_yes),
      b = sum(.data[[exposure]] == exp_yes & .data[[outcome]] == out_no),
      c = sum(.data[[exposure]] == exp_no  & .data[[outcome]] == out_yes),
      d = sum(.data[[exposure]] == exp_no  & .data[[outcome]] == out_no),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(
      m1 = a + c,  # cases
      m0 = b + d,  # controls
      n1 = a + b,  # exposed
      n0 = c + d,  # unexposed
      # Stratum-specific OR and CI (Woolf method)
      or = (a * d) / (b * c),
      se_ln_or = sqrt(1/a + 1/b + 1/c + 1/d),
      or_lower = exp(log(or) - z * se_ln_or),
      or_upper = exp(log(or) + z * se_ln_or)
    ) |>
    filter(a * d > 0 | b * c > 0)

  n_total_strata <- n_distinct(data$strata)
  n_informative <- nrow(results)

  # Stata mhodds formulas for binary exposure (from mhodds.ado lines 197-269):
  # Using 2x2 table notation: a=d1 (exposed cases), b=h1 (exposed controls),
  #                          c=d2 (unexposed cases), d=h2 (unexposed controls)
  # q = a * d / n  (equivalent to d1 * h2 / pt)
  # r = b * c / n  (equivalent to d2 * h1 / pt)
  # u = q - r
  # v = n1 * n0 * m1 * m0 / (n^2 * (n - 1))  (equivalent to p1*p2*dt*ht/((pt-1)*pt^2))
  #
  # OR = q / r = sum(q) / sum(r) for stratified
  # chi2 = u^2 / v = (sum(u))^2 / sum(v) for stratified
  # CI: ef = exp(z * sqrt(v / (q * r))), then [OR/ef, OR*ef]
  #     For stratified: ef = exp(z * sqrt(sum(v) / (sum(q) * sum(r))))

  # Stratum-specific components
  q_i <- results$a * results$d / results$n
  r_i <- results$b * results$c / results$n
  u_i <- q_i - r_i
  v_i <- as.numeric(results$n1) * as.numeric(results$n0) *
         as.numeric(results$m1) * as.numeric(results$m0) /
         (as.numeric(results$n)^2 * (results$n - 1))

  # Pooled estimates
  q_total <- sum(q_i)
  r_total <- sum(r_i)
  u_total <- sum(u_i)
  v_total <- sum(v_i)

  # MH OR = q / r
  or_mh <- q_total / r_total

  # Chi-squared test: chi2 = u^2 / v
  chi2_mh <- u_total^2 / v_total
  p_value <- pchisq(chi2_mh, df = 1, lower.tail = FALSE)

  # Confidence interval: ef = exp(z * sqrt(v / (q * r)))
  # CI = [OR/ef, OR*ef] which equals exp(ln(OR) +/- z * sqrt(v/(q*r)))
  ef <- exp(z * sqrt(v_total / (q_total * r_total)))
  ci_lower <- or_mh / ef
  ci_upper <- or_mh * ef

  # Standard error of ln(OR): derived from CI formula
  # ln(OR) +/- z * sqrt(v/(q*r)) => se = sqrt(v/(q*r))
  se_ln_or <- sqrt(v_total / (q_total * r_total))

  # Output header
  cat("\n")
  cat("Mantel-Haenszel Odds Ratio Analysis\n")
  cat("====================================\n")
  cat("Exposure:", exposure, "\n")
  cat("  comparing", exposure, "=", exp_yes, "vs.", exposure, "=", exp_no, "\n")
  cat("Outcome:", outcome, "\n")
  if (!is.null(strata_vars)) {
    cat("Stratified by:", paste(strata_vars, collapse = ", "), "\n")
  }
  cat("Confidence level:", conf.level * 100, "%\n")
  cat("\n")

  # Calculate overall 2x2 table for display
  total_a <- sum(results$a)
  total_b <- sum(results$b)
  total_c <- sum(results$c)
  total_d <- sum(results$d)

  cat("=== 2x2 Table ===\n\n")
  table_display <- data.frame(
    Exposure = c("Exposed", "Unexposed", "Total"),
    Cases = c(total_a, total_c, total_a + total_c),
    Controls = c(total_b, total_d, total_b + total_d),
    Total = c(total_a + total_b, total_c + total_d, total_a + total_b + total_c + total_d)
  )
  print(table_display, row.names = FALSE)
  cat("\n")

  # Stratified output (only if strata provided)
  if (!is.null(strata_vars)) {

    cat("=== Stratum-Specific Odds Ratios ===\n\n")

    stratum_output <- results |>
      select(strata, a, b, c, d, or, or_lower, or_upper) |>
      mutate(
        across(c(or, or_lower, or_upper), \(x) round(x, 4))
      )

    print(as.data.frame(stratum_output), row.names = FALSE)

    cat("\n=== Mantel-Haenszel Pooled Estimate ===\n\n")
    cat("Note:", n_informative, "of", n_total_strata,
        "strata contribute information\n\n")

    # Homogeneity test (Stata mhodds formula - lines 248-258)
    # Using q_i, r_i, v_i computed above
    Q_homog <- sum((q_i * r_total - r_i * q_total)^2 / (q_total * r_total * v_i))
    df_homog <- n_informative - 1
    p_homog <- pchisq(Q_homog, df = df_homog, lower.tail = FALSE)

  } else {
    cat("=== Crude Odds Ratio ===\n\n")
    Q_homog <- NA
    df_homog <- NA
    p_homog <- NA
  }

  cat("OR:", round(or_mh, 4), "\n")
  cat(conf.level * 100, "% CI: ", round(ci_lower, 4), " - ", round(ci_upper, 4), "\n", sep = "")
  cat("SE(ln OR):", round(se_ln_or, 4), "\n")

  cat("\n=== Test of Association ===\n\n")
  cat("Chi2(1):", round(chi2_mh, 4), "\n")
  cat("p-value:", format_pvalue(p_value), "\n")

  # Homogeneity test (only if stratified)
  if (!is.null(strata_vars)) {
    cat("\n=== Test of Homogeneity ===\n\n")
    cat("Q statistic:", round(Q_homog, 3), "(df =", df_homog, ")\n")
    cat("p-value:", format_pvalue(p_homog), "\n")
  }

  # Return results
  invisible(list(
    OR_MH = or_mh,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    se_ln_or = se_ln_or,
    chi2 = chi2_mh,
    p_value = p_value,
    n_total_strata = if (!is.null(strata_vars)) n_total_strata else NA,
    n_informative = if (!is.null(strata_vars)) n_informative else NA,
    strata_data = if (!is.null(strata_vars)) results else NULL,
    homogeneity_Q = Q_homog,
    homogeneity_df = df_homog,
    homogeneity_p = p_homog
  ))
}


#' Mantel-Haenszel Trend Test for Odds Ratios (internal)
#'
#' Calculates a score test for trend in odds ratios across ordered exposure levels.
#' This is called internally by mh_or when the exposure has more than 2 levels.
#' Based on the Mantel extension test for trend. Matches Stata's mhodds command.
#'
#' @param data A data frame.
#' @param exposure Name of numeric exposure variable.
#' @param outcome Name of binary outcome variable.
#' @param strata_vars Name of stratification variable(s). When provided, computes
#'   u and v within each stratum, then sums across strata for the pooled estimate.
#'   A test of homogeneity is also performed.
#' @param conf.level Confidence level.
#'
#' @return A list with trend test results.
#' @keywords internal
mh_or_trend <- function(data, exposure, outcome, strata_vars = NULL, conf.level = 0.95) {

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

  # Detect outcome coding
  out_levels <- sort(unique(data[[outcome]]))
  out_yes <- out_levels[2]
  out_no <- out_levels[1]

  # Create stratum variable if stratification requested
  if (!is.null(strata_vars)) {
    data$strata <- interaction(data[, strata_vars, drop = FALSE], drop = TRUE)
  } else {
    data$strata <- "Overall"
  }

  # Aggregate by stratum and exposure level
  # Stata mhodds trend test formulas (from mhodds.ado lines 273-343):
  # Within each stratum:
  # dt = sum(d) = total cases in stratum
  # pt = sum(n) = total sample size in stratum
  # d1 = sum(d * x) = weighted sum of cases by exposure score
  # p1 = sum(n * x) = weighted sum of total by exposure score
  # p2 = sum(n * x^2) = weighted sum of total by exposure score squared
  #
  # u = d1 - (dt * p1 / pt)  - score statistic
  # v = dt * (pt - dt) * (p2 - p1*p1/pt) / (pt * (pt - 1))  - variance WITH finite pop correction
  #
  # For stratified: sum u and v across strata
  # OR = exp(sum(u)/sum(v))
  # CI: exp(ln(OR) +/- z/sqrt(sum(v)))

  # Compute stratum-specific u and v
  strata_stats <- data |>
    group_by(strata, .data[[exposure]]) |>
    summarise(
      cases = sum(.data[[outcome]] == out_yes),
      controls = sum(.data[[outcome]] == out_no),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(x = .data[[exposure]]) |>
    group_by(strata) |>
    summarise(
      dt = sum(cases),           # total cases in stratum
      pt = sum(n),               # total sample size in stratum
      d1 = sum(cases * x),       # sum(d * x)
      p1 = sum(n * x),           # sum(n * x)
      p2 = sum(n * x^2),         # sum(n * x^2)
      .groups = "drop"
    ) |>
    mutate(
      u = d1 - (dt * p1 / pt),
      v = dt * (pt - dt) * (p2 - p1 * p1 / pt) / (pt * (pt - 1))
    ) |>
    filter(v > 0 & !is.na(v))  # Only informative strata

  # Overall exposure-level frequencies (for display)
  agg <- data |>
    group_by(.data[[exposure]]) |>
    summarise(
      cases = sum(.data[[outcome]] == out_yes),
      controls = sum(.data[[outcome]] == out_no),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(
      x = .data[[exposure]],
      prop_cases = cases / n
    )

  # Pooled estimates (sum across strata)
  u_total <- sum(strata_stats$u)
  v_total <- sum(strata_stats$v)

  # Chi-squared test statistic: chi2 = u^2 / v
  chi2_trend <- u_total^2 / v_total
  p_trend <- pchisq(chi2_trend, df = 1, lower.tail = FALSE)

  # Odds ratio estimate: OR = exp(u/v)
  ln_or <- u_total / v_total
  or_trend <- exp(ln_or)

  # Confidence interval: exp(ln_or +/- z / sqrt(v))
  ci_lower <- exp(ln_or - z / sqrt(v_total))
  ci_upper <- exp(ln_or + z / sqrt(v_total))

  # Standard error of ln(OR): se = 1/sqrt(v)
  se_ln_or <- 1 / sqrt(v_total)

  # Homogeneity test (only if stratified with multiple informative strata)
  n_total_strata <- n_distinct(data$strata)
  n_informative <- nrow(strata_stats)

  if (!is.null(strata_vars) && n_informative > 1) {
    # Homogeneity statistic: sum(u_i^2/v_i) - (sum(u_i))^2/sum(v_i)
    Q_homog <- sum(strata_stats$u^2 / strata_stats$v) - (u_total^2 / v_total)
    df_homog <- n_informative - 1
    p_homog <- pchisq(Q_homog, df = df_homog, lower.tail = FALSE)
  } else {
    Q_homog <- NA
    df_homog <- NA
    p_homog <- NA
  }

  # Output
  cat("\n")
  cat("Score Test for Trend of Odds Ratios\n")
  cat("====================================\n")
  cat("Outcome:", outcome, "(Case =", out_yes, ", Control =", out_no, ")\n")
  cat("Exposure:", exposure, "(", length(unique(agg$x)), "levels )\n")
  if (!is.null(strata_vars)) {
    cat("Stratified by:", paste(strata_vars, collapse = ", "), "\n")
  }
  cat("Confidence level:", conf.level * 100, "%\n")
  cat("\n")

  cat("=== Exposure-Specific Frequencies ===\n\n")

  freq_output <- agg |>
    select(x, cases, controls, n) |>
    mutate(
      pct_cases = round(cases / n * 100, 1)
    )

  names(freq_output) <- c(exposure, "Cases", "Controls", "Total", "% Cases")
  print(as.data.frame(freq_output), row.names = FALSE)

  if (!is.null(strata_vars)) {
    cat("\n=== Stratum-Specific Score Statistics ===\n\n")
    cat("Note:", n_informative, "of", n_total_strata,
        "strata contribute information\n\n")

    strata_output <- strata_stats |>
      mutate(
        or_i = exp(u / v),
        across(c(u, v, or_i), \(x) round(x, 4))
      ) |>
      select(strata, dt, pt, u, v, or_i)

    names(strata_output) <- c("Stratum", "Cases", "Total", "u", "v", "OR")
    print(as.data.frame(strata_output), row.names = FALSE)

    cat("\n=== Mantel-Haenszel Pooled Trend Estimate ===\n\n")
  } else {
    cat("\n=== Mantel-Haenszel Trend Estimate ===\n\n")
  }

  cat("Odds ratio per unit increase:", round(or_trend, 4), "\n")
  cat(conf.level * 100, "% CI: ", round(ci_lower, 4), " - ", round(ci_upper, 4), "\n", sep = "")
  cat("SE(ln OR):", round(se_ln_or, 4), "\n")
  cat("\n")
  cat("Note: The odds ratio is for a one-unit increase in", exposure, "\n")

  cat("\n=== Score Test for Trend ===\n\n")
  cat("Chi2(1):", round(chi2_trend, 2), "\n")
  cat("p-value:", format_pvalue(p_trend), "\n")

  # Homogeneity test (only if stratified)
  if (!is.null(strata_vars) && n_informative > 1) {
    cat("\n=== Test of Homogeneity ===\n\n")
    cat("Q statistic:", round(Q_homog, 3), "(df =", df_homog, ")\n")
    cat("p-value:", format_pvalue(p_homog), "\n")
  } else if (!is.null(strata_vars) && n_informative <= 1) {
    cat("\nNote: Too few informative strata to test for homogeneity.\n")
  }

  # Return results
  invisible(list(
    exposure_levels = agg,
    OR_MH = or_trend,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    se_ln_or = se_ln_or,
    chi2 = chi2_trend,
    p_value = p_trend,
    n_total_strata = if (!is.null(strata_vars)) n_total_strata else NA,
    n_informative = if (!is.null(strata_vars)) n_informative else NA,
    strata_data = if (!is.null(strata_vars)) strata_stats else NULL,
    homogeneity_Q = Q_homog,
    homogeneity_df = df_homog,
    homogeneity_p = p_homog,
    type = "trend"
  ))
}


#' Stratified Mantel-Haenszel Rate Ratio
#'
#' Calculates the Mantel-Haenszel pooled incidence rate ratio (IRR) for stratified
#' cohort or follow-up data. This function aggregates events and person-time by
#' stratum and exposure status, computes stratum-specific rate ratios, and provides
#' the pooled estimate with confidence intervals using the Greenland-Robins variance
#' estimator.
#'
#' @param data A data frame containing the variables for analysis.
#' @param event Character string specifying the name of the event/case count variable.
#'   Should contain numeric values representing the number of events (typically 0 or 1
#'   for individual-level data, or counts for aggregated data).
#' @param exposure Character string specifying the name of the exposure variable.
#'   Should be a binary variable coded as 0 (unexposed) and 1 (exposed).
#' @param time Character string specifying the name of the person-time variable.
#'   Should contain numeric values representing the follow-up time (e.g., person-years).
#' @param strata Character string or vector specifying the stratification variable(s).
#'   If \code{NULL} (the default), a crude (unstratified) rate ratio is calculated.
#' @param conf.level Numeric value between 0 and 1 specifying the confidence level
#'   for the confidence interval. Default is 0.95 (95% CI).
#'
#' @return Invisibly returns a list with the following components:
#' \describe{
#'   \item{strata}{Data frame with stratum-specific results (NULL if unstratified).}
#'   \item{irr}{The Mantel-Haenszel pooled incidence rate ratio.}
#'   \item{se}{Standard error of the log rate ratio.}
#'   \item{ci}{Numeric vector of length 2 containing the confidence interval bounds.}
#'   \item{homogeneity_Q}{Q statistic for the test of homogeneity (NA if unstratified).}
#'   \item{homogeneity_df}{Degrees of freedom for the homogeneity test.}
#'   \item{homogeneity_p}{P-value for the test of homogeneity.}
#' }
#'
#' @details
#' The function assumes that the exposure variable is coded as 0 (unexposed) and
#' 1 (exposed). Events and person-time are aggregated within each stratum-exposure
#' combination.
#'
#' The pooled rate ratio is calculated using the Mantel-Haenszel method, and the
#' variance is estimated using the Greenland-Robins formula. When stratification
#' variables are provided, a test of homogeneity (Q statistic) is also performed
#' to assess whether the stratum-specific rate ratios are consistent.
#'
#' A warning is issued if any stratum has zero events in either the exposed or
#' unexposed group, as this can lead to unreliable confidence intervals.
#'
#' @seealso \code{\link{mh_or}} for stratified odds ratios,
#'   \code{\link{mh_analysis}} for a unified interface.
#'
#' @examples
#' # Using the included mortality cohort dataset
#' data(mortality_cohort)
#'
#' # Crude rate ratio (no stratification)
#' crude_result <- mh_analysis(mortality_cohort, "treatment", "death",
#'                             time = "person_years", measure = "irr")
#'
#' # Stratified by disease severity
#' strat_result <- mh_analysis(mortality_cohort, "treatment", "death",
#'                             time = "person_years", strata = "severity", measure = "irr")
#'
#' # Access the rate ratio
#' strat_result$estimate
#'
#' @keywords internal
stmh_r <- function(data, event, exposure, time, strata = NULL, conf.level = 0.95) {

  # Check required variables exist
  required_vars <- c(event, exposure, time)
  if (!is.null(strata)) required_vars <- c(required_vars, strata)

  if (!all(required_vars %in% names(data))) {
    stop("One or more variables not found in data")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

 # Determine if exposure is binary or multi-level
  exp_values <- sort(unique(data[[exposure]]))
  n_levels <- length(exp_values)

  if (n_levels < 2) {
    stop("Exposure variable must have at least 2 unique values.")
  }

  is_binary <- (n_levels == 2)

  # For multi-level exposures, check if numeric for trend test
  if (!is_binary) {
    if (!is.numeric(data[[exposure]])) {
      stop("Multi-level exposure must be numeric for trend test.\n",
           "Found ", n_levels, " levels in '", exposure, "': ",
           paste(head(exp_values, 5), collapse = ", "),
           if(n_levels > 5) "..." else "",
           "\nConvert to numeric scores (e.g., 1, 2, 3) for trend analysis.")
    }
    # Run trend test
    return(stmh_r_trend(data, event, exposure, time, strata, conf.level))
  }

  # Detect event coding

  event_values <- sort(unique(data[[event]]))
  event_yes <- event_values[length(event_values)]  # highest value = event
  event_no <- event_values[1]                       # lowest value = no event

  # Binary exposure: check if coded as 0/1
  exp_levels <- sort(exp_values)
  exp_yes <- exp_levels[2]  # Store original labels before recoding
 exp_no <- exp_levels[1]

  if (!all(exp_values %in% c(0, 1))) {
    # Convert to 0/1
    data[[exposure]] <- as.integer(data[[exposure]] == exp_levels[2])
  }

  # Print variable coding
  cat("Variable coding detected:\n")
  cat("  Exposure:", exposure, "-> Exposed =", exp_yes, ", Unexposed =", exp_no, "\n")
  cat("  Event:", event, "-> Event =", event_yes, ", No event =", event_no, "\n\n")

  # If no strata, create a dummy stratum
  if (is.null(strata)) {
    data$.stratum <- "Overall"
    strata_var <- ".stratum"
  } else {
    strata_var <- strata
  }

  # Aggregate by stratum and exposure
  agg <- data |>
    group_by(across(all_of(c(strata_var, exposure)))) |>
    summarise(
      events = sum(.data[[event]], na.rm = TRUE),
      pyears = sum(.data[[time]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = all_of(exposure),
      values_from = c(events, pyears),
      names_glue = "{.value}_{.name}"
    )

  # Rename columns (assuming exposure coded 1/0)
  names(agg) <- c("stratum", "d0", "d1", "py0", "py1")

  agg <- agg |>
    select(stratum, d1, py1, d0, py0)

  # Warn about zero cells
  if (any(agg$d1 == 0) | any(agg$d0 == 0)) {
    warning("Zero events in some strata. CIs may be unreliable.")
  }

  # Stratum-specific estimates
  agg <- agg |>
    mutate(
      rate1 = d1 / py1,
      rate0 = d0 / py0,
      rr = rate1 / rate0,
      ln_rr = log(rr),
      # Correct variance for stratum-specific IRR: Var(ln(IRR)) = 1/d1 + 1/d0
      var_ln_rr = 1/d1 + 1/d0,
      se_ln_rr = sqrt(var_ln_rr),
      ci_lower = exp(ln_rr - z * se_ln_rr),
      ci_upper = exp(ln_rr + z * se_ln_rr)
    )

  # Pooled MH estimate
  T_i <- agg$py1 + agg$py0
  num <- sum((agg$d1 * agg$py0) / T_i)
  den <- sum((agg$d0 * agg$py1) / T_i)
  rr_mh <- num / den

  # Greenland-Robins variance
  var_mh <- sum((agg$d1 + agg$d0) * agg$py1 * agg$py0 / T_i^2) / (num * den)
  se_mh <- sqrt(var_mh)

  ci_mh_lower <- exp(log(rr_mh) - z * se_mh)
  ci_mh_upper <- exp(log(rr_mh) + z * se_mh)

  # Output header
  cat("\n")
  cat("Mantel-Haenszel Rate Ratio Analysis\n")
  cat("====================================\n")
  cat("Event:", event, "\n")
  cat("Exposure:", exposure, "\n")
  cat("  comparing", exposure, "=", exp_yes, "vs.", exposure, "=", exp_no, "\n")
  if (!is.null(strata)) {
    cat("Stratified by:", strata, "\n")
  }
  cat("Confidence level:", conf.level * 100, "%\n")
  cat("\n")

  # Calculate overall rates by exposure (for display)
  # With confidence intervals using Stata strate method:
  # SE(ln(rate)) = 1/sqrt(D), CI = exp(ln(rate) +/- z * SE)
  total_d1 <- sum(agg$d1)
  total_d0 <- sum(agg$d0)
  total_py1 <- sum(agg$py1)
  total_py0 <- sum(agg$py0)

  rate0 <- total_d0 / total_py0
  rate1 <- total_d1 / total_py1

  # CI for rates (Stata strate method)
  se_ln_rate0 <- ifelse(total_d0 > 0, 1 / sqrt(total_d0), NA)
  se_ln_rate1 <- ifelse(total_d1 > 0, 1 / sqrt(total_d1), NA)

  rate0_ci_lower <- ifelse(total_d0 > 0, exp(log(rate0) - z * se_ln_rate0), NA)
  rate0_ci_upper <- ifelse(total_d0 > 0, exp(log(rate0) + z * se_ln_rate0), NA)
  rate1_ci_lower <- ifelse(total_d1 > 0, exp(log(rate1) - z * se_ln_rate1), NA)
  rate1_ci_upper <- ifelse(total_d1 > 0, exp(log(rate1) + z * se_ln_rate1), NA)

  cat("=== Exposure-Specific Rates ===\n\n")
  rate_table <- data.frame(
    Exposure = c(exp_no, exp_yes),
    Events = c(total_d0, total_d1),
    `Person-time` = round(c(total_py0, total_py1), 1),
    Rate = round(c(rate0, rate1), 4),
    `CI lower` = round(c(rate0_ci_lower, rate1_ci_lower), 4),
    `CI upper` = round(c(rate0_ci_upper, rate1_ci_upper), 4),
    check.names = FALSE
  )
  names(rate_table)[1] <- exposure
  print(rate_table, row.names = FALSE)
  cat("\n")

  # Stratified output (only if strata provided)
  if (!is.null(strata)) {

    cat("=== Stratum-Specific Rate Ratios ===\n\n")

    stratum_output <- agg |>
      select(stratum, d1, py1, d0, py0, rr, ci_lower, ci_upper) |>
      mutate(
        across(c(rr, ci_lower, ci_upper), \(x) round(x, 4)),
        across(c(py1, py0), \(x) round(x, 1))
      )

    print(as.data.frame(stratum_output), row.names = FALSE)

    cat("\n=== Mantel-Haenszel Pooled Estimate ===\n\n")

    # Homogeneity test (Stata stmh formula - lines 186-194 of stmh.ado)
    # q_i = d1_i * py0_i / T_i
    # r_i = d0_i * py1_i / T_i
    # v_i = (d1_i + d0_i) * py1_i * py0_i / T_i^2
    # Q_total = sum(q_i), R_total = sum(r_i)
    # Homogeneity stat = sum((q_i * R_total - r_i * Q_total)^2 / (Q_total * R_total * v_i))

    q_i <- agg$d1 * agg$py0 / T_i
    r_i <- agg$d0 * agg$py1 / T_i
    v_i <- (agg$d1 + agg$d0) * agg$py1 * agg$py0 / T_i^2

    Q_total <- sum(q_i)  # This equals 'num' computed above
    R_total <- sum(r_i)  # This equals 'den' computed above

    # Homogeneity Q statistic (Breslow-Day type test)
    Q <- sum((q_i * R_total - r_i * Q_total)^2 / (Q_total * R_total * v_i))
    df_q <- nrow(agg) - 1
    p_homog <- pchisq(Q, df = df_q, lower.tail = FALSE)

  } else {
    cat("=== Crude Rate Ratio ===\n\n")
    Q <- NA
    df_q <- NA
    p_homog <- NA
  }

  cat("IRR:", round(rr_mh, 4), "\n")
  cat("95% CI:", round(ci_mh_lower, 4), "-", round(ci_mh_upper, 4), "\n")
  cat("SE(ln IRR):", round(se_mh, 4), "\n")

  # MH Chi-squared test for association (score test) - matches Stata
  # Under H0: IRR = 1, E(d1) = (d1 + d0) * py1 / T
  # Var(d1) = (d1 + d0) * py1 * py0 / T^2
  observed_d1 <- sum(agg$d1)
  expected_d1 <- sum((agg$d1 + agg$d0) * agg$py1 / T_i)
  var_d1 <- sum((agg$d1 + agg$d0) * agg$py1 * agg$py0 / T_i^2)

  chi2_test <- (observed_d1 - expected_d1)^2 / var_d1
  p_chi2 <- pchisq(chi2_test, df = 1, lower.tail = FALSE)

  cat("\n=== Test of Association ===\n\n")
  cat("Chi2(1):", round(chi2_test, 4), "\n")
  cat("p-value:", format_pvalue(p_chi2), "\n")

  # Homogeneity test (only if stratified)
  if (!is.null(strata)) {
    cat("\n=== Test of Homogeneity ===\n\n")
    cat("Q statistic:", round(Q, 3), "(df =", df_q, ")\n")
    cat("p-value:", format_pvalue(p_homog), "\n")
  }

  # Return results
  invisible(list(
    strata = if (!is.null(strata)) agg else NULL,
    irr = rr_mh,
    se = se_mh,
    ci = c(ci_mh_lower, ci_mh_upper),
    chi2 = chi2_test,
    chi2_p = p_chi2,
    homogeneity_Q = Q,
    homogeneity_df = df_q,
    homogeneity_p = p_homog
  ))
}


#' Mantel-Haenszel Trend Test for Rate Ratios (internal)
#'
#' Calculates a score test for trend in rates across ordered exposure levels.
#' This is called internally by stmh_r when the exposure has more than 2 levels.
#' Matches Stata's stmh command exactly.
#'
#' @param data A data frame.
#' @param event Name of event variable.
#' @param exposure Name of numeric exposure variable.
#' @param time Name of person-time variable.
#' @param strata Name of stratification variable(s). When provided, computes
#'   u and v within each stratum, then sums across strata for the pooled estimate.
#'   A test of homogeneity is also performed.
#' @param conf.level Confidence level.
#'
#' @return A list with trend test results.
#' @keywords internal
stmh_r_trend <- function(data, event, exposure, time, strata = NULL, conf.level = 0.95) {

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

  # Create stratum variable if stratification requested
  if (!is.null(strata)) {
    data$.strata <- interaction(data[, strata, drop = FALSE], drop = TRUE)
  } else {
    data$.strata <- "Overall"
  }

  # Stata stmh trend test formulas (from stmh.ado lines 209-273):
  # Within each stratum:
  # dt = sum(d) = total events in stratum
  # yt = sum(py) = total person-time in stratum
  # d1 = sum(d * x) = weighted sum of events by exposure score
  # y1 = sum(py * x) = weighted sum of person-time by exposure score
  # y2 = sum(py * x^2) = weighted sum of person-time by exposure score squared
  #
  # u = d1 - dt * y1 / yt     - score statistic
  # v = dt * (y2 - y1^2/yt) / yt  - variance (NO finite population correction!)
  #
  # For stratified: sum u and v across strata
  # RR = exp(sum(u)/sum(v))
  # CI: exp(ln(RR) +/- z/sqrt(sum(v)))

  # Compute stratum-specific u and v
  strata_stats <- data |>
    group_by(.strata, .data[[exposure]]) |>
    summarise(
      d = sum(.data[[event]], na.rm = TRUE),
      py = sum(.data[[time]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(x = .data[[exposure]]) |>
    group_by(.strata) |>
    summarise(
      dt = sum(d),           # total events in stratum
      yt = sum(py),          # total person-time in stratum
      d1 = sum(d * x),       # sum(d * x)
      y1 = sum(py * x),      # sum(py * x)
      y2 = sum(py * x^2),    # sum(py * x^2)
      .groups = "drop"
    ) |>
    mutate(
      u = d1 - dt * y1 / yt,
      v = dt * (y2 - y1^2 / yt) / yt
    ) |>
    filter(v > 0 & !is.na(v))  # Only informative strata

  # Overall exposure-level frequencies (for display)
  agg <- data |>
    group_by(.data[[exposure]]) |>
    summarise(
      d = sum(.data[[event]], na.rm = TRUE),
      py = sum(.data[[time]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      x = .data[[exposure]],
      rate = d / py
    )

  # Pooled estimates (sum across strata)
  u_total <- sum(strata_stats$u)
  v_total <- sum(strata_stats$v)

  # Chi-squared test statistic: chi2 = u^2 / v
  chi2_trend <- u_total^2 / v_total
  p_trend <- pchisq(chi2_trend, df = 1, lower.tail = FALSE)

  # Rate ratio estimate: RR = exp(u/v)
  ln_rr <- u_total / v_total
  rr_trend <- exp(ln_rr)

  # Confidence interval: exp(ln_rr +/- z / sqrt(v))
  ci_lower <- exp(ln_rr - z / sqrt(v_total))
  ci_upper <- exp(ln_rr + z / sqrt(v_total))

  # Standard error of ln(RR): se = 1/sqrt(v)
  se_ln_rr <- 1 / sqrt(v_total)

  # Homogeneity test (only if stratified with multiple informative strata)
  n_total_strata <- n_distinct(data$.strata)
  n_informative <- nrow(strata_stats)

  if (!is.null(strata) && n_informative > 1) {
    # Homogeneity statistic: sum(u_i^2/v_i) - (sum(u_i))^2/sum(v_i)
    Q_homog <- sum(strata_stats$u^2 / strata_stats$v) - (u_total^2 / v_total)
    df_homog <- n_informative - 1
    p_homog <- pchisq(Q_homog, df = df_homog, lower.tail = FALSE)
  } else {
    Q_homog <- NA
    df_homog <- NA
    p_homog <- NA
  }

  # Output
  cat("\n")
  cat("Score Test for Trend of Rates\n")
  cat("==============================\n")
  cat("Event:", event, "\n")
  cat("Exposure:", exposure, "(", length(unique(agg$x)), "levels )\n")
  if (!is.null(strata)) {
    cat("Stratified by:", paste(strata, collapse = ", "), "\n")
  }
  cat("Confidence level:", conf.level * 100, "%\n")
  cat("\n")

  cat("=== Exposure-Specific Rates ===\n\n")

  # Add CIs for rates using Stata strate method: SE(ln(rate)) = 1/sqrt(D)
  rate_output <- agg |>
    select(x, d, py, rate) |>
    mutate(
      se_ln_rate = ifelse(d > 0, 1 / sqrt(d), NA),
      rate_ci_lower = ifelse(d > 0, exp(log(rate) - z * se_ln_rate), NA),
      rate_ci_upper = ifelse(d > 0, exp(log(rate) + z * se_ln_rate), NA),
      rate_per_1000 = round(rate * 1000, 2),
      ci_lower_1000 = round(rate_ci_lower * 1000, 2),
      ci_upper_1000 = round(rate_ci_upper * 1000, 2),
      py = round(py, 1)
    ) |>
    select(x, d, py, rate_per_1000, ci_lower_1000, ci_upper_1000)

  names(rate_output) <- c(exposure, "Events", "Person-time", "Rate/1000", "CI lower", "CI upper")
  print(as.data.frame(rate_output), row.names = FALSE)

  if (!is.null(strata)) {
    cat("\n=== Stratum-Specific Score Statistics ===\n\n")
    cat("Note:", n_informative, "of", n_total_strata,
        "strata contribute information\n\n")

    strata_output <- strata_stats |>
      mutate(
        rr_i = exp(u / v),
        across(c(u, v, rr_i), \(x) round(x, 4)),
        yt = round(yt, 1)
      ) |>
      select(.strata, dt, yt, u, v, rr_i)

    names(strata_output) <- c("Stratum", "Events", "Person-time", "u", "v", "RR")
    print(as.data.frame(strata_output), row.names = FALSE)

    cat("\n=== Mantel-Haenszel Pooled Trend Estimate ===\n\n")
  } else {
    cat("\n=== Mantel-Haenszel Trend Estimate ===\n\n")
  }

  cat("Rate ratio per unit increase:", round(rr_trend, 4), "\n")
  cat(conf.level * 100, "% CI: ", round(ci_lower, 4), " - ", round(ci_upper, 4), "\n", sep = "")
  cat("\n")
  cat("Note: The rate ratio is an approximation to the rate ratio\n")
  cat("      for a one-unit increase in", exposure, "\n")

  cat("\n=== Score Test for Trend ===\n\n")
  cat("Chi2(1):", round(chi2_trend, 2), "\n")
  cat("p-value:", format_pvalue(p_trend), "\n")

  # Homogeneity test (only if stratified)
  if (!is.null(strata) && n_informative > 1) {
    cat("\n=== Test of Homogeneity ===\n\n")
    cat("Q statistic:", round(Q_homog, 3), "(df =", df_homog, ")\n")
    cat("p-value:", format_pvalue(p_homog), "\n")
  } else if (!is.null(strata) && n_informative <= 1) {
    cat("\nNote: Too few informative strata to test for homogeneity.\n")
  }

  # Return results
  invisible(list(
    exposure_levels = agg,
    irr = rr_trend,
    se = se_ln_rr,
    ci = c(ci_lower, ci_upper),
    chi2_trend = chi2_trend,
    p_trend = p_trend,
    n_total_strata = if (!is.null(strata)) n_total_strata else NA,
    n_informative = if (!is.null(strata)) n_informative else NA,
    strata_data = if (!is.null(strata)) strata_stats else NULL,
    homogeneity_Q = Q_homog,
    homogeneity_df = df_homog,
    homogeneity_p = p_homog,
    type = "trend"
  ))
}


#' Mantel-Haenszel Analysis
#'
#' A unified interface for performing Mantel-Haenszel stratified analyses,
#' supporting both odds ratios and incidence rate ratios . This function provides a consistent interface and returns standardised output.
#'
#' @param data A data frame containing the variables for analysis.
#' @param exposure Character string specifying the name of the exposure variable.
#' @param outcome Character string specifying the name of the outcome variable.
#'   For odds ratio analysis (\code{measure = "or"}), this is the case/control
#'   indicator. For rate ratio analysis (\code{measure = "irr"}), this is the
#'   event/failure indicator.
#' @param time Character string specifying the name of the person-time variable.
#'   Required when \code{measure = "irr"}. Ignored when \code{measure = "or"}.
#' @param strata Character string or vector specifying the stratification variable(s).
#'   Optional for both measures (crude analysis if not provided).
#' @param measure Character string specifying the type of analysis: \code{"or"} for
#'   odds ratio (default) or \code{"irr"} for incidence rate ratio.
#' @param conf.level Numeric value between 0 and 1 specifying the confidence level.
#'   Default is 0.95 (95% CI).
#'
#' @return Invisibly returns an object of class \code{"mh_analysis"} containing:
#' \describe{
#'   \item{measure}{Character string describing the measure type.}
#'   \item{estimate}{The pooled Mantel-Haenszel estimate.}
#'   \item{ci}{Numeric vector of length 2 with confidence interval bounds.}
#'   \item{se_ln}{Standard error of the log-transformed estimate.}
#'   \item{test_stat}{Test statistic (homogeneity Q for stratified analyses).}
#'   \item{test_name}{Name of the test performed.}
#'   \item{test_df}{Degrees of freedom for the homogeneity test.}
#'   \item{p_value}{P-value from the homogeneity test.}
#'   \item{chi2}{Chi-squared test statistic for association.}
#'   \item{chi2_p}{P-value for the association test.}
#'   \item{conf.level}{The confidence level used.}
#'   \item{strata_data}{Data frame with stratum-specific results.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @details
#' This function serves as a wrapper around \code{\link{mh_or}} and
#' \code{\link{stmh_r}}, providing a unified interface for both types of analysis.
#' The output is returned as an S3 object with a custom print method for clean
#' display.
#'
#' For odds ratio analysis (\code{measure = "or"}):
#' \itemize{
#'   \item \code{outcome} is the case/control indicator (required)
#'   \item \code{strata} is optional (crude analysis if not provided)
#'   \item \code{time} is ignored with a warning if provided
#'   \item Uses \code{\link{mh_or}} internally
#' }
#'
#' For rate ratio analysis (\code{measure = "irr"}):
#' \itemize{
#'   \item \code{outcome} is the event/failure indicator (required)
#'   \item \code{time} is required (person-time at risk)
#'   \item \code{strata} is optional (crude analysis if not provided)
#'   \item Uses \code{\link{stmh_r}} internally
#' }
#'
#' @seealso \code{\link{mh_or}} for direct odds ratio calculation,
#'   \code{\link{stmh_r}} for direct rate ratio calculation,
#'   \code{\link{print.mh_analysis}} for the print method.
#'
#' @examples
#' # Example for odds ratio analysis using lung cancer data
#' data(lung_cancer_cc)
#'
#' or_result <- mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case",
#'                          strata = "age_group", measure = "or")
#' print(or_result)
#'
#' # Example for rate ratio analysis using mortality cohort data
#' data(mortality_cohort)
#'
#' irr_result <- mh_analysis(mortality_cohort, exposure = "treatment", outcome = "death",
#'                           time = "person_years", strata = "severity", measure = "irr")
#' print(irr_result)
#'
#' @export
mh_analysis <- function(data,
                        exposure,
                        outcome,
                        time = NULL,
                        strata = NULL,
                        measure = c("or", "irr"),
                        conf.level = 0.95) {

  measure <- match.arg(measure)

  # Validate inputs based on measure type
  if (measure == "or") {

    if (!is.null(time)) {
      warning("'time' is ignored for odds ratio analysis")
    }

    # Call mh_or function
    result <- mh_or(
      data = data,
      exposure = exposure,
      outcome = outcome,
      strata_vars = strata,
      conf.level = conf.level
    )

    # Standardise output
    out <- list(
      measure = "Odds Ratio",
      estimate = result$OR_MH,
      ci = c(result$CI_lower, result$CI_upper),
      se_ln = result$se_ln_or,
      test_stat = result$homogeneity_Q,
      test_name = "Homogeneity Q",
      test_df = result$homogeneity_df,
      p_value = result$homogeneity_p,
      chi2 = result$chi2,
      chi2_p = result$p_value,
      conf.level = conf.level,
      strata_data = result$strata_data,
      n_strata = result$n_total_strata,
      n_informative = result$n_informative,
      stratified = !is.null(strata),
      call = match.call()
    )

  } else if (measure == "irr") {

    if (is.null(time)) {
      stop("'time' must be specified for rate ratio (measure = 'irr')")
    }

    # Call stmh_r function (outcome is the event variable for IRR)
    result <- stmh_r(
      data = data,
      event = outcome,
      exposure = exposure,
      time = time,
      strata = strata,
      conf.level = conf.level
    )

    # Standardise output
    out <- list(
      measure = "Incidence Rate Ratio",
      estimate = result$irr,
      ci = result$ci,
      se_ln = result$se,
      test_stat = result$homogeneity_Q,
      test_name = "Homogeneity Q",
      test_df = result$homogeneity_df,
      p_value = result$homogeneity_p,
      chi2 = result$chi2,
      chi2_p = result$chi2_p,
      conf.level = conf.level,
      strata_data = result$strata,
      stratified = !is.null(strata),
      call = match.call()
    )
  }

  class(out) <- "mh_analysis"
  invisible(out)
}

#' Print Method for mh_analysis Objects
#'
#' Formats and displays the results of a Mantel-Haenszel analysis in a clean,
#' human-readable format.
#'
#' @param x An object of class \code{"mh_analysis"}, typically returned by
#'   \code{\link{mh_analysis}}.
#' @param ... Additional arguments passed to print methods (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @seealso \code{\link{mh_analysis}} for creating mh_analysis objects.
#'
#' @examples
#' # Using the included lung cancer case-control dataset
#' data(lung_cancer_cc)
#'
#' result <- mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case",
#'                       strata = "age_group", measure = "or")
#' print(result)
#'
#' @export
print.mh_analysis <- function(x, ...) {

  cat("\n")
  cat("Mantel-Haenszel", x$measure, "\n")
  cat(rep("=", 40), "\n", sep = "")
  cat("\n")
  cat(x$measure, ":", round(x$estimate, 4), "\n")
  cat(x$conf.level * 100, "% CI: ", round(x$ci[1], 4), " - ", round(x$ci[2], 4), "\n", sep = "")
  cat("SE(ln):", round(x$se_ln, 4), "\n")

  if (!is.na(x$test_stat)) {
    cat("\n")
    cat(x$test_name, ":", round(x$test_stat, 4))
    if (!is.null(x$test_df)) cat(" (df =", x$test_df, ")")
    cat("\n")
    cat("p-value:", format_pvalue(x$p_value), "\n")
  }

  invisible(x)
}
