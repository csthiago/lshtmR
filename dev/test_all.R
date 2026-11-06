# =============================================================================
# lshtmR - test script
#
# Run from the package root (the folder with DESCRIPTION):
#   Rscript dev/test_all.R
# or in RStudio: source("dev/test_all.R")
#
# Each check prints PASS/FAIL. Expected values are computed independently
# (by hand from 2x2 tables / rates, or with stats::mantelhaen.test), not
# copied from the package output.
# =============================================================================

suppressMessages(devtools::load_all(quiet = TRUE))

n_pass <- 0; n_fail <- 0
check <- function(name, ok) {
  ok <- isTRUE(ok)
  cat(if (ok) "PASS " else "FAIL ", name, "\n", sep = "")
  if (ok) n_pass <<- n_pass + 1 else n_fail <<- n_fail + 1
}
near <- function(a, b, tol = 1e-6) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))
# run silently (no printed output, no warnings)
quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressWarnings(expr)))
  out
}
section <- function(x) cat("\n---", x, "---\n")

lc <- lung_cancer_cc
mc <- mortality_cohort

# -----------------------------------------------------------------------------
section("Odds ratio - binary exposure")
# -----------------------------------------------------------------------------
tab <- table(exp = lc$smoking, case = lc$case)
or_crude <- (tab["1", "1"] * tab["0", "0"]) / (tab["1", "0"] * tab["0", "1"])
r <- quiet(mh_analysis(lc, "smoking", "case"))
check("crude OR = ad/bc", near(r$estimate, or_crude))
check("crude OR < CI upper and > CI lower", r$ci[1] < r$estimate && r$estimate < r$ci[2])

mh <- mantelhaen.test(table(lc$smoking, lc$case, lc$age_group), correct = FALSE)
r <- quiet(mh_analysis(lc, "smoking", "case", strata = "age_group"))
# mantelhaen.test uses table order (0,1) x (0,1) -> OR of the same direction
check("stratified MH OR = mantelhaen.test", near(r$estimate, mh$estimate))
check("stratified MH chi2 = mantelhaen.test", near(r$chi2, unname(mh$statistic)))
check("homogeneity test reported (df = strata - 1)", r$test_df == 2 && !is.na(r$test_stat))

r2 <- quiet(mh_analysis(lc, "smoking", "case", strata = c("age_group", "sex")))
mh2 <- mantelhaen.test(table(lc$smoking, lc$case, interaction(lc$age_group, lc$sex)), correct = FALSE)
check("two strata variables OR = mantelhaen.test", near(r2$estimate, mh2$estimate))

r <- quiet(mh_analysis(lc, "smoking", "case", reference = 1))
check("reference = 1 gives 1/OR", near(r$estimate, 1 / or_crude))

d <- lc; d$smk <- factor(ifelse(d$smoking == 1, "yes", "no"))
r <- quiet(mh_analysis(d, "smk", "case"))
check("factor exposure gives same OR", near(r$estimate, or_crude))

r90 <- quiet(mh_analysis(lc, "smoking", "case", conf.level = 0.90))
r95 <- quiet(mh_analysis(lc, "smoking", "case"))
check("90% CI narrower than 95% CI", r90$ci[1] > r95$ci[1] && r90$ci[2] < r95$ci[2])

# -----------------------------------------------------------------------------
section("Rate ratio - binary exposure")
# -----------------------------------------------------------------------------
d1 <- sum(mc$death[mc$treatment == 1]); py1 <- sum(mc$person_years[mc$treatment == 1])
d0 <- sum(mc$death[mc$treatment == 0]); py0 <- sum(mc$person_years[mc$treatment == 0])
irr_crude <- (d1 / py1) / (d0 / py0)
r <- quiet(mh_analysis(mc, "treatment", "death", time = "person_years", measure = "irr"))
check("crude IRR = rate1/rate0", near(r$estimate, irr_crude))

mh_irr <- function(dat, strata) {
  s <- split(dat, interaction(dat[strata], drop = TRUE))
  num <- den <- 0
  for (x in s) {
    a <- sum(x$death[x$treatment == 1]); p1 <- sum(x$person_years[x$treatment == 1])
    b <- sum(x$death[x$treatment == 0]); p0 <- sum(x$person_years[x$treatment == 0])
    num <- num + a * p0 / (p1 + p0); den <- den + b * p1 / (p1 + p0)
  }
  num / den
}
r <- quiet(mh_analysis(mc, "treatment", "death", time = "person_years",
                       strata = "severity", measure = "irr"))
check("stratified MH IRR = hand formula", near(r$estimate, mh_irr(mc, "severity")))

r <- quiet(mh_analysis(mc, "treatment", "death", time = "person_years",
                       strata = c("severity", "age_cat"), measure = "irr"))
check("two strata variables IRR (used to error)", !is.null(r) && near(r$estimate, mh_irr(mc, c("severity", "age_cat"))))

# first stratum has only exposed people -> columns used to be swapped
d <- mc; d$s <- ifelse(d$severity == "Mild", "A", "B")
d <- d[!(d$s == "A" & d$treatment == 0), ]
r <- quiet(mh_analysis(d, "treatment", "death", time = "person_years", strata = "s", measure = "irr"))
sB <- r$strata_data[r$strata_data$stratum == "B", ]
check("stratum counts not swapped (d1 = exposed events)",
      sB$d1 == sum(d$death[d$s == "B" & d$treatment == 1]) &&
      sB$d0 == sum(d$death[d$s == "B" & d$treatment == 0]))
check("stratified IRR correct when a stratum lacks unexposed", near(r$estimate, mh_irr(d, "s")))

d <- mc; d$death <- d$death + 1
r <- quiet(mh_analysis(d, "treatment", "death", time = "person_years", measure = "irr"))
check("event coded 1/2 gives same IRR", near(r$estimate, irr_crude))
d$death <- factor(ifelse(mc$death == 1, "yes", "no"))
r <- quiet(mh_analysis(d, "treatment", "death", time = "person_years", measure = "irr"))
check("factor event gives same IRR", near(r$estimate, irr_crude))

d <- mc; d$person_years[1:50] <- NA
r <- quiet(mh_analysis(d, "treatment", "death", time = "person_years", measure = "irr"))
rc <- quiet(mh_analysis(d[!is.na(d$person_years), ], "treatment", "death", time = "person_years", measure = "irr"))
check("missing person-time rows are dropped", near(r$estimate, rc$estimate))
w <- tryCatch({capture.output(mh_analysis(d, "treatment", "death", time = "person_years", measure = "irr")); NULL},
              warning = function(w) conditionMessage(w))
check("missing person-time gives a warning", grepl("person-time", w))

o <- capture.output(mh_analysis(mc, "treatment", "death", time = "person_years",
                                measure = "irr", conf.level = 0.9))
check("IRR CI label follows conf.level", any(grepl("^90% CI", o)) && !any(grepl("^95% CI", o)))

d <- mc; d$death[d$severity == "Mild"] <- 0
r <- quiet(mh_analysis(d, "treatment", "death", time = "person_years", strata = "severity", measure = "irr"))
check("homogeneity ignores strata without events (no NaN)", !is.nan(r$test_stat) && r$test_df == 1)

# -----------------------------------------------------------------------------
section("Multi-level exposure: each level vs reference")
# -----------------------------------------------------------------------------
r <- quiet(mh_analysis(lc, "age_group", "case", strata = "sex"))
check("type is multilevel", identical(r$type, "multilevel"))
check("reference defaults to first level", r$reference == "40-54")
check("one comparison per non-reference level", nrow(r$comparisons) == 2)
sub <- lc[lc$age_group %in% c("40-54", "70+"), ]; sub$x <- as.integer(sub$age_group == "70+")
rb <- quiet(mh_analysis(sub, "x", "case", strata = "sex"))
check("70+ vs 40-54 = binary analysis on subset", near(r$estimate[["70+"]], rb$estimate))

r <- quiet(mh_analysis(mc, "severity", "death", time = "person_years",
                       strata = "age_cat", measure = "irr", reference = "Moderate"))
sub <- mc[mc$severity %in% c("Moderate", "Severe"), ]; sub$treatment <- as.integer(sub$severity == "Severe")
check("IRR Severe vs Moderate = hand MH formula", near(r$estimate[["Severe"]], mh_irr(sub, "age_cat")))
check("reference argument respected", r$reference == "Moderate" && !("Moderate" %in% r$comparisons$level))

d <- lc; d$ag <- as.integer(d$age_group)
r <- quiet(mh_analysis(d, "ag", "case"))
check("numeric multi-level -> comparisons by default", identical(r$type, "multilevel"))
r <- quiet(mh_analysis(d, "ag", "case", trend = TRUE))
check("trend = TRUE -> score test for trend", identical(r$type, "trend") && !is.null(r$chi2))

# hand score test for trend (unstratified, Stata mhodds formula)
x <- d$ag; y <- d$case; N <- length(y); D <- sum(y)
u <- sum(y * x) - D * sum(x) / N
v <- D * (N - D) * (sum(x^2) - sum(x)^2 / N) / (N * (N - 1))
check("OR trend estimate = exp(u/v)", near(r$estimate, exp(u / v)))
check("OR trend chi2 = u^2/v", near(r$chi2, u^2 / v))

mc2 <- mc; mc2$sev <- as.integer(mc2$severity)
r <- quiet(mh_analysis(mc2, "sev", "death", time = "person_years", measure = "irr", trend = TRUE))
x <- mc2$sev; dd <- mc2$death; py <- mc2$person_years
u <- sum(dd * x) - sum(dd) * sum(py * x) / sum(py)
v <- sum(dd) * (sum(py * x^2) - sum(py * x)^2 / sum(py)) / sum(py)
check("IRR trend estimate = exp(u/v)", near(r$estimate, exp(u / v)))
check("IRR trend chi2 returned (used to be NULL)", near(r$chi2, u^2 / v))

# -----------------------------------------------------------------------------
section("Printing, verbose and summary")
# -----------------------------------------------------------------------------
run <- capture.output(r <- mh_analysis(lc, "age_group", "case", strata = "sex"))
check("print(object) == output while running", identical(run, capture.output(print(r))))
run <- capture.output(r <- mh_analysis(lc, "smoking", "case", verbose = FALSE))
check("verbose = FALSE prints nothing", length(run) == 0)
check("verbose = FALSE still stores full output", length(r$output) > 0)
r <- quiet(mh_analysis(lc, "age_group", "case", strata = "sex", verbose = FALSE))
s <- capture.output(summary(r))
check("summary(multilevel) shows only the table", any(grepl("each level vs. reference", s)) &&
        !any(grepl("2x2 Table", s)))
r <- quiet(mh_analysis(d, "ag", "case", trend = TRUE, verbose = FALSE))
s <- capture.output(summary(r))
check("summary(trend) shows score test chi2/p", any(grepl("Score test for trend", s)) && any(grepl("Chi2", s)))

# -----------------------------------------------------------------------------
section("Missing data reporting and p-value formatting")
# -----------------------------------------------------------------------------
d <- lc; d$smoking[1:50] <- NA
o <- suppressWarnings(capture.output(mh_analysis(d, "smoking", "case")))
check("missing data shown as a real percentage (10%)", any(grepl("50 rows \\(10%\\)", o)))
check("p-value never in scientific notation", format_pvalue(0.0002) == "0.0002")
check("p-value < 0.0001 shown as <0.0001", format_pvalue(1e-6) == "<0.0001")

# -----------------------------------------------------------------------------
section("strate")
# -----------------------------------------------------------------------------
s <- quiet(strate(mc, "death", "person_years", strata = "severity"))
rate_mild <- sum(mc$death[mc$severity == "Mild"]) / sum(mc$person_years[mc$severity == "Mild"])
check("rate = events / person-time", near(s$rate[s$severity == "Mild"], rate_mild))
D <- sum(mc$death[mc$severity == "Mild"])
check("CI = exp(ln(rate) +/- 1.96/sqrt(D))",
      near(s$ci_lower[s$severity == "Mild"], rate_mild * exp(-qnorm(0.975) / sqrt(D))))
d <- mc; d$severity[1:3] <- NA
s <- quiet(strate(d, "death", "person_years", strata = c("severity", "age_cat")))
check("strate with 2 strata + missing strata (used to error)", !is.null(s) && sum(s$events) == sum(d$death[-(1:3)]))

# -----------------------------------------------------------------------------
section("Input errors")
# -----------------------------------------------------------------------------
err <- function(expr) inherits(tryCatch(quiet(expr), error = function(e) e), "error")
check("irr without time -> error", err(mh_analysis(mc, "treatment", "death", measure = "irr")))
check("unknown reference -> error", err(mh_analysis(lc, "age_group", "case", reference = "99")))
check("non-numeric trend -> error", err(mh_analysis(lc, "age_group", "case", trend = TRUE)))
check("unknown variable -> error", err(mh_analysis(lc, "nope", "case")))

# -----------------------------------------------------------------------------
cat("\n=============================\n")
cat("Passed:", n_pass, " Failed:", n_fail, "\n")
cat("=============================\n")
if (n_fail > 0 && !interactive()) quit(status = 1)
