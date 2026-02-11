#' Simulated Case-Control Study Data for Lung Cancer
#'
#' A simulated dataset from a case-control study examining the association
#' between smoking and lung cancer, stratified by age group and sex.
#' This dataset is designed for demonstrating the \code{\link{mh_or}} function.
#'
#' @format A data frame with 500 rows and 5 variables:
#' \describe{
#'   \item{id}{Unique participant identifier.}
#'   \item{case}{Case status: 1 = lung cancer case, 0 = control.}
#'   \item{smoking}{Smoking exposure: 1 = smoker, 0 = non-smoker.}
#'   \item{age_group}{Age group: "40-54", "55-69", or "70+".}
#'   \item{sex}{Sex: "Male" or "Female".}
#' }
#'
#' @details
#' The data were simulated with a true odds ratio of approximately 3 for the
#' association between smoking and lung cancer. Smoking prevalence varies by
#' age group and sex, and the case probability increases with age.
#'
#' @source Simulated data for educational purposes.
#'
#' @examples
#' # Load the data
#' data(lung_cancer_cc)
#'
#' # View structure
#' str(lung_cancer_cc)
#'
#' # Calculate stratified odds ratio by age group
#' mh_or(lung_cancer_cc, "smoking", "case", "age_group")
#'
#' # Calculate stratified odds ratio by sex
#' mh_or(lung_cancer_cc, "smoking", "case", "sex")
#'
#' # Stratify by both age group and sex
#' mh_or(lung_cancer_cc, "smoking", "case", c("age_group", "sex"))
#'
"lung_cancer_cc"


#' Simulated Cohort Study Data for Mortality
#'
#' A simulated dataset from a cohort study examining the effect of a treatment
#' on mortality, stratified by disease severity and age category.
#' This dataset is designed for demonstrating the \code{\link{stmh_r}} function.
#'
#' @format A data frame with 600 rows and 6 variables:
#' \describe{
#'   \item{id}{Unique participant identifier.}
#'   \item{death}{Death outcome: 1 = died during follow-up, 0 = survived.}
#'   \item{treatment}{Treatment exposure: 1 = treated, 0 = untreated.}
#'   \item{person_years}{Follow-up time in person-years.}
#'   \item{severity}{Disease severity: "Mild", "Moderate", or "Severe".}
#'   \item{age_cat}{Age category: "Young", "Middle", or "Old".}
#' }
#'
#' @details
#' The data were simulated with a true incidence rate ratio of approximately 0.6
#' for the protective effect of treatment on mortality. Mortality rates increase
#' with disease severity and age, and treatment probability decreases with
#' severity.
#'
#' @source Simulated data for educational purposes.
#'
#' @examples
#' # Load the data
#' data(mortality_cohort)
#'
#' # View structure
#' str(mortality_cohort)
#'
#' # Calculate crude rate ratio (no stratification
#' stmh_r(mortality_cohort, "death", "treatment", "person_years")
#'
#' # Calculate stratified rate ratio by severity
#' stmh_r(mortality_cohort, "death", "treatment", "person_years", strata = "severity")
#'
#' # Calculate stratified rate ratio by age category
#' stmh_r(mortality_cohort, "death", "treatment", "person_years", strata = "age_cat")
#'
"mortality_cohort"
