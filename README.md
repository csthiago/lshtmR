# lshtmR

You can install using

``` r
devtools::install_github("csthiago/lshtmR")
```

Minimal example

``` r
library(lshtmR)
data("lung_cancer_cc")
data("mortality_cohort")
```

# Odds Ratio

``` r


mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case", measure = "or",
            strata="age_group")

mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case", measure = "or")
```

# Incidence rates

``` r
strate(
  data = mortality_cohort,
  event = "death",
  time = "person_years",
  strata = "age_cat",
  per = 1000,

)
```

# Incidence Rate Ratio

``` r
mh_analysis(mortality_cohort, exposure = "treatment", outcome = "death",
            time = "person_years", measure = "irr")

mh_analysis(mortality_cohort, exposure = "treatment", outcome = "death",
            strata = "age_cat",
            time = "person_years", measure = "irr")
```
