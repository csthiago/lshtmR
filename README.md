# lshtmR

You can install using

```r
devtools::install_github("csthiago/lshtmR")
```

Minimal example
```r
library(lshtmR)
data("lung_cancer_cc")
data("mortality_cohort")

# Odds Ratio

mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case", measure = "or",
            strata="age_group")

mh_analysis(lung_cancer_cc, exposure = "smoking", outcome = "case", measure = "or")
            
# Incidence Rate Ratio
mh_analysis(mortality_cohort, exposure = "treatment", outcome = "death",
            time = "person_years", measure = "irr")

mh_analysis(mortality_cohort, exposure = "treatment", outcome = "death",
            strata = "age_cat",
            time = "person_years", measure = "irr")
```

