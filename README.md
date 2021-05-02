
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Calphad2.5

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/AObaied/Calphad2.5.svg?branch=master)](https://travis-ci.com/AObaied/Calphad2.5)
<!-- badges: end -->

The goal of Calphad2.5 is to Genarate Cp and S desriptions below room
temperature, using only two input parameters (Cp298, S298)

## Installation

You can install the released version of Calphad2.5 from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Calphad2.5")
```

## Example - Pure Si

Temp dependent solution:

``` r
library(Calphad2.5)
Calphad_2.5("Si", CP298 = 20,   S298 = 18.82)
#>                 Element             CP100_model             CP200_model 
#>                    "Si"      "7.34868397325056"      "16.4242540835715" 
#>             CP298_model             CP298_Input              S100_model 
#>      "19.9999999661793"                    "20"               "3.17488" 
#>              S200_model              S298_model              S298_Input 
#>               "11.4971"              "18.82345"                 "18.82" 
#>                      Td            a1_T_dep_sol            b1_T_dep_sol 
#>      "641.718628386998"                "531.25"     "0.370513595126608" 
#>                     NOA                  S-diff                    a_S1 
#>                     "1"          "2.0146299925"                     "0" 
#>                    b_S1                    c_S1                    d_S1 
#>  "7.29932299502528e-07"   "-9.121418869777e-07"  "1.31804839681616e-05" 
#>                    e_S1                    f_S1                    g_S1 
#> "-4.58764191863264e-08"   "7.0982874224285e-10" "-8.55373783396784e-12" 
#>                    a_S2                    b_S2                    c_S2 
#>     "-53.8253728702997"     "0.684898243599751"  "-0.00293030703966731" 
#>                    d_S2                    e_S2                    f_S2 
#>  "6.31673816821527e-06" "-5.47783983112289e-09"      "1863.32441236062" 
#>                    g_S2 
#>     "-24148.9542846084"
```

## Example - Pure Al

Temp dependent solution:

``` r
library(Calphad2.5)
Calphad_2.5("Al", CP298 = 24.209,   S298 = 28.275)
#>                 Element             CP100_model             CP200_model 
#>                    "Al"      "12.9552049682052"      "21.4638289939939" 
#>             CP298_model             CP298_Input              S100_model 
#>      "24.2102679753391"                "24.209"               "6.89097" 
#>              S200_model              S298_model              S298_Input 
#>              "19.12283"              "28.29019"                "28.275" 
#>                      Td            a1_T_dep_sol            b1_T_dep_sol 
#>                 "402.1"                "531.25"     "0.370513595126608" 
#>              b_sol_coef                     NOA                  S-diff 
#>   "0.00468156404862815"                     "1"     "-11.5624408966875" 
#>                    a_S1                    b_S1                    c_S1 
#>                     "0"  "6.10338468887655e-06" "-8.65978636340368e-06" 
#>                    d_S1                    e_S1                    f_S1 
#>  "3.22049618964679e-05" "-2.25021948965693e-07"   "9.1769878865714e-09" 
#>                    g_S1                    a_S2                    b_S2 
#> "-1.33722632643667e-10"     "-30.1127795226621"      "0.66334222900953" 
#>                    c_S2                    d_S2                    e_S2 
#>  "-0.00357194238978573"  "9.25860273960745e-06" "-9.34627783268443e-09" 
#>                    f_S2                    g_S2 
#>      "327.801832098728"      "3830.35805582564"
```
