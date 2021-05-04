
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
```

<table class="kable_wrapper">
<caption>
Calphad 2.5 model results
</caption>
<tbody>
<tr>
<td>

| Quantity   | Input value | Model Result |
|:-----------|:------------|-------------:|
| CP\_100    | \-          |     7.348684 |
| CP\_200    | \-          |    16.424254 |
| CP\_298.15 | 20          |    20.000000 |
| S\_100     | \-          |     3.174880 |
| S\_200     | \-          |    11.497100 |
| S\_298.15  | 18.82       |    18.823450 |

</td>
<td>

| Parameter       |     Value |
|:----------------|----------:|
| Td              | 641.71863 |
| a1\_T\_dep\_sol | 531.25000 |
| b1\_T\_dep\_sol |   0.37051 |
| b\_sol\_coef    |   0.00000 |
| Number of atoms |   1.00000 |

</td>
</tr>
</tbody>
</table>

## Example - Pure Al

Temp dependent solution:

``` r
library(Calphad2.5)
Calphad_2.5("Al",  CP298 = 24.209,  S298 = 28.275)
```

<table class="kable_wrapper">
<caption>
Calphad 2.5 model results
</caption>
<tbody>
<tr>
<td>

| Quantity   | Input value | Model Result |
|:-----------|:------------|-------------:|
| CP\_100    | \-          |     12.95520 |
| CP\_200    | \-          |     21.46383 |
| CP\_298.15 | 24.209      |     24.21027 |
| S\_100     | \-          |      6.89097 |
| S\_200     | \-          |     19.12283 |
| S\_298.15  | 28.275      |     28.29019 |

</td>
<td>

| Parameter       |     Value |
|:----------------|----------:|
| Td              | 402.10000 |
| a1\_T\_dep\_sol |   0.00000 |
| b1\_T\_dep\_sol |   0.00000 |
| b\_sol\_coef    |   0.00468 |
| Number of atoms |   1.00000 |

</td>
</tr>
</tbody>
</table>
