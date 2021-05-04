
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

Temp. dependent solution:

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
<table>
<thead>
<tr>
<th style="text-align:left;">
Quantity
</th>
<th style="text-align:left;">
Input value
</th>
<th style="text-align:right;">
Model Result
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CP\_100
</td>
<td style="text-align:left;">

-   </td>
    <td style="text-align:right;">
    7.348684
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    CP\_200
    </td>
    <td style="text-align:left;">

    -   </td>
        <td style="text-align:right;">
        16.424254
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        CP\_298.15
        </td>
        <td style="text-align:left;">
        20
        </td>
        <td style="text-align:right;">
        20.000000
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        S\_100
        </td>
        <td style="text-align:left;">

        -   </td>
            <td style="text-align:right;">
            3.174880
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            S\_200
            </td>
            <td style="text-align:left;">

            -   </td>
                <td style="text-align:right;">
                11.497100
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                S\_298.15
                </td>
                <td style="text-align:left;">
                18.82
                </td>
                <td style="text-align:right;">
                18.823450
                </td>
                </tr>
                </tbody>
                </table>

</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
Parameter
</th>
<th style="text-align:right;">
Value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Td
</td>
<td style="text-align:right;">
641.71863
</td>
</tr>
<tr>
<td style="text-align:left;">
a1\_T\_dep\_sol
</td>
<td style="text-align:right;">
531.25000
</td>
</tr>
<tr>
<td style="text-align:left;">
b1\_T\_dep\_sol
</td>
<td style="text-align:right;">
0.37051
</td>
</tr>
<tr>
<td style="text-align:left;">
b\_sol\_coef
</td>
<td style="text-align:right;">
0.00000
</td>
</tr>
<tr>
<td style="text-align:left;">
Number of atoms
</td>
<td style="text-align:right;">
1.00000
</td>
</tr>
</tbody>
</table>
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
<table>
<thead>
<tr>
<th style="text-align:left;">
Quantity
</th>
<th style="text-align:left;">
Input value
</th>
<th style="text-align:right;">
Model Result
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CP\_100
</td>
<td style="text-align:left;">

-   </td>
    <td style="text-align:right;">
    12.95520
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    CP\_200
    </td>
    <td style="text-align:left;">

    -   </td>
        <td style="text-align:right;">
        21.46383
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        CP\_298.15
        </td>
        <td style="text-align:left;">
        24.209
        </td>
        <td style="text-align:right;">
        24.21027
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">
        S\_100
        </td>
        <td style="text-align:left;">

        -   </td>
            <td style="text-align:right;">
            6.89097
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">
            S\_200
            </td>
            <td style="text-align:left;">

            -   </td>
                <td style="text-align:right;">
                19.12283
                </td>
                </tr>
                <tr>
                <td style="text-align:left;">
                S\_298.15
                </td>
                <td style="text-align:left;">
                28.275
                </td>
                <td style="text-align:right;">
                28.29019
                </td>
                </tr>
                </tbody>
                </table>

</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
Parameter
</th>
<th style="text-align:right;">
Value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Td
</td>
<td style="text-align:right;">
402.10000
</td>
</tr>
<tr>
<td style="text-align:left;">
a1\_T\_dep\_sol
</td>
<td style="text-align:right;">
0.00000
</td>
</tr>
<tr>
<td style="text-align:left;">
b1\_T\_dep\_sol
</td>
<td style="text-align:right;">
0.00000
</td>
</tr>
<tr>
<td style="text-align:left;">
b\_sol\_coef
</td>
<td style="text-align:right;">
0.00468
</td>
</tr>
<tr>
<td style="text-align:left;">
Number of atoms
</td>
<td style="text-align:right;">
1.00000
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>
