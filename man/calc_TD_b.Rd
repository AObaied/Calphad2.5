% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b_solution.R
\name{calc_TD_b}
\alias{calc_TD_b}
\title{Obtain Td and b for b_solution}
\usage{
calc_TD_b(L1, L2)
}
\arguments{
\item{L1}{parameter}

\item{L2}{parameter}
}
\value{
Solve the 2 systems of equations to obtain Td and b if there are no solution obtained,
then the element with the problem will be saved into a Csv file named (elements_with_problems.Csv)
}
\description{
Obtain Td and b for b_solution
}
