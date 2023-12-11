# Purpose: Run TMscoreAlign shiny app
# Author: Kevin Zhu
# Date: 12.10.2023
# Version: 1.0.0
# Bugs and Issues: N/A

#' Launch Shiny App For Package TMscoreAlign
#'
#' A function that launches the shiny app for this package.
#' The shiny app offers users an interactive interface that allows them to
#' conduct protein structure alignment utilizing the TM-score metric.
#' Supplementary results, including the actual TM-score, RMSD, and
#' transformation matrix, are also provided to the user.
#'
#' @return No return value but a shiny page is opened.
#'
#' @examples
#' \dontrun{
#' runTMscoreAlign()
#' }
#'
#' @references
#' Zhang, Y., and Skolnick, J. (2004). Scoring function for automated assessment
#' of protein structure template quality. \emph{Proteins, Structure, Function,
#' and Bioinformatics}, 57(4), 702–710.
#' \href{https://doi.org/10.1002/prot.20264}{Link}
#'
#' Zhang, Y., and Skolnick, J. (2005). TM-align: a protein structure alignment
#' algorithm based on the TM-score. \emph{Nucleic Acids Res}, 33(7), 2302-2309.
#' \href{https://doi.org/10.1093\%2Fnar\%2Fgki524}{Link}
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom colourpicker colourInput
runTMscoreAlign <- function() {
  appDir <- system.file("shiny-scripts", package = "TMscoreAlign")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
