# Purpose: Run TMscoreAlign shiny app
# Author: Kevin Zhu
# Date: 12.10.2023
# Version: 1.0.0
# Bugs and Issues: N/A

#' Launch Shiny App For Package TMscoreAlign
#'
#' A function that launches the shiny app for this package.
#' The shiny app provides an interface giving the user the ability to perform
#' structural alignment to perform protein structure alignment using the
#' TM-score metric. Extra results, such as the actual TM-score, RMSD, and
#' transformation matrix are also provided to the user.
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
#' @export
#' @importFrom shiny runApp
#' @importFrom colourpicker colourInput
runTMscoreAlign <- function() {
  appDir <- system.file("shiny-scripts", package = "TMscoreAlign")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
