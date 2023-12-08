#' Launch Shiny App For Package TMscoreAlign
#'
#' A function that launches the shiny app for this package.
#' TODO: finish documentation for this
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' runTMscoreAlign()
#' }
#'
#' @author Kevin Zhu, \email{kevqyzhu@gmail.com} # is this needed?
#'
#' @references
#' TODO: finish documentation for this
#'
#' @export
#' @import xtable
#' @importFrom shiny runApp
#' @importFrom colourpicker colourInput
runTMscoreAlign <- function() {
  appDir <- system.file("shiny-scripts", package = "TMscoreAlign")
  shiny::runApp(appDir, display.mode = "normal")
  # shiny::runApp(appDir, display.mode = "showcase")
  return()
}
