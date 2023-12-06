library(shiny)
library(r3dmol)
library(colourpicker)

ui <- pageWithSidebar(
  headerPanel("TMscoreAlign: Protein Structure Alignment using TM-score"),
  sidebarPanel(
    fileInput(
      inputId = "upload_pdb1",
      label = "Upload PDB File 1",
      accept = c(".pdb")
    ),
    fileInput(
      inputId = "upload_pdb2",
      label = "Upload PDB File 2",
      accept = c(".pdb")
    ),
    colourpicker::colourInput(
      inputId = "chain1_color",
      label = "Select Chain 1 Color",
      closeOnClick = TRUE,
      value = "#636efa"  # Blue
    ),
    colourpicker::colourInput(
      inputId = "chain2_color",
      label = "Select Chain 2 Color",
      closeOnClick = TRUE,
      value = "#ff7f0e"  # Orange
    ),
    textInput(
      inputId = "chain1_id",
      label = "Enter Chain ID for PDB File 1",
      value = "A"
    ),
    textInput(
      inputId = "chain2_id",
      label = "Enter Chain ID for PDB File 2",
      value = "A"
    ),
    actionButton(
      inputId = "run_button",
      label = "Run"
    )
  ),
  mainPanel(
    r3dmolOutput(outputId = "r3dmol", height = "700px"),
  )
)

server <- function(input, output, session) {
  uploaded_pdb1 <- reactive({
    req(input$upload_pdb1)
    inFile <- input$upload_pdb1
    return(inFile$datapath)
  })

  uploaded_pdb2 <- reactive({
    req(input$upload_pdb2)
    inFile <- input$upload_pdb2
    return(inFile$datapath)
  })

  observeEvent(input$run_button, {
    chain1 <- isolate(input$chain1_id)
    chain2 <- isolate(input$chain2_id)

    alignment <- get_alignment(uploaded_pdb1(), uploaded_pdb2(),
                               chain1, chain2, method = "alignment",
                               optimize = TRUE
    )
    alignment <- optimize_alignment(alignment, maxit = 400, restart = TRUE)

    isolate(print(get_tmscore(alignment)))

    outfile <- tempfile(fileext = '.pdb')

    write_pdb(alignment, outputfile = outfile, appended = TRUE,
              uploaded_pdb1(), uploaded_pdb2(), chain1, chain2
    )

    chain1_color <- isolate(input$chain1_color)
    chain2_color <- isolate(input$chain2_color)

    output$r3dmol <- renderR3dmol({
      visualize_alignment_pdb(outfile, chain1 = chain1_color,
                              chain2 = chain2_color
      )
    })

  })

}

shinyApp(ui, server)
