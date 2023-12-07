library(shiny)
library(r3dmol)
library(colourpicker)
library(shinycssloaders)

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
    actionButton(
      inputId = "run_button",
      label = "Run"
    )
  ),
  mainPanel(
    tabsetPanel(
      id="inTabset",
      tabPanel("About",
               tags$h3("Description:"),
               uiOutput("description"),
               tags$h3("How to Use:"),
               uiOutput("how_to"),
               actionButton(
                 inputId = "tutorial_button",
                 label = "View Tutorial"
               ),
               ),
      tabPanel("Tutorial",
               tags$h3("Instructions:"),
               uiOutput("instructions"),
               # tags$h3("Demo:"),
               # uiOutput("demo"),
               # actionButton(
               #   inputId = "demo_button",
               #   label = "Run Demo"
               # ),
               ),
      tabPanel("Visualization",
               r3dmolOutput(outputId = "r3dmol", height = "700px")
               ),
      tabPanel("Results",
               h4("TM-Score:"),
               textOutput(outputId = "tmscore_output"),
               h4("RMSD:"),
               textOutput(outputId = "rmsd_output"),
               ),
      )
    )
  )

server <- function(input, output, session) {
  output$description <- renderUI(
    HTML("<p> This is a Shiny App for the TMscoreAlign R package
          (Zhu, 2023). <code>TMscoreAlign</code> is an R package that
          implements the TM-align program in R <a href=
          'https://doi.org/10.1002/prot.20264'>
          (Zhang & Skolnick, 2005)</a>. TM-align is a method for aligning
          proteins based on TM-score (Template Modelling score), which
          calculates topological similarity between two protein structures
          <a href='https://doi.org/10.1093/nar/gki524'>
          (Zhang & Skolnick, 2004)</a>. TM-score improves upon existing
          means of structural comparison such as RMSD (root-mean-square
          deviation) as it is length-independent and more sensitive to
          global similarities between structures rather than local
          deviations. This package is targeted for structural biologists
          who may use this package to investigate different conformations
          of the same protein, informing a structural basis of protein
          functions. </p>"
    )
  )
  output$how_to <- renderUI(
    HTML("<p> You can run <code>browseVignettes('TMscoreAlign')</code> to get
          more information about using this package. You can learn how
          to use this Shiny app by going over the tutorial:</p>"
         )
    )
  output$instructions <- renderUI(
    HTML("<ol type='1'><li><h5>Upload PDBs:</h5></li>
          <p> First, you nead to read in two PDB files. These will
          contain the protein structures that you will perform
          structural alignment on. Use <b>Upload PDB File 1</b> and
          <b>Upload PDB File 2</b> to upload your PDB files. You can find
          sample data here: <a href=
          'https://github.com/kevqyzhu/TMscoreAlign/blob/main/inst/extdata/1LNIA_decoy1_4.pdb'>
          sample_pdb1</a> and
          <a href=
          'https://github.com/kevqyzhu/TMscoreAlign/blob/main/inst/extdata/1LNIA_decoy2_180.pdb'>
          sample_pdb2</a>. You can get more information about this data by
          running <code>?sample_pdb1</code> and <code>?sample_pdb2</code>.</p>

          <li><h5>Specify Chain IDs:</h5></li>
          <p> PDB files can have multiple chains, so you must also specify for
          each PDB file which chain you want to use in this analysis. Use
          <b>Enter Chain ID for PDB File 1</b> and
          <b>Enter Chain ID for PDB File 2</b> to specify which chain
          you want to use for your first and second PDB files
          respectively.</p>

          <li><h5>Specify Chain Colors:</h5></li>
          <p> The final alignment structure will have two chains,
          one for each PDB you uploaded as inputs. Use
          <b>Select Chain 1 Color</b> and
          <b>Select Chain 2 Color</b> to specify the color you want to
          visualize for your first and second chains respectively.</p>

          <li><h5>Run and Visualize:</h5></li>
          <p> After pressing <b>Run</b>, you will be able to visualize your
          alignments and results, including TM-Score and RMSD.</p>
          </ol>"
         )
    )
  output$demo <- renderUI(
    HTML("<p> Press the button below to run the demo. This will load in the
    aforementioned sample PDB files and perform the alignment."
    )
  )
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

  observeEvent(input$tutorial_button, {
    updateTabsetPanel(session, "inTabset", selected = "Tutorial")
  })

  # observeEvent(input$demo_button, {
  #   demoFilePath <- system.file("extdata", "1LNIA_decoy1_4.pdb",
  #                               package="TMscoreAlign")
  #   if (is.null(input$upload_pdb1)) {
  #     demoFilePath
  #     }
  #   })

  observeEvent(input$run_button, {
    updateTabsetPanel(session, "inTabset", selected = "Visualization")
    chain1 <- isolate(input$chain1_id)
    chain2 <- isolate(input$chain2_id)

    alignment <- get_alignment(uploaded_pdb1(), uploaded_pdb2(),
                               chain1, chain2, method = "alignment",
                               optimize = FALSE
    )
    alignment <- optimize_alignment(alignment, maxit = 900, restart = TRUE)

    outfile <- tempfile(fileext = '.pdb')

    write_pdb(alignment, outputfile = outfile, appended = TRUE,
              uploaded_pdb1(), uploaded_pdb2(), chain1, chain2
    )

    output$r3dmol <- renderR3dmol({
      visualize_alignment_pdb(outfile, chain1 = input$chain1_color,
                              chain2 = input$chain2_color
      )
    })

    # Display additional results in the "Results" tab
    output$tmscore_output <- renderText({
      paste(round(get_tmscore(alignment), 3), "\n")
      })

    output$rmsd_output <- renderText({
      paste(round(get_rmsd(alignment), 3), "\n")
      })
  })
}

shinyApp(ui, server)
