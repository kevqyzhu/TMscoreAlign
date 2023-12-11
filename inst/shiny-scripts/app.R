library(shiny)
library(r3dmol)
library(colourpicker)

ui <- pageWithSidebar(
  headerPanel("TMscoreAlign: Protein Structure Alignment using TM-score"),
  sidebarPanel(
    fileInput(
      inputId = "upload_pdb1",
      placeholder = "sample_pdb1",
      label = "Upload PDB File 1",
      accept = c(".pdb")
      ),
    fileInput(
      inputId = "upload_pdb2",
      placeholder = "sample_pdb2",
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
      tabPanel("Visualization",
               r3dmol::r3dmolOutput(outputId = "r3dmol", height = "700px")
               ),
      tabPanel("Results",
               h4("TM-Score:"),
               textOutput(outputId = "tmscore_output"),
               h4("RMSD:"),
               textOutput(outputId = "rmsd_output"),
               h4("Transformation Matrix:"),
               uiOutput('matrix'),
               h4("Interpreting Results:"),
               uiOutput("explain_results"),
               ),
      tabPanel("Tutorial",
               tags$h3("Instructions:"),
               uiOutput("instructions"),
               ),
      )
    )
  )

server <- function(input, output, session) {
  output$description <- renderUI(
    HTML("<p> This is a Shiny App for the TMscoreAlign R package
          (Zhu, 2023). <code>TMscoreAlign</code> is an R package that
          implements the TM-align program in R <a href=
          'https://doi.org/10.1093/nar/gki524'>
          (Zhang & Skolnick, 2005)</a>. TM-align is a method for aligning
          proteins based on TM-score (Template Modelling score), which
          calculates topological similarity between two protein structures
          <a href='https://doi.org/10.1002/prot.20264'>
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
          <p> First, the user needs to read in two PDB files. These files will
          contain the protein structures for which structural alignment will be
          performed. The user should utilize <b>Upload PDB File 1</b> and
          <b>Upload PDB File 2</b> to upload the PDB files. Example data has
          already been provided:
          <a href=
          'https://github.com/kevqyzhu/TMscoreAlign/blob/main/inst/extdata/1LNIA_decoy1_4.pdb'>
          sample_pdb1</a> and
          <a href=
          'https://github.com/kevqyzhu/TMscoreAlign/blob/main/inst/extdata/1LNIA_decoy2_180.pdb'>
          sample_pdb2</a>. Users can obtain more information about this data by
          running <code>?sample_pdb1</code> and <code>?sample_pdb2</code>.</p>

          <li><h5>Specify Chain IDs:</h5></li>
          <p> Users need to consider that PDB files might contain multiple
          chains. Consequently, for each PDB file, it is necessary to specify
          the chain intended for use in the analysis. Utilize <b>Enter Chain ID
          for PDB File 1</b> and <b>Enter Chain ID for PDB File 2</b> to
          designate the desired chain for the first and second PDB files,
          respectively.</p>

          <li><h5>Specify Chain Colors:</h5></li>
          <p> The resulting alignment structure will consist of two chains,
          each corresponding to the PDB file uploaded as input. Utilize
          <b>Select Chain 1 Color</b> and <b>Select Chain 2 Color</b> to
          determine the color for visualizing the first and second chains,
          respectively.</p>

          <li><h5>Run and Visualize:</h5></li>
          <p> Upon initiating the <b>Run</b> process, users will gain the
          ability to visualize their alignments and results, including
          TM-Score, RMSD, and the transformation matrix.</p>
          </ol>"
         )
    )

  default_output <- renderText({paste("----", "\n")})
  default_matrix <- renderTable({
    M <- matrix(NA, nrow = 3, ncol = 4)
    rownames(M) <- c('1','2','3')
    colnames(M) <- c('Rotation 1','Rotation 2','Rotation 3', 'Translation')
    return(M)
    }, rownames = TRUE)

  output$tmscore_output <- default_output
  output$rmsd_output <- default_output
  output$matrix <- default_matrix

  output$explain_results <- renderUI(
    HTML("<p> <b>TM-score</b> (Template Modelling score) calculates topological
         similarity between two protein structures. It improves upon existing
         means of structural comparison as it is length-independent and more
         sensitive to global similarities between structures rather than local
         deviations. According to <a href='https://doi.org/10.1002/prot.20264'>
         Zhang & Skolnick (2004)</a>: </p>

         <p style='margin-left: 40px'> 0.0 < TM-score < 0.17 indicates that the
         proteins have random structural similarity <br>
         0.5 < TM-score < 1.00 indicates that the proteins are likely in the
         same structural group</p>

         <p> <b>RMSD</b> (root-mean-square deviation) is the most commonly
         used measure to assess structural similarity between two proteins. For
         identical structures, it is 0, with its value increasing as the two
         structures become more different. </p>

         <p> The <b>Transformation Matrix</b> represents the mathematical
         transformation (a combination of rotation and translation) applied to
         one protein to align it with the other structure. The translation and
         rotation parameters are calculated via optimizing for maximum
         TM-score. </p>"
         )
    )
  uploaded_pdb1 <- reactive({
    if (is.null(input$upload_pdb1)) {
      return (system.file("extdata", "1LNIA_decoy1_4.pdb",
                          package="TMscoreAlign"
                          )
              )
      } else {
        req(input$upload_pdb1)
        inFile <- input$upload_pdb1
        return(inFile$datapath)
        }
    })

  uploaded_pdb2 <- reactive({
    if (is.null(input$upload_pdb2)) {
      return (system.file("extdata", "1LNIA_decoy2_180.pdb",
                          package="TMscoreAlign"
                          )
              )
      } else {
        req(input$upload_pdb2)
        inFile <- input$upload_pdb2
        return(inFile$datapath)
        }
    })

  observeEvent(input$tutorial_button, {
    updateTabsetPanel(session, "inTabset", selected = "Tutorial")
    })

  observeEvent(input$run_button, {
    chain1 <- isolate(input$chain1_id)
    chain2 <- isolate(input$chain2_id)

    tryCatch({
      alignment <- get_alignment(uploaded_pdb1(), uploaded_pdb2(),
                                 chain1, chain2, method = "alignment",
                                 optimize = FALSE
                                 )
      alignment <- optimize_alignment(alignment, maxit = 900, restart = TRUE)

      outfile <- tempfile(fileext = '.pdb')

      write_pdb(alignment, outputfile = outfile, appended = TRUE,
                uploaded_pdb1(), uploaded_pdb2(), chain1, chain2
                )

      updateTabsetPanel(session, "inTabset", selected = "Visualization")

      output$r3dmol <- r3dmol::renderR3dmol({
        visualize_alignment_pdb(outfile,
                                chain1 = input$chain1_color,
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

      output$matrix <- renderTable({
        values <- alignment$values
        M <- get_matrix(values)[-4,]
        rownames(M) <- c('1','2','3')
        colnames(M) <- c('Rotation 1','Rotation 2','Rotation 3', 'Translation')
        M
        }, rownames = TRUE)

      }, error = function(e) {
        showNotification((paste(e)), type="error")
        output$r3dmol <- NULL
        output$tmscore_output <- default_output
        output$rmsd_output <- default_output
        output$matrix <- default_matrix
        }
      )
  })
}

shinyApp(ui, server)
