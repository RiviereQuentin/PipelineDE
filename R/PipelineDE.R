PipelineDE <- function(...) {

  ui <- fluidPage(
    titlePanel("Differential expression analysis in R"),
    sidebarLayout(
      sidebarPanel(
        div("blabla", style = "background-color: white; font-size: 16px; color: white"),
        div("Hi Epiplant!", style = "background-color: white; font-size: 20px; color: orange; font-weight: bold"),
        div("blabla", style = "background-color: white; font-size: 16px; color: white"),
        div("To run the differential expression analysis with egdgeR and DESeq2, you will need:", style = "background-color: white; font-size: 20px"),
        div("   - the bam files associated with your experiment,", style = "background-color: white; font-size: 20px; font-weight: bold; padding:0 0 0 20px"),
        div("   - the gtf file representing your reference genome.", style = "background-color: white; font-size: 20px; font-weight: bold; padding:0 0 0 20px"),
        div("blabla", style = "background-color: white; font-size: 16px; color: white"),
        div("The bam files contain the (ordered) aligned sequences obtained with each sample.", style = "background-color: white; font-size: 16px"),
        div("Important:", style = "background-color: white; font-size: 16px; color: orange"),
        div("   - the bam files of the control condition must contain \'ctrl\' or \'control\' (case insensitive) in their names.", style = "background-color: white; font-size: 16px; padding:0 0 0 20px"),
        div("   - if you have paired-end read data, the mates need to be aligned separately. The names of the bam files must contain \'_R1_\' and \'_R2_\' to refer to the first and second mates, respectively.", style = "background-color: white; font-size: 16px; padding:0 0 0 20px"),
        div("blabla", style = "background-color: white; font-size: 16px; color: white"),
        div("The gtf file allows to locate the transcripts and their exons. Use exactly the same gtf file that was used to align the sequences and generate the bam files.", style = "background-color: white; font-size: 16px"),
        div("blabla", style = "background-color: white; font-size: 16px; color: white"),
        br(),
        div("1. Place all the bam files and the gtf file in a folder and upload this folder:", style = "background-color: white; font-size: 20px; font-weight: bold; padding:0 0 0 20px"),
        br(),
        shinyDirButton(id = "directory", label = "Upload...", title = "Data directory"),
        br(),
        br(),
        div("2. Determine the false-discovery rate and log2 fold-change thresholds to identify the up- and down-regulated genes:", style = "background-color: white; font-size: 20px; font-weight: bold; padding:0 0 0 20px"),
        br(),
        numericInput("fdr", "FDR threshold:", 0.05, min = 0, max = 1),
        numericInput("log2fc", "Log2FC threshold (in absolute value):", 1, min = 0, max = 5),
        div("3. Set whether you have paired-end read data or not:", style = "background-color: white; font-size: 20px; font-weight: bold; padding:0 0 0 20px"),
        br(),
        checkboxInput('isPairedEnd', label = "Paired-end reads?", value = TRUE),
        br(),
        div("We are ready! Let\'s go!", style = "background-color: white; font-size: 20px; color: orange"),
        br(),
        actionButton("run", "Run DE analysis!")
      ),
    mainPanel(
      div("Download the results once the output table is printed below", style = "background-color: white; font-size: 20px; font-weight: bold; color: orange"),
      div("blabla", style = "background-color: white; font-size: 16px; color: white"),
      downloadButton("download", "Download results"),
      uiOutput("description"),
      dataTableOutput("DEres")
    )
  )
  )

  server <- function(input, output, session) {
    root <- getwd()
    root <- unlist(strsplit(root, "/"))[c(1,2)]
    root <- paste0(root[1], "/", root[2])
    shinyDirChoose(input, "directory", roots=c(wd=root))
    directory <- eventReactive(input$run, {input$directory})
    fdr <- eventReactive(input$run, {input$fdr})
    log2fc <- eventReactive(input$run, {input$log2fc})
    isPairedEnd <- eventReactive(input$run, {input$isPairedEnd})

    DEres <- eventReactive(input$run, {
      DEinr(DEdirectory = file.path(root, paste0(unlist(directory()$path), collapse = "/")),
                   isPairedEnd = isPairedEnd(), pval.threshold = fdr(), log2FC.threshold = log2fc())
      })

    output$DEres <- renderDataTable({
      DEres()
    })

    output$download <- downloadHandler(filename = "DEG_results.tsv",
                                       content = function(file){
                                         write.table(DEres(),
                                                     file,
                                                     sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE
                                                     )
                                       })

    output$description <- renderUI({
    DE <- DEres()
    tagList(br(),
             div("Description of the output table:", style = "background-color: white; font-size: 16px"),
             div("- The first columns give the names of the differentially expressed genes and the expression values (in TPM) in each sample", style = "background-color: white; font-size: 16px; padding: 0 0 0 20px"),
             div("- The four last columns provide the log2 fold-change (Log2FC) and false-discovery rates (FDR) obtained with edgeR and DESeq2", style = "background-color: white; font-size: 16px; padding: 0 0 0 20px"),
             br(),
             div("Remark:", style = "background-color: white; font-size: 16px; font-weight: bold"),
             div("With edgeR, Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests were applied; with DESeq2, Negative Binomial GLM fitting and Wald significance tests were performed.", style = "background-color: white; font-size: 16px"),
            br()
             )

    })
  }

  shinyApp(ui, server, options = list(launch.browser = TRUE),...)
}
