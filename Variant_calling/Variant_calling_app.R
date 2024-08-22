install.packages("BiocManager",ask = FALSE)
install.packages("shiny",ask = FALSE)
install.packages("shinyFiles",ask = FALSE)
BiocManager::install("AllelicImbalance",ask = FALSE)
BiocManager::install("Biostrings",ask = FALSE)


library(shiny) ## v1.7.5.1
library(AllelicImbalance) ## v1.36.0
library(shinyFiles) ## v0.9.3
library(Biostrings) ## v2.66.0

# Define UI for data upload app ----
ui <- fluidPage(
  navbarPage("Characterization of Variants",
             tags$head(
               tags$style(HTML("
                    .navbar .navbar-nav {float: right;
                                         font-size: 25px;
                                         background-color: #7B8EDD ;
                                         height:100px;
                                         display: flex;
                                        align-items: center;}
                    .navbar.navbar-default.navbar-static-top{
                                        background-color: #7B8EDD ;
                                        }
                    .navbar .navbar-header {float: left;
                                            align-items: center
                                            }
                    .navbar-default .navbar-brand { color:black;
                                                    font-size:45px;
                                                    font-weight:bold;
                                                    background-color: #7B8EDD ;
                                                    display: flex;
                                                    align-items: center;
                                                    vertical-align: center;}
                "))
             ),
             
             tabPanel("Tool",     
                      tags$head(
                        tags$style(HTML("
          body {
            background-color: #D1DAFE 
          }
          .title-container h1{
            font-size: 50px; /* Adjust the font size as needed */
            font-weight: bold; /* Optionally adjust font weight */
          }
        "))
                      ),
                      p("This tool was created to perform variant calling on human samples. 
         It is based on the AllelicImbalance package available in R with Bioconductor. 
         At the moment, it is only available for 3 genes : TP53, IDH1 and IDH2.
         For TP53, this tool only looks into the 166 most commons mutations referenced in TP53 COSMIC database (V99).
         For IDH1, the mutations taken into account are R132H, R132C, R132G, R132L and R132S and for IDH2 : R140Q, R172G, R172K, R172M, R172S, R172T.
         Before any new request, please reload the page.
         ",style = "color:black;font-size:15px;"),
                      p("To use the tool, the first step is to select a folder that contains at least 2 BAM files (file.bam). Then, the gene of interest has to be specified.
         Finally, when the search button is clicked, the result can take up to several minutes to appear depending on the number of files in the folder.",
                        style = "color:black;font-size:15px;"),
                      p("Once the result table is displayed, you have the possibility to download it (Download Results Table) and/or to download the raw results
         to explore all SNPs found on the gene of interest (Download Raw Results)."
                        ,style = "color:black;font-size:15px;"),
                      tags$a(href = "https://github.com/alizeeln/Stage_M2", "Find the code of the app here, GitHub Repository"),
                      # Sidebar layout with input and output definitions ----
                      sidebarLayout(
                        
                        # Sidebar panel for inputs ----
                        sidebarPanel(
                          
                          # Input: Select a file ----
                          shinyDirButton(id = 'dir', 
                                         label = 'Folder select', 
                                         title = 'Please select the folder containing your BAM'),
                          
                          selectInput(inputId = "species",
                                      label = "Species", 
                                      choices = c("Human" = "human")),
                          
                          selectInput(inputId = "mutations",
                                      label = "Mutations",
                                      choices = c("TP53" = "tp53", "IDH1" = "idh1", "IDH2" = "idh2")),
                          
                          
                          actionButton("search", "Search"),
                          # Horizontal line ----
                          tags$hr(),
                          
                        ),
                        
                        # Main panel for displaying outputs ----
                        mainPanel(
                          actionButton("download_raw", "Download Raw Results"),
                          actionButton("download_table", "Download Results Table"),
                          # Output: Data file ----
                          tableOutput("contents"),
                          
                          tabPanel("logo",
                                   tags$div(align = "left",
                                            style = "position: fixed;
                                                bottom: 0;
                                                width: 100%;
                                                padding: 10px;
                                                background-color: #D1DAFE;",
                                            img(src = "ircm.jpeg", height = "50px"),
                                            img(src = "icm.png", height = "50px"),
                                            img(src = "univ_montp.png", height = "50px"),
                                            img(src = "inserm_logo.png", height = "50px")
                                   )
                          )
                        )
                      )
             ),
             
             tabPanel("Description",
                      p("Hello,",style = "color:black;font-size:20px;font-weight:bold;"),
                      p("This tool was developped by Alizée LANON during her internship (March to July 2024) in Laurent Le CAM's Team (IRCM)"),
                      p("To implement this application, I developped a parser from the output of the AllelicImbalance package in R
            and compared it to a .csv file extracted from the Cosmic database. There is one .csv file per gene at the moment.
             As the output of the AllelicImbalance package is a list of genomic coordinates as illustrated in Figure 5 (one element of the list) on which there is a SNP, the first 
             check that is made is to evaluate if the genomic coordinate of a given element is found in the Cosmic database."
                      ),
                      p(" If a given SNP's genomic coordinate is not found in the 
           cosmic file as illustrated in figure 6 (An example of a cosmic mutation on IDH1 gene), then, it is not processed further as a possible mutation.
           In case of a positive match, it is indicated as a possible mutation which initiates the next step of the analysis.  Then, for each element of the list 
           of possible mutation, we parse each sample to check if it is mutated. It is considered mutated if there are 3 or more reads having the mutated allele.
           However, all variants are listed in the Raw Results output."
                      ),
                      p("A percentage will be calculated between the number of reads that have the mutant nucleotide and the total
           of reads at this position for the given sample. Before calculating a ratio or displaying the sample as mutated, a last verification 
           is made. The algorithm checks that the nucleotide(s) implicated in the mutation observed are the same than in the reference database 
           CDS Mutation column."
                      ),
                      p("The tool takes into consideration the strand on which a gene is located."
                      ),
                      img(src = "exemple_cosmic_package.png", height = "200px"),
                      p("The tool was validated on more than 30 samples originating from 3 differents cohorts.",
                      ),
                      tabPanel("logo",
                               tags$div(align = "center",
                                        style = "position: fixed;
                                                bottom: 0;
                                                width: 100%;
                                                padding: 10px;
                                                background-color: #D1DAFE;",
                                        img(src = "ircm.jpeg", height = "50px"),
                                        img(src = "icm.png", height = "50px"),
                                        img(src = "univ_montp.png", height = "50px"),
                                        img(src = "inserm_logo.png", height = "50px")
                               )
                      )
             )
  )
)


# Define server logic to read selected file ----
server <- function(input, output) {
  print("In server")
  options(shiny.maxRequestSize = 5*1024^3)
  
  notification_id <- reactiveVal(NULL)
  
  shinyFiles::shinyDirChoose(
    input = input,
    id = 'dir',
    defaultPath = "",
    roots = c(home = if (Sys.info()['sysname'] == "Windows") "C:/" else "~"),
    filetypes = c('', 'bam') ### or 'bigWig', "tsv", "csv", "bw"
  )
  
  global <- reactiveValues(datapath = getwd())
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    global$datapath
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath(c(home = if (Sys.info()['sysname'] == "Windows") "C:/" else "~"))
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                 
                 system(command= paste0("ls ",global$datapath))
                 
               })
  observeEvent(input$search,{
    print("In ObserveEvent")
    output$contents <- renderTable({
      req(input$dir)
      id <- showNotification("Processing...", duration = NULL)
      notification_id(id)
      
      observeEvent(input$dir,{
        print("Before conditions")
        if (input$species == "human" & input$mutations == "tp53") {
          searchArea <- GRanges(seqnames = c("17"), ranges = IRanges(7666402,7689550))
          strand <- "reverse"
        }
        else if (input$species == "human" & input$mutations == "idh1") {
          searchArea <- GRanges(seqnames = c("2"), ranges = IRanges(208234227,208257143))
          strand <- "reverse"
        }
        else if (input$species == "human" & input$mutations == "idh2") {
          searchArea <- GRanges(seqnames = c("15"), ranges = IRanges(90081045,90104468))
          strand <- "reverse"
        }
        print(searchArea)
        print("Before bams")
        directory <- paste(unlist(unique(c(c(home = if (Sys.info()['sysname'] == "Windows") "C:/" else "~"),input$dir$path))),sep='/')
        directory <- directory[directory != ""]
        directory <- gsub('"', '', directory)
        directory <- paste(directory, collapse = "/")
        bams <- directory
        
        ## create .bai
        print("Before reads")
        reads <- impBamGAL(bams, searchArea, verbose = FALSE)
        ## look at heterozygous positions
        print("Before heterozygotePositions")
        heterozygotePositions <- scanForHeterozygotes(reads, verbose = FALSE)
        ## create table of heterozygous positions
        print("Before countList")
        countList <- getAlleleCounts(reads, heterozygotePositions, verbose = FALSE)
        print(countList)
        
        print("Before Cosmic Importation")
        if (input$species == "human" & input$mutations == "tp53") {
          cosmic <- read.csv("Cosmic_ref/Gene_mutations_p53_COSMIC_2_1.csv", header = TRUE, sep = ";")
          cosmic <- cosmic
          strand <- "reverse"
          print("tp53")
        }
        else if (input$species == "human" & input$mutations == "idh1") {
          cosmic <- read.csv("Cosmic_ref/Gene_mutations_IDH1_COSMIC.csv", header = TRUE, sep = ";")
          strand <- "reverse"
          print("idh1")
        }
        else if (input$species == "human" & input$mutations == "idh2") {
          cosmic <- read.csv("Cosmic_ref/Gene_mutations_IDH2_COSMIC.csv", header = TRUE, sep = ";")
          strand <- "reverse"
          print("idh2")
        }
        
        filtered_countList <- list()
        sample_names <- list()
        
        result_table <- data.frame(
          Sample_Name = character(),
          Profile = character(),
          Mutation = character(),
          Percentage = numeric(),
          Number_of_wt_reads = numeric(),
          Number_of_mut_reads = numeric(),
          stringsAsFactors = FALSE
        )
        
        print("Before loop over countList results")
        for (e in names(countList)) {
          # Check if the name is present in cosmic$Coordinate
          if (e %in% cosmic$Coordinate) {
            # If the name is present, add the corresponding matrix to the filtered countList
            filtered_countList[[e]] <- countList[[e]]
            # Import cosmic mutations
            results <- subset(cosmic, Coordinate == e)
            # Read the R output of Allelic Imbalance
            for (i in seq_along(filtered_countList)) {
              # Extract the matrix from each element
              matrix_data <- filtered_countList[[i]]
              
              zero_counts <- numeric()
              
              zero_counts_df <- data.frame(column_name = character(), zero_count = numeric(), stringsAsFactors = FALSE)
              
              sample_names <- rownames(matrix_data)
              
              if (strand == "reverse"){
                nuc <- colnames(matrix_data)
                nuc2 <- paste(nuc, collapse = "")
                dna_string <- rev(reverseComplement(DNAString(nuc2)))
                colname <- as.character(dna_string)
                col <- strsplit(colname, split = "")
                colnames(matrix_data) <- col[[1]]
              }
              else if (strand == "forward"){
                colnames(matrix_data) <- colnames(matrix_data)
              }
              for (o in 1:nrow(matrix_data)) {
                # Iterate over columns
                col_hete_homo <- list()
                for (p in 1:ncol(matrix_data)) {
                  if (matrix_data[o,p] > 2){
                    col_hete_homo <- append(col_hete_homo, colnames(matrix_data)[p])
                  }
                }
                ## if mutation is homozygous for this sample
                if (length(col_hete_homo)  == 1 ) {
                  nuc_case <- col_hete_homo[1]
                  # Si des lignes correspondent à la condition, extrayez les valeurs de la colonne "CDS"
                  if (nrow(results) > 0) {
                    for (q in 1:nrow(results)){
                      if (results$Coordinate[q] == names(filtered_countList[i])){
                        CDS_mut <- results$CDS.Mutation[q]
                        nuc_mut_before <- sub(".*([ACGT])>[^>]*$", "\\1", CDS_mut)
                        nuc_mut_after <- sub(".*\\b([^>]*)>", "\\1", CDS_mut)
                        if (nuc_mut_after == nuc_case){
                          new_row <- data.frame(
                            Sample_Name = rownames(matrix_data)[o],
                            Profile = "Mutated",
                            Mutation = results$AA.Mutation[q],
                            Percentage = "100%",
                            Number_of_wt_reads = "No",
                            Number_of_mut_reads = matrix_data[o,nuc_mut_after]
                          )
                          result_table <- rbind(result_table , new_row)
                          print("samples mutated")
                        }
                      }
                    }
                  }
                }
                ## if mutation is heterozygous
                else if (length(col_hete_homo)  == 2 ){
                  #nuc_case <- col_hete_homo[1]
                  # Si des lignes correspondent à la condition, extrayez les valeurs de la colonne "CDS"
                  if (nrow(results) > 0) {
                    for (q in 1:nrow(results)){
                      if (results$Coordinate[q] == names(filtered_countList[i])){
                        CDS_mut <- results$CDS.Mutation[q]
                        nuc_mut_before <- sub(".*([ACGT])>[^>]*$", "\\1", CDS_mut)
                        nuc_mut_after <- sub(".*\\b([^>]*)>", "\\1", CDS_mut)
                        count_ref_col <- matrix_data[o,nuc_mut_before]
                        if (nuc_mut_after %in% col_hete_homo && nuc_mut_before %in% col_hete_homo){
                          new_row <- data.frame(
                            Sample_Name = rownames(matrix_data)[o],
                            Profile = "Mutated",
                            Mutation = results$AA.Mutation[q],
                            Percentage = round((matrix_data[o,nuc_mut_after]/(count_ref_col + matrix_data[o,nuc_mut_after]))*100,4),
                            Number_of_wt_reads = count_ref_col,
                            Number_of_mut_reads = matrix_data[o,nuc_mut_after]
                          )
                          result_table <- rbind(result_table , new_row)
                          print("samples mutated")
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        if (length(sample_names) == 0){
          print("No samples have the chosen mutation")
        } else {
          for (s in sample_names){
            if (!s %in% result_table$Sample_Name){
              new_row <- data.frame(
                Sample_Name = s,
                Profile = "No mutations detected amongst the one chosen (see description)",
                Mutation = "No",
                Percentage = "No",
                Number_of_wt_reads = "No",
                Number_of_mut_reads = "No"
              )
              result_table <- rbind(result_table , new_row)
              print("No samples have mutations on IDH1 gene")
            }
          }
        }
        result_table <- unique(result_table)
        if (length(result_table) == 0){
          new_row <- data.frame(
          Sample_Name = "None",
          Profile = "No mutations detected on the samples",
          Mutation = "No",
          Percentage = "No",
          Number_of_wt_reads = "No",
          Number_of_mut_reads = "No"
          )
          result_table <- rbind(result_table , new_row)
        }
        
          output$contents <- renderTable({
            removeNotification(notification_id())
            if (length(result_table) > 0){
              result_table
            } else if(length(result_table == 0)) {
              textOutput("No samples have mutations on this gene")
            }
          })
      })
    })
  })
  output$download_table <- downloadHandler(
      filename = function() {
        "Variant_calling_result_table.csv"
      },
      content = function(file) {
        write.csv(result_table, file)
      }
    )
  
  output$download_raw <- downloadHandler(
      filename = function() {
        "SNP_samples.csv"
      },
      content = function(file) {
        write.csv(countList, file)
      }
    )
}
# Create Shiny app ----
shinyApp(ui, server)