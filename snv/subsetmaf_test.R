
# console test works well -------------------------------------------------

ct <- reactive({c("LUAD")}) # pre define a object stored a cancer type.
query = as.expression("cancer_types %in% isolate(ct())") # define a query for filter
maftools::subsetMaf(mc3_pass, query = query, mafObj = T) -> gene_list_maf # subsetMaf filter

gene_list_maf %>%
  maftools::getClinicalData() # check the filter results
Tumor_Sample_Barcode cancer_types
1: TCGA-05-4244-01A-01D-1105-08         LUAD
2: TCGA-05-4249-01A-01D-1105-08         LUAD
3: TCGA-05-4250-01A-01D-1105-08         LUAD
4: TCGA-05-4382-01A-01D-1931-08         LUAD
5: TCGA-05-4384-01A-01D-1753-08         LUAD
---                                          
  508: TCGA-NJ-A55O-01A-11D-A25L-08         LUAD
509: TCGA-NJ-A55R-01A-11D-A25L-08         LUAD
510: TCGA-NJ-A7XG-01A-12D-A397-08         LUAD
511: TCGA-O1-A52J-01A-11D-A25L-08         LUAD
512: TCGA-S2-AA1A-01A-12D-A397-08         LUAD


# shiny active value lead an error ----------------------------------------

library(shiny)
mc3_pass <- readr::read_rds(file.path(config$database, "TCGA", "snv", "snv_mutation_mc3_public.pass.filtered_maf-4cancers.rds.gz"))

Lung_choice <- list(
  "Lung Adenocarcinoma(LUAD)" = "LUAD",
  "Lung Squamous Cell Carcinoma(LUSC)" = "LUSC"
)
ui <- fluidPage(
  checkboxGroupInput(inputId = "select",label = "Select cancers",
                     choices = Lung_choice),
  actionButton("go",label = "Go"),
  shiny::textOutput(outputId = "table")
)

server <- function(input, output, session) {
  ct <- reactive({
    input$select
  })
  observeEvent(input$go,{
    print(ct())
    print(class(ct()))
    ct_2 <- as.character(ct())
    print(ct_2)
    # query = as.expression("cancer_types %in% ct()") # define a query for filter
    query = as.expression("cancer_types %in% ct_2") # ct(),ct_2,input$select all lead an same error in subsetMaf(...,query = query)
    print(query)
    
    maftools::subsetMaf(mc3_pass, query = query, mafObj = T) -> gene_list_maf # subsetMaf filter
    print("maf done!")
    gene_list_maf %>%
      maftools::getClinicalData() %>%
      dplyr::select(cancer_types) %>%
      unique() %>% t() %>%
      as.vector()->.x# check the filter results
    print(".x done!")
    output$table <- renderText({
      return(.x)
    })
    print("out done!")
    print(.x)
  })
  
}

shinyApp(ui, server)
