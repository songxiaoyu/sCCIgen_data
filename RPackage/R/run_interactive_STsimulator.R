#' Run interactive STsimulator
#'
#' @import dplyr
#' @return A simulated dataset. Either as individual csv files with counts, cell features, or a Giotto object.
#'
#' @export
#'
run_interactive_STsimulator <- function() {

  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("STsimulator"),
    miniUI::miniContentPanel(
      shiny::tabsetPanel(
        type = "tabs",

        shiny::tabPanel("Download a pre-simulated dataset",
                        shiny::radioButtons(inputId = "dataset",
                                            label = "Select a pre-simulated spatial transcriptomics dataset to download.
                                            These pre-simulated data include 5% spatial differential expressed genes in one cell
                                            type, 5% interaction changed genes in one cell type that are associated with the
                                            proximity of another cell type, and 5% interaction changed genes in cell-type pair that
                                            are associated with each other in the neighboring cells.",
                                            choices = c("Simulated data for 4,751 genes in 4,000 cells of 6 cell types
                                                        in two regions on a unit square, using normal breast data profiled by
                                                        snRNAseq as reference." = "example1",
                                                        "Simulated data for 10,000 genes in 500 spots of 6 cell types, using
                                                        mouse brain data profiled by SeqFISH+ as reference." = "example2",
                                                        "Simulated data for 550 genes in 10,000 cells of 6 cell types, using
                                                        human ovarian cancer data profiled by MERFISH as reference." = "example3"),
                                            width = "100%",
                                            selected = character(0)
                        ), # ends radioButtons for dataset

                        shiny::downloadButton("downloadCount", "Download Counts"),
                        shiny::downloadButton("downloadMeta", "Download Metadata"),
                        shiny::downloadButton("downloadExpr", "Download Expression Pattern"),
                        shiny::downloadButton("downloadParameter", "Download Parameter file")
        ), # ends Download a pre-simulated dataset tab

        shiny::tabPanel("Run a simulation",
                        shiny::tabsetPanel(
                          type = "tabs",
                          shiny::tabPanel("Step 1: Select your dataset",
                                          shiny::radioButtons(inputId = "inputdata",
                                                              label = "What data do you want to use for simulation? If the data
                                                       is not yet in your working directory, use the Download buttons
                                                       to get the files in your computer.",
                                                       choices = c("Decoy data 1: It includes (1) count matrix for 10 genes by
                                                                   1000 cells of 2 cell types, and (2) cell feature matrix for
                                                                   annotated cell type." = "fake1",
                                                                   "Decoy data 2: It includes (1) count matrix for 10 genes by
                                                                   1000 cells of 2 cell types, and (2) cell feature matrix for
                                                                   annotated cell type and spatial coordinate." = "fake2",
                                                                   "Decoy data 3: It includes (1) count matrix for 10 genes by
                                                                   1000 cells of 2 cell types, and (2) cell feature matrix for
                                                                   annotated cell type, spatial coordinate, and region." = "fake3",
                                                                   "Normal human breast snRNAseq data: It includes (1) count matrix
                                                                   for 4751 genes by 5990 cells of 6 cell types (epithelial cell,
                                                                   adipocyte, fibroblast, endothelial cell, immune (myeloid) and
                                                                   muscle), and (2) cell feature matrix for annotated cell type.
                                                                   PMID: 35549429" = "snRNAseq_breast_033023",
                                                                   "Normal mouse brain SeqFISH+ data: It includes (1) count matrix
                                                                   for 10,000 genes by 511 cells of 6 cell types (excitatory neuron,
                                                                   interneuron, astrocyte, microglia, oligodendrocyte and
                                                                   endothelial cells), and (2)cell feature matrix including cell
                                                                   type annotation and spatial coordinate on 2D (x, y).
                                                                   PMID: 35549429" = "SeqFishPlusCortex_033023",
                                                                   "Ovarian cancer MERFISH data: It includes (1) count matrix
                                                                   for 550 genes by 355,633 cells of 6 cell types (tumor, adipose,
                                                                   endo-immune, macro-immune, macrophage, and others), and (2) cell
                                                                   feature matrix including cell type annotation and spatial
                                                                   coordinate on 2D. Data source: vizgen" = "MERFISH_OV",
                                                                   "I want to use my own dataset" = "user_input"),
                                                       width = "100%",
                                                       selected = character(0)
                                          ), # ends radioButtons for inputdata

                                          shiny::uiOutput("downloadExpression"),

                                          shiny::uiOutput("downloadFeature"),

                                          shiny::HTML(strrep(htmltools::br(), 2)),

                                          shiny::uiOutput("userinputexpression"),

                                          shiny::uiOutput("userinputcellfeature"),

                          ), # ends Step 1: Select your dataset

                          shiny::tabPanel("Step 2: Select the simulation parameters",

                                          shiny::textInput(inputId = "inputdir",
                                                           label = "Where is your input data located? By default, we will use your working directory. If you want to specify a diferent folder, type the path here. Note that the folder must be under your working directory.",
                                                           width = "100%"),

                                          shiny::uiOutput("ask_simulate_cells"),

                                          shiny::uiOutput("ask_numbercells"),

                                          shiny::uiOutput("ask_numberregions"),

                                          shiny::uiOutput("ask_custom_cell_type_prop"),

                                          shiny::uiOutput("customcellproportionsSelection"),

                                          shiny::uiOutput("button_customcellproportions"),

                                          shiny::uiOutput("ask_custom_interaction"),

                                          shiny::uiOutput("custominteractionSelection"),

                                          shiny::uiOutput("button_custominteraction"),

                                          shiny::uiOutput("ask_even_distribution"),

                                          shiny::uiOutput("ask_windowmethod"),

                                          shiny::uiOutput("ask_overlapcutoff"),

                                          # Parameters for expression profiles

                                          shiny::numericInput(inputId = "depthratio",
                                                              label = "Select the sequencing depth ratio between
                                                              simulated and reference data. By default, simulated data
                                                              have the same depth as the input data; an alteration of this
                                                              parameter is a way to generate batch effects for different simulations.",
                                                              value = 1,
                                                              width = "100%"),

                                          shiny::uiOutput("ask_model_per_region"),

                                          shiny::uiOutput("ask_custom_model_per_region"),

                                          shiny::uiOutput("button_custom_model_per_region"),

                                          shiny::radioButtons(inputId = "mimiccorrelation",
                                                              label = "Do you want to mimic gene-gene correlation of
                                                              the reference data? Select 'Yes' only if the gene-gene
                                                              correlation is pre-estimated for each cell type.
                                                              STsimulator provides functions/pipelines for estimating
                                                              and saving gene-gene correlation, which can be served
                                                              as input here.",
                                                              choices = c("No, simulate independent genes" = FALSE,
                                                                          "Yes" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("mimiccorrelationSelection"),

                                          shiny::radioButtons(inputId = "spatialpatterns",
                                                              label = "Would you like to add spatial patterns?",
                                                              choices = c("No" = FALSE,
                                                                          "Yes" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("spatialpatternsSelection"),

                                          shiny::uiOutput("button_spatialpatterns"),

                                          shiny::uiOutput("ask_addcellcellinteraction"),

                                          shiny::uiOutput("addcellcellinteractionSelection"),

                                          shiny::uiOutput("button_addcellcellinteraction"),

                                          shiny::uiOutput("ask_addcellcellinteractionexpression"),

                                          shiny::uiOutput("addcellcellinteractionexpressionSelection"),

                                          shiny::uiOutput("button_addcellcellinteractionexpression"),

                                          shiny::radioButtons(inputId = "multicell",
                                                              label = "What resolution would you like to use to simulate your data?",
                                                              choices = c("single-cell resolution" = FALSE,
                                                                          "multi-cell resolution" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("multicellSelection"),

                                          shiny::numericInput(inputId = "ndatasets",
                                                              label = "Select the number of simulated data sets",
                                                              value = 1,
                                                              width = "100%"),

                                          shiny::numericInput(inputId = "seed",
                                                              label = "Specify an umbrella seed for reproducible simulation.
                                                              Seed for each data set will be generated using the umbrella seed.",
                                                              value = 1234,
                                                              width = "100%"),

                                          shiny::textInput(inputId = "outdir",
                                                           label = "Provide the path to save the output: [default = working_directory/output_files/]",
                                                           value = "output_files",
                                                           width = "100%"),

                                          shiny::textInput(inputId = "outname",
                                                           label = "Type the name of the parameter file.",
                                                           value = "parameter_file.csv",
                                                           width = "100%"),

                                          shiny::downloadButton("downloadparams", "Create and Download my parameter file")

                          ), # ends Step 2: Select the simulation parameters

                          shiny::tabPanel("Step 3: Run the simulation",

                                          shiny::textInput(inputId = "inputparamfile",
                                                           label = "By default, we will look for a parameters_file.tsv file. If your file has a different name, type it here:",
                                                           width = "100%",
                                                           value = "parameters_file.tsv"),

                                          shiny::radioButtons(inputId = "createGiotto",
                                                              label = "By the default, the simulation will create the expression, metadata, and pattern output files. Optionally, you can also create a Giotto object.",
                                                              choices = c("Run the simulation and create the default output files" = FALSE,
                                                                          "Run the simulation and create a Giotto object" = TRUE),
                                                              width = "100%"),

                                          shiny::uiOutput("createGiottoSelection"),

                                          shiny::actionButton(inputId = "runsimulation",
                                                              label = "Run simulation",
                                                              width = "100%")

                          ) # ends Step 3: Run the simulation

                        ) # ends tabsetPanel inside Run a simulation


        ) # ends Run a simulation
      ) # ends tabsetPanel



    ) # ends miniContentPanel
  ) # ends ui

  server <- function(input, output,session) {

    # Download pre-simulated data
    shiny::observeEvent(input$dataset, {

      output$downloadCount <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_count.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/STsimulator/raw/main/DockerUI/example_data/",input$dataset,"_count.tsv"),
                        destfile = con)
        }
      )

      output$downloadMeta <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_meta.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/STsimulator/raw/main/DockerUI/example_data/",input$dataset,"_meta.tsv"),
                        destfile = con)
        }
      )

      output$downloadExpr <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_expr_pattern.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/STsimulator/raw/main/DockerUI/example_data/",input$dataset,"_expr_pattern.tsv"),
                        destfile = con)
        }
      )

      output$downloadParameter <- shiny::downloadHandler(
        filename = paste0(input$dataset,"_parameter.tsv"),
        content = function(con) {
          download.file(url = paste0("https://github.com/songxiaoyu/STsimulator/raw/main/DockerUI/example_data/",input$dataset,"_parameter.tsv"),
                        destfile = con)
        }
      )

    }) # ends download pre-simulated data


    # Download Rdata
    shiny::observeEvent(input$inputdata, {

      if (input$inputdata != "user_input") {

        # show download buttons
        output$downloadExpression <-  renderUI({
          shiny::downloadButton("downloadselectedexpression", "Download selected expression data")
        })

        output$downloadFeature <-  renderUI({
          shiny::downloadButton("downloadselectedcellfeature", "Download selected cellfeature data")
        })

        output$downloadselectedexpression <- shiny::downloadHandler(
          filename = paste0(input$inputdata,"_expr.Rdata"),
          content = function(con) {
            download.file(url = paste0("https://github.com/josschavezf/STsimulator/raw/main/DockerUI/InputData/expression_data/",input$inputdata,"_expr.Rdata"),
                          destfile = con)
          })

        output$downloadselectedcellfeature <- shiny::downloadHandler(
          filename = paste0(input$inputdata,"_cellfeature.Rdata"),
          content = function(con) {
            download.file(url = paste0("https://github.com/josschavezf/STsimulator/raw/main/DockerUI/InputData/cell_feature_data/",input$inputdata,"_cellfeature.Rdata"),
                          destfile = con)
          })
      }
    }) # ends download data

    # Create parameters file

    path_to_input_dir <- shiny::reactiveVal()

    shiny::observeEvent(input$inputdir, {
      if(is.null(input$inputdir)) {
        path_to_input_dir(getwd())
      } else {path_to_input_dir(here::here(input$inputdir)) }
    })

    expression_data_file <- shiny::reactiveVal()
    expression_data_file_type <- shiny::reactiveVal()
    cell_feature_data_file <- shiny::reactiveVal()
    expression_data_cell_types <- shiny::reactiveVal()
    expr_data <- shiny::reactiveVal()
    cellfeature_data <- shiny::reactiveVal()
    ncol_feature_data <- shiny::reactiveVal()

    shiny::observeEvent(input$inputdata, {

      if(input$inputdata == "user_input") {

        output$userinputexpression <- renderUI({
          shiny::textInput(inputId = "user_expression",
                           label = "If you're using your own dataset, type the name of the expression file (e.g. expression.Rdata, .tsv, or .csv). Upload a G gene by N cell matrix for expression count data. Instruction: Row names should be unique identifiers of genes.",
                           width = "100%")
        })

        output$userinputcellfeature <- renderUI({
          shiny::textInput(inputId = "user_cellfeature",
                           label = "If you're using your own odataset, type the name of the cellfeature file (e.g. cellfeature.Rdata, .tsv, or .csv). Upload a N by K matrix for cell feature data. Instruction: Column and row names are not expected. Cells should be in the same order as the uploaded expression count data. The first column is required, which should be cell type annotation. Other columns are optional.  If spatial coordinates on x, y axes are provided, they should be the 2-3 columns of the input cell feature data. If a spatial region variable is provided, it should be the 4th column of the data.",
                           width = "100%")
        })

        shiny::observeEvent(input$user_expression, {
          expression_data_file(input$user_expression)
        })

        shiny::observeEvent(input$user_cellfeature, {
          cell_feature_data_file(input$user_cellfeature)
        })

      } else {
        expression_data_file(paste0(input$inputdata,"_expr.Rdata"))
        cell_feature_data_file(paste0(input$inputdata,"_cellfeature.Rdata"))
      }

      if(file.exists(fs::path(path_to_input_dir(), expression_data_file())) &
         file.exists(fs::path(path_to_input_dir(), cell_feature_data_file())) ) {

        expression_data_file_type(tools::file_ext(expression_data_file()))
        if(expression_data_file_type() == "Rdata" ||
           expression_data_file_type() == "RData") {

          load(fs::path(path_to_input_dir(), expression_data_file()), x <- new.env())
          expr_data(get(ls(x), envir = x))
          rm(x)

          load(fs::path(path_to_input_dir(), cell_feature_data_file()), x <- new.env())
          y = get(ls(x), envir = x)
          rm(x)

          colnames(y)[1] = "cell_type"
          cellfeature_data(y)

        } else if(expression_data_file_type() == "csv") {

          expr_data(read.csv(fs::path(path_to_input_dir(), expression_data_file()),
                             row.names = 1))

          cellfeature_data(read.csv(fs::path(path_to_input_dir(), cell_feature_data_file()),
                                    row.names = 1))

        } else if(expression_data_file_type() == "tsv") {

          expr_data(read.delim(fs::path(path_to_input_dir(), expression_data_file()),
                               row.names = 1))

          cellfeature_data(read.delim(fs::path(path_to_input_dir(), cell_feature_data_file()),
                                      row.names = 1))
        }

        expression_data_cell_types(paste(unique(sort(unique(cellfeature_data()[,1]) )),
                                         collapse = ","))

        ncol_feature_data(ncol(cellfeature_data()))

      }

    })

    simulate_spatial_data <- shiny::reactiveVal()

    window_method <- shiny::reactiveVal()

    num_regions <- shiny::reactiveVal()
    custom_cell_type_proportions <- shiny::reactiveVal()
    cell_type_proportions <- shiny::reactiveVal()
    custom_props <- shiny::reactiveVal()
    custom_cell_location_interactions <- shiny::reactiveVal()
    cell_location_interactions_df <- shiny::reactiveVal()
    cell_even_distribution <- shiny::reactiveVal()
    region_specific_model <- shiny::reactiveVal()

    shiny::observeEvent(ncol_feature_data(), {

      if(ncol_feature_data() > 1) {

        output$ask_simulate_cells <- shiny::renderUI({
          shiny::radioButtons(inputId = "simulatecells",
                              label = "Do you want to simulate new cells?",
                              choices = c("Do not simulate new cells, but simulate new expression data
                                          for existing cells" = FALSE,
                                          "I want to simulate new cells" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$simulatecells, {

          simulate_spatial_data(input$simulatecells)

          if(input$simulatecells == TRUE) {

            output$ask_windowmethod <- shiny::renderUI({
              shiny::radioButtons(inputId = "windowmethod",
                                  label = "Select the method for determining the window on existing ST data",
                                  choices = c("network" = "network",
                                              "rectangle" = "rectangle",
                                              "convex" = "convex",
                                              "convex2" = "convex2",
                                              "convex3" = "convex3",
                                              "convex5" = "convex5"))
            })

            shiny::observeEvent(input$windowmethod, {
              window_method(input$windowmethod)
            })

          }
        })

        if(ncol_feature_data() > 3) {
          unique_regions = unique(sort(cellfeature_data()[,4]))

          num_regions(length(unique_regions))

          region_specific_model("NULL")

        } else {
          output$ask_numberregions <- shiny::renderUI({
            shiny::numericInput(inputId = "numberregions",
                                label = "Enter the number of regions (suggested: 1-10)",
                                value = 2,
                                width = "100%")
          })

          shiny::observeEvent(input$numberregions, {
            num_regions(input$numberregions)

            if(input$numberregions > 1) {
              output$ask_model_per_region <- shiny::renderUI({
                shiny::radioButtons(inputId = "model_per_region",
                                    label = "Do you want to model the input expression data
                                  separately for each region?",
                                  choices = c("No" = FALSE,
                                              "Yes" = TRUE),
                                  width = "100%")
              })

              shiny::observeEvent(input$model_per_region, {

                if(input$model_per_region == TRUE) {
                  output$ask_custom_model_per_region <- shiny::renderUI({
                    shiny::textInput(inputId = "custom_model_per_region",
                                     label = "Provide the model of input expression data for each region in the format: <><> (e.g. )",
                                     width = "100%")
                  })

                  output$button_custom_model_per_region <- shiny::renderUI({
                    shiny::actionButton("save_custom_model_per_region", "Use these expression models")
                  })

                  shiny::observeEvent(input$save_custom_model_per_region, {
                    region_specific_model(input$custom_model_per_region)
                  })

                } else {region_specific_model("NULL")}
              })
            }
          })

        }

      } else { # when ncol == 1
        simulate_spatial_data(TRUE)

        output$ask_numberregions <- shiny::renderUI({
          shiny::numericInput(inputId = "numberregions",
                              label = "Enter the number of regions (suggested: 1-10)",
                              value = 2,
                              width = "100%")
        })

        shiny::observeEvent(input$numberregions, {
          num_regions(input$numberregions)

          if(input$numberregions > 1) {
            output$ask_model_per_region <- shiny::renderUI({
              shiny::radioButtons(inputId = "model_per_region",
                                  label = "Do you want to model the input expression data
                                  separately for each region?",
                                  choices = c("No" = FALSE,
                                              "Yes" = TRUE),
                                  width = "100%")
            })

            shiny::observeEvent(input$model_per_region, {

              if(input$model_per_region == TRUE) {
                output$ask_custom_model_per_region <- shiny::renderUI({
                  shiny::textInput(inputId = "custom_model_per_region",
                                   label = "Provide the model of input expression data for each region in the format: <><> (e.g. )",
                                   width = "100%")
                })

                output$button_custom_model_per_region <- shiny::renderUI({
                  shiny::actionButton("save_custom_model_per_region", "Use these expression models")
                })

                shiny::observeEvent(input$save_custom_model_per_region, {
                  region_specific_model(input$custom_model_per_region)
                })

              } else {region_specific_model("NULL")}
            })
          }
        })

        output$ask_custom_cell_type_prop <- shiny::renderUI({
          shiny::radioButtons(inputId = "customcellproportions",
                              label = "Choose an option for cell-type proportions in each region.",
                              choices = c("Use all equals to cell-type proportions of the input expression data" = FALSE,
                                          "I want to input custom cell-type proportions" = TRUE),
                              width = "100%"
          )
        })

        shiny::observeEvent(input$customcellproportions, {

          custom_props(input$customcellproportions)

          if(input$customcellproportions == TRUE) { # if TRUE, use custom cell proportions

            output$customcellproportionsSelection <- shiny::renderUI({
              shiny::textInput(inputId = "cellproportions",
                               label = "Enter the custom cell-type proportions in a format <region>,<cell_type>,<proportion>. Separate multiple entries by blank space (e.g '1,cell_type_A,0.25 1,cell_type_B,0.75 2,cell_type_C,1.0'). NOTE: Cell type proportions for a given region must sum up to 1.",
                               width = "100%")
            })

            output$button_customcellproportions <- shiny::renderUI({
              shiny::actionButton("save_customcellproportions", "Use these cell proportions")
            })

            shiny::observeEvent(input$save_customcellproportions, {

              shiny::observeEvent(input$cellproportions, {
                x = unlist(stringr::str_split(string = input$cellproportions,
                                              pattern = " "))

                n_total = length(x)
                cell_type_proportions(data.frame(parameters = paste0("cell_type_proportion_", 1:n_total),
                                                 value = x))
              })

            })

          } else {
            x = cellfeature_data() %>%
              count(cell_type) %>%
              mutate(proportions = round(n/sum(n),3),
                     proportions2 = paste0(cell_type,",",proportions))

            n_regions = input$numberregions
            n_total = n_regions*length(x$proportions2)

            cell_type_proportions(data.frame(parameters = paste0("cell_type_proportion_", 1:n_total),
                                             value = paste0(base::sort(rep(1:n_regions),length(x$proportions2)), ",", x$proportions2)) )
          }

        })

        output$ask_custom_interaction <- shiny::renderUI({
          shiny::radioButtons(inputId = "custominteraction",
                              label = "Users can select cell type pairs and determine the strength. Strength < 0 indicates cell-cell inhibition; and strength > 0 indicates cell-cell attraction. Would you like to specify the cell-cell location interaction (inhibition/attraction)?",
                              choices = c("Use the default cell-cell location interaction (inhibition/attraction)" = FALSE,
                                          "Specify cell-cell location interaction (inhibition/attraction)" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$custominteraction, {

          custom_cell_location_interactions(input$custominteraction)

          if(input$custominteraction == TRUE) {

            output$custominteractionSelection <- shiny::renderUI({
              shiny::textInput(inputId = "locationinteraction",
                               label = "To specify cell-cell location interaction, enter the cell type pairs followed by the interaction level (suggested value -2 to 2) in a format <cell_type_A>-<cell_type_B>,<value>. Separate the entries by blank space (e.g. 'cell_type_A-cell_type_B,1.2 cell_type_B-cell_type_C,-0.8')",
                               width = "100%")
            })

            output$button_custominteraction <- shiny::renderUI({
              shiny::actionButton("save_custominteraction", "Use these interactions")
            })

            shiny::observeEvent(input$save_custominteraction, {

              shiny::observeEvent(input$locationinteraction, {
                x = unlist(stringr::str_split(input$locationinteraction,
                                              pattern = " "))

                cell_location_interactions_df(data.frame(parameters = paste0("cell_interaction_",
                                                                             1:length(x)),
                                                         value = x))
              })

            })

          }
        })

        output$ask_even_distribution <- shiny::renderUI({
          shiny::sliderInput(inputId = "evendistribution",
                             label = "Adjusts cells to be evenly distributed on a slide. 0 = uneven distribution, 1 = even distribution.",
                             min = 0,
                             max = 1,
                             value = 0,
                             width = "100%")
        })

        shiny::observeEvent(input$evendistribution, {
          cell_even_distribution(input$evendistribution)
        })

      }
    })

    num_simulated_cells <- shiny::reactiveVal()

    cell_overlap_cutoff <- shiny::reactiveVal()


    shiny::observeEvent(simulate_spatial_data(), {

      if(simulate_spatial_data() == TRUE) {

        output$ask_numbercells <- shiny::renderUI({
          shiny::numericInput(inputId = "numbercells",
                              label = "Enter the number of simulated cells",
                              value = 10000,
                              width = "100%")
        })

        shiny::observeEvent(input$numbercells, {
          num_simulated_cells(input$numbercells)
        })

        output$ask_overlapcutoff <- shiny::renderUI({
          shiny::sliderInput(inputId = "overlapcutoff",
                             label = "Select the overlap cutoff. Cells closer than this cutoff will be considered overlapping to each other and all but one will be removed (range: 0-0.1 of the slide length/width)",
                             min = 0,
                             max = 0.1,
                             value = 0.02,
                             width = "100%")
        })

        shiny::observeEvent(input$overlapcutoff, {
          cell_overlap_cutoff(input$overlapcutoff)
        })

      }

    })

    # Parameters for expression profiles

    expr_depth_ratio <- shiny::reactiveVal()

    shiny::observeEvent(input$depthratio, {
      expr_depth_ratio(input$depthratio)
    })

    gene_cor <- shiny::reactiveVal()
    copula_input <- shiny::reactiveVal()

    shiny::observeEvent(input$mimiccorrelation, {

      gene_cor(input$mimiccorrelation)

      if(input$mimiccorrelation == TRUE) {
        output$mimiccorrelationSelection <- shiny::renderUI({
          shiny::textInput(inputId = "genecorfile",
                           label = "Provide a path for uploading a pre-estimated gene-gene correlation file:",
                           width = "100%")
        })

        shiny::observeEvent(input$genecorfile, {
          copula_input(input$genecorfile)
        })

      }
    })

    spatialpatterns_df <- shiny::reactiveVal()

    shiny::observeEvent(input$spatialpatterns, {

      if(input$spatialpatterns == TRUE) {
        output$spatialpatternsSelection <- shiny::renderUI({
          shiny::textInput(inputId = "inputspatialpatterns",
                           label = "Specify the spatial patterns in the format
                           <region>,<cell_type>,<gene_id(optional)>,<gene_proportion>,
                           <Mean effect at log(count) scale>,<SD of effect at log(count) scale>.
                           Separate the entries by blank space (e.g. '1,cell_type_A,NULL,0.1,0.5,0 1,cell_type_A,gene_A,0.2,0.5,0.1')",
                           width = "100%")
        })

        output$button_spatialpatterns <- shiny::renderUI({
          shiny::actionButton("save_spatialpatterns", "Use these spatial patterns")
        })

        shiny::observeEvent(input$save_spatialpatterns, {

          shiny::observeEvent(input$inputspatialpatterns, {
            x_vector = unlist(stringr::str_split(input$inputspatialpatterns,
                                                 pattern = " "))

            x_df <- data.frame(parameters = character(),
                               value = character())

            for (i in 1:length(x_vector)) {
              x = data.frame(parameters = c(paste0("spatial_pattern_",i,"_region"),
                                            paste0("spatial_pattern_",i,"_cell_type"),
                                            paste0("spatial_pattern_",i,"_gene_id"),
                                            paste0("spatial_pattern_",i,"_gene_prop"),
                                            paste0("spatial_pattern_",i,"_mean"),
                                            paste0("spatial_pattern_",i,"_sd")),
                             value = unlist(stringr::str_split(x_vector[i],
                                                               pattern = ",")) )

              x_df = rbind(x_df, x)
            }

            spatialpatterns_df(x_df)

          })

        })

      }
    })


    cellcellinteractions_df <- shiny::reactiveVal()
    cellcellinteractionsexpression_df <- shiny::reactiveVal()


    shiny::observeEvent(input$numberregions, {

      if(input$numberregions > 1) {

        output$ask_addcellcellinteraction <- shiny::renderUI({
          shiny::radioButtons(inputId = "addcellcellinteraction",
                              label = "Would you like to add cell-cell interactions - expression
                                                              associated with cell-cell distance?",
                              choices = c("No" = FALSE,
                                          "Yes" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$addcellcellinteraction, {

          if(input$addcellcellinteraction == TRUE) {

            output$addcellcellinteractionSelection <- shiny::renderUI({
              shiny::textInput(inputId = "cellcellinteractions",
                               label = "Specify the <region>,<perturbed cell>,<adjacent cell>,<interaction
                           distance threshold (default 0.1)>,<gene ID (optional)>,<proportion of
                           the gene or proportion of genes with the pattern>,<Mean effect at log(count)
                           scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>.
                           Separate the entries by blank space (e.g. '1,cell_A,cell_B,0.1,NULL,0.1,0.5,0 2,cell_B,cell_C,0.2,gene_A,0.2,0.6,0.1').
                           NOTE: Adjacent cell type cannot be the same as peturbed cell type.",
                           width = "100%")
            })

            output$button_addcellcellinteraction <- shiny::renderUI({
              shiny::actionButton("save_addcellcellinteraction", "Use these cell-cell interactions")
            })

            shiny::observeEvent(input$save_addcellcellinteraction, {

              shiny::observeEvent(input$cellcellinteractions, {
                x_vector = unlist(stringr::str_split(input$cellcellinteractions,
                                                     pattern = " "))

                x_df <- data.frame(parameters = character(),
                                   value = character())

                for (i in 1:length(x_vector)) {
                  x = data.frame(parameters = c(paste0("spatial_int_dist_",i,"_region"),
                                                paste0("spatial_int_dist_",i,"_cell_type_perturbed"),
                                                paste0("spatial_int_dist_",i,"_cell_type_adj"),
                                                paste0("spatial_int_dist_",i,"_dist_cutoff"),
                                                paste0("spatial_int_dist_",i,"_gene_id1"),
                                                paste0("spatial_int_dist_",i,"_gene_prop"),
                                                paste0("spatial_int_dist_",i,"_mean"),
                                                paste0("spatial_int_dist_",i,"_sd")),
                                 value = unlist(stringr::str_split(x_vector[i],
                                                                   pattern = ",")) )

                  x_df = rbind(x_df, x)
                }

                cellcellinteractions_df(x_df)
              })
            })

          }
        })

        output$ask_addcellcellinteractionexpression <- shiny::renderUI({
          shiny::radioButtons(inputId = "addcellcellinteractionexpression",
                              label = "Would you like to add cell-cell interactions - expression
                                      associated with expression of neighboring cells??",
                              choices = c("No" = FALSE,
                                          "Yes" = TRUE),
                              width = "100%")
        })

        shiny::observeEvent(input$addcellcellinteractionexpression, {

          if(input$addcellcellinteractionexpression == TRUE) {

            output$addcellcellinteractionexpressionSelection <- shiny::renderUI({
              shiny::textInput(inputId = "cellcellinteractionsexpression",
                               label = "Specify the <region>,<perturbed cell>,<adjacent cell>,
                           <interaction distance threshold (default 0.1)>,<gene ID 1 (optional)>,
                           <gene ID 2 (optional)>,<proportion of genes with the pattern>,
                           <Bidirectional association (TRUE or FALSE)>,<Mean effect at log(count)
                           scale (default = 0.5)>,<SD of effect at log(count) scale (default = 0)>.
                           Separate the entries by blank space (e.g. '1,cell_A,cell_B,0.1,NULL,NULL,0.1,TRUE,0.5,0 2,cell_B,cell_C,0.2,gene_A,gene_B,0.2,FALSE,0.6,0.1').
                           NOTE: Adjacent cell type cannot be the same as peturbed cell type.",
                           width = "100%")
            })

            output$button_addcellcellinteractionexpression <- shiny::renderUI({
              shiny::actionButton("save_addcellcellinteractionexpression", "Use these cell-cell interactions")
            })

            shiny::observeEvent(input$save_addcellcellinteractionexpression, {

              shiny::observeEvent(input$cellcellinteractionsexpression, {
                x_vector = unlist(stringr::str_split(input$cellcellinteractionsexpression,
                                                     pattern = " "))

                x_df <- data.frame(parameters = character(),
                                   value = character())

                for (i in 1:length(x_vector)) {
                  x = data.frame(parameters = c(paste0("spatial_int_expr_",i,"_region"),
                                                paste0("spatial_int_expr_",i,"_cell_type_perturbed"),
                                                paste0("spatial_int_expr_",i,"_cell_type_adj"),
                                                paste0("spatial_int_expr_",i,"_dist_cutoff"),
                                                paste0("spatial_int_expr_",i,"_gene_id1"),
                                                paste0("spatial_int_expr_",i,"_gene_id2"),
                                                paste0("spatial_int_expr_",i,"_gene_prop"),
                                                paste0("spatial_int_expr_",i,"_bidirectional"),
                                                paste0("spatial_int_expr_",i,"_mean"),
                                                paste0("spatial_int_expr_",i,"_sd")),
                                 value = unlist(stringr::str_split(x_vector[i],
                                                                   pattern = ",")) )

                  x_df = rbind(x_df, x)
                }

                cellcellinteractionsexpression_df(x_df)
              })
            })

          }
        })

      }
    })


    num_spots <- shiny::reactiveVal()

    shiny::observeEvent(input$multicell, {

      if(input$multicell == TRUE) {

        output$multicellSelection <- renderUI({
          shiny::numericInput(inputId = "nspots",
                              label = "Specify the number of spots",
                              value = 1,
                              width = "100%")
        })

        shiny::observeEvent(input$nspots, {
          num_spots(input$nspots)
        })
      } else {num_spots("NULL")}
    })

    num_simulated_datasets <- shiny::reactiveVal()

    shiny::observeEvent(input$ndatasets, {
      num_simulated_datasets(input$ndatasets)
    })

    parent_simulation_seed <- shiny::reactiveVal()
    simulation_seed_for_each_dataset <- shiny::reactiveVal()

    shiny::observeEvent(input$seed, {
      parent_simulation_seed(input$seed)
      simulation_seed_for_each_dataset(input$seed)
    })

    path_to_output_dir <- shiny::reactiveVal()

    shiny::observeEvent(input$outdir, {

      if(input$outdir == "output_files") {
        path_to_output_dir(here::here("output_files"))
      } else { path_to_output_dir(input$outdir) }
    })

    output_name <- shiny::reactiveVal()

    shiny::observeEvent(input$outname, {
      output_name(input$outname)
    })



    # export parameters file
    output$downloadparams <- shiny::downloadHandler(
      filename = output_name(),
      content = function(con) {

        param_df <- data.frame(parameters = c("path_to_input_dir",
                                              "expression_data_file",
                                              "cell_feature_data_file",
                                              "expression_data_file_type",
                                              "expression_data_cell_types",
                                              "simulate_spatial_data"
        ),
        value = c(path_to_input_dir(),
                  expression_data_file(),
                  cell_feature_data_file(),
                  expression_data_file_type(),
                  expression_data_cell_types(),
                  simulate_spatial_data()
        )
        )

        if(simulate_spatial_data() == TRUE) {
          df = data.frame(parameters = c("num_simulated_cells",
                                         "cell_overlap_cutoff"),
                          value = c(num_simulated_cells(),
                                    cell_overlap_cutoff()
                          )
          )

          param_df = rbind(param_df, df)
        }

        if(ncol_feature_data() > 1 & simulate_spatial_data() == TRUE) {
          param_df = rbind(param_df, c("window_method", window_method()))
        }

        if(!is.null(num_regions())) {
          param_df = rbind(param_df, c("num_regions", num_regions()))

          if(num_regions() > 1) {
            param_df = rbind(param_df, c("region_specific_model",
                                         region_specific_model()))
          }
        }

        if(ncol_feature_data() == 1) {


          param_df = rbind(param_df, c("custom_cell_type_proportions",
                                       custom_props()))

          param_df = rbind(param_df, cell_type_proportions())

          param_df = rbind(param_df, c("custom_cell_location_interactions",
                                       custom_cell_location_interactions() ))

          if(custom_cell_location_interactions() == TRUE) {
            param_df = rbind(param_df, cell_location_interactions_df())
          }

          param_df = rbind(param_df,c("cell_even_distribution",
                                      cell_even_distribution() ))
        }


        # simulation parameters for expression profiles
        param_df = rbind(param_df, c("expr_depth_ratio",
                                     expr_depth_ratio() ))

        param_df = rbind(param_df, c("gene_cor",
                                     gene_cor() ))

        if(gene_cor() == TRUE) {
          param_df = rbind(param_df, c("copula_input",
                                       copula_input() ))
        }


        if(!is.null(spatialpatterns_df())) {
          param_df = rbind(param_df, spatialpatterns_df())
        }


        if(!is.null(num_regions()) & num_regions() > 1) {
          if(!is.null(cellcellinteractions_df())) {
            param_df = rbind(param_df, cellcellinteractions_df())
          }

          if(!is.null(cellcellinteractionsexpression_df())) {
            param_df = rbind(param_df, cellcellinteractionsexpression_df())
          }

        } else {
          param_df = rbind(param_df, c("region", "NULL"))
        }

        param_df = rbind(param_df, c("num_spots",
                                     num_spots()))

        param_df = rbind(param_df, c("num_simulated_datasets",
                                     num_simulated_datasets()))

        param_df = rbind(param_df, c("parent_simulation_seed",
                                     parent_simulation_seed()))

        param_df = rbind(param_df, c("simulation_seed_for_each_dataset",
                                     simulation_seed_for_each_dataset()))

        param_df = rbind(param_df, c("path_to_output_dir",
                                     path_to_output_dir()))

        param_df = rbind(param_df, c("output_name",
                                     output_name()))



        write.table(param_df,
                    file = con,
                    quote = FALSE,
                    row.names = FALSE,
                    sep = "\t")
      }
    ) # ends export parameters file

    # Run simulation

    parameter_file <- shiny::reactiveVal()
    shiny::observeEvent(input$inputparamfile, {
      parameter_file(input$inputparamfile)
    })

    create_Giotto <- shiny::reactiveVal()
    giotto_folder <- shiny::reactiveVal()

    shiny::observeEvent(input$createGiotto, {
      create_Giotto(input$createGiotto)

      if(input$createGiotto) {
        output$createGiottoSelection <- renderUI({
          shiny::textInput(inputId = "giottofolder",
                           label = "Specify the folder name to export the Giotto object",
                           value = "giottobject",
                           width = "100%")
        })

        shiny::observeEvent(input$giottofolder, {
          giotto_folder(input$giottofolder)
        })
      }


      shiny::observeEvent(input$runsimulation, {

        ParaSimulation(input = parameter_file())

        if(create_Giotto()) {

          x_param = ParaDigest(parameter_file())

          x_expression = read.delim(paste0(x_param$path_to_output_dir,
                                           x_param$output_name,
                                           "_count_1.tsv"),
                                    row.names = 1)

          x_meta = read.delim(paste0(x_param$path_to_output_dir,
                                     x_param$output_name,
                                     "_meta_1.tsv"))

          x_spatlocs = x_meta[,c("Cell", "x.loc", "y.loc")]
          colnames(x_spatlocs) = c("cell_ID", "sdimx", "sdimy")

          x_meta = x_meta[,c("Cell", "annotation", "region")]
          colnames(x_meta)[1] = "cell_ID"

          x = Giotto::createGiottoObject(expression = x_expression,
                                         spatial_locs = x_spatlocs)

          x = Giotto::addCellMetadata(x,
                                      new_metadata = x_meta)

          Giotto::saveGiotto(x, foldername = giotto_folder(), overwrite = TRUE)
        }

        shiny::stopApp()
      })

    }) # ends run simulation



    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
  }

  shiny::runGadget(ui, server)
}

