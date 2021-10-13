# This is a Shiny web application. You can run the application by clicking

library(shiny)
library(shinyjs)
library(jpeg)
library(ggplot2)
library(plotly)
library(DT)

# Source functions used
source("functions.R")

# Metadata dataframe set as NULL at the beginning to avoid showing error
if (! base::exists("metadata_df")) metadata_df <- NULL
if (! base::exists("img")) img_obj <- NULL
# metadata_df <- base::readRDS("../tonsil_atlas/spatial_transcriptomics/misc/metadata_esvq52_nluss5.rds")
# img_obj <- jpeg::readJPEG("../tonsil_atlas/spatial_transcriptomics/01-spaceranger/img/esvq52_nluss5_V19S23-039_C1.jpg")

# img <- jpeg::readJPEG(img_path)
# Define UI for application that draws a histogram
ui <- shiny::fluidPage(

    # Application title
    shiny::titlePanel("Pathology annotation"),

    # Sidebar with a slider input for number of bins 
    shiny::sidebarLayout(
        shiny::sidebarPanel(width = 2,
            # Load Metadata
            shiny::fileInput(inputId = "metadata",
                      label = "Metadata + 2D embedings (coord_x, coord_y)",
                      multiple = FALSE),
            
            # Load Expression data
            shiny::fileInput(inputId = "img_fn",
                      label = "Image",
                      multiple = FALSE),
            
            # Enter value to group by
            shiny::selectizeInput("groupby", "Coloring feature:",
                           selected = NULL,
                           choices = NULL,
                           options = list(create = TRUE),
                           # multiple = TRUE,
                           multiple = FALSE),
            shiny::actionButton(inputId = "apply_groupby", label = "Update Grouping"),
            
            # Slider to determine point size of UMAP plots
            # sliderInput("size", "Dot size:",
            #             min = 0,
            #             max = 5,
            #             value = 3,
            #             step = 0.05),
            # Slider to determine point size of UMAP plots
            # sliderInput("alpha", "Dot transparency:",
            #             min = 0,
            #             max = 1,
            #             value = 1,
            #             step = 0.1)
        ),

        # Show a plot of the generated distribution
        shiny::mainPanel(
            # fluidRow(
                     # column(6,shiny::plotOutput('img_plot')),
                     # column(6,plotly::plotlyOutput('sel_plot'))
            # ),
            # shiny::plotOutput("img_plot", width = "50%"),
            plotly::plotlyOutput("sel_plot", width = "100%", height = 800),
            DT::DTOutput("barcode_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # Setting maximum file size to 8GB
    options(shiny.maxRequestSize = 8000 * 1024 ^ 2)
    
    observeEvent(input$apply_checkbox, {
        shinyjs::toggle("interactive_filter")
    })
    
    ##############################
    ##### Update input fields ####
    ##############################
    shiny::observe({
        
        # Load data
        # Read marker list
        file1 <- input$metadata
        if (base::is.null( file1 ) ) { return() }
        metadata_df <<- base::readRDS(file1$datapath)
        
        file2 <- input$img_fn
        if( base::is.null( file2 ) ) { return() }
        # Load image as jpg
        img_obj <<- jpeg::readJPEG(file2$datapath)
        
        # Subset character/factor columns with <= 50 unique values
        feat_ch1 <- base::sapply(metadata_df, function(x) base::length(base::unique(x)) <= 50)
        feat_ch2 <- base::sapply(base::colnames(metadata_df), function(x) !base::is.numeric(metadata_df[, x]))
        feat_ch <- base::colnames(metadata_df)[feat_ch1 & feat_ch2]
        
        # Subset numeric metadata columns
        feat_num1 <- base::sapply(base::colnames(metadata_df), function(x) base::is.numeric(metadata_df[, x]))
        feat_num <- base::colnames(metadata_df)[feat_num1]
        
        # Join metadata variables of interest subset
        feat_sub <- c(feat_ch, feat_num)
        
        # Update groupby selection
        shiny::updateSelectizeInput(session,
                             inputId = "groupby",
                             choices = c("", 
                                         feat_sub),
                             selected = feat_sub[1],
                             # consider using server-side selectize for massively improved performance. See the Details section of the ?selectizeInput help topic.
                             # The server-side selectize input uses R to process searching, and R will return the filtered data to selectize.
                             server = TRUE)
    })
    
    ########################################
    ######## Setting reactive events #######
    ########################################
    groupby_var <- shiny::eventReactive(input$apply_groupby, {
        input$groupby
    })
    
    dfInput <- shiny::reactive({
        ## modifying metadata to see if this helps
        metadata_df[base::order(metadata_df[, groupby_var()]), ]
    })
    
    # Lasso selection plot
    output$sel_plot <- plotly::renderPlotly({
        
        # Read data from reactive observed slots
        metadata_df <- dfInput()
        # if (!bases::is.null(metadata_df)) {
        pp <- plotly::plot_ly(
            data = metadata_df,
            x = metadata_df[, "coord_x"],
            y = metadata_df[, "coord_y"],
            # color = metadata_df[, groupby_var()],
            color = metadata_df[, "Spatial_snn_res.0.8"],
            size = 0.75,
            alpha = 0.5
        ) %>%
        plotly::layout(
            dragmode = "lasso",
            yaxis = base::list(
                scaleanchor = "x",
                scaleratio = 1),
            images = base::list(
                # source = base64enc::dataURI(file = "../tonsil_atlas/spatial_transcriptomics/01-spaceranger/img/esvq52_nluss5_V19S23-039_C1.jpg"),
                source = plotly::raster2uri(grDevices::as.raster(img_obj)),
                x = 0, y = -base::nrow(img_obj),
                sizex = base::ncol(img_obj), sizey = -base::nrow(img_obj),
                xref = "x", yref = "y",
                xanchor = "left", yanchor = "bottom",
                # https://plotly.com/r/reference/#Layout_and_layout_style_objects
                sizing = "contain",
                layer = "below"
            )
        )
            # } else pp <- NULL
            
            pp
        
        })
    
    # Lasso selection table
    output$barcode_table <- DT::renderDT(
        # It is important here to set server = FALSE so that when we save the table in CSV format i saves ALL the entries. By default it is TRUE and it only saves those entries currently shown in the table!
        # Note that your users might run into performance and memory issues using server=FALSE if your data table is very large.
        # https://stackoverflow.com/questions/50508854/button-extension-to-download-all-data-or-only-visible-data
        server = FALSE,
        {
            # Read data from reactive observed slots
            metadata_df <- dfInput()
            
            metadata_df <- metadata_df %>%
                # Round and make as character so we can then join with the d dataframe representing the coordinates in the shiny app.
                dplyr::mutate(
                    coord_x = base::as.character(base::round(coord_x, 6)),
                    coord_y = base::as.character(base::round(coord_y, 6))
                )
            
            
            d <- plotly::event_data(event = "plotly_selected")
            if (!base::is.null(d)) {
                d <- d %>%
                    # we need to round the numbers to 6 since those are the coord that the metadata gives, shiny app returns up tp 14 decimals so the left join doesn't work
                    dplyr::mutate(
                        x = base::as.character(base::round(x, 6)),
                        y = base::as.character(base::round(y, 6))
                    ) %>%
                    # dplyr::mutate(y = -y) %>%
                    dplyr::left_join(metadata_df,
                                     by = c("y" = "coord_y",
                                            "x" = "coord_x")) %>%
                    dplyr::select("pointNumber", "x", "y", "barcode")
                
                # filename <- glue::glue("{sel_gene()}-shinyapp")
                filename <- glue::glue("selection-shinyapp")
                # Return datatable with the csv option to save the table directly
                # Download options following this -> https://rstudio.github.io/DT/003-tabletools-buttons.html
                # Also this -> https://github.com/rstudio/DT/issues/409
                DT::datatable(
                    data = d,
                    extensions = c("Buttons"),
                    options = base::list(
                        dom = 'Bfrtip',
                        buttons =  base::list(base::list(extend = "csv",
                                                         filename = filename))
                        # buttons = list(list(
                        #   extend = 'collection',
                        #   buttons = c('csv', 'excel', 'pdf'),
                        #   text = 'Download'
                        # ))
                        # list("copy", "print", list(
                        #   extend = "collection",
                        #   buttons = c("csv", "excel", "pdf"),
                        #   filename = filename,
                        #   text = "Download"
                        # ))
                    )
                )
            }
        })}

# Run the application 
shiny::shinyApp(ui = ui, server = server)
