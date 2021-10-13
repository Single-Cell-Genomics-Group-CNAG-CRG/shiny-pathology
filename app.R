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
if (! exists("metadata_df")) metadata_df <- NULL
if (! exists("img")) img <- NULL
metadata_df <- readRDS("../tonsil_atlas/spatial_transcriptomics/misc/metadata_esvq52_nluss5.rds")
img <- jpeg::readJPEG("../tonsil_atlas/spatial_transcriptomics/01-spaceranger/img/esvq52_nluss5_V19S23-039_C1.jpg")

# img <- jpeg::readJPEG(img_path)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Pathology annotation"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width = 2,
            # Load Metadata
            fileInput(inputId = "metadata",
                      label = "Metadata + 2D embedings (coord_x, coord_y)",
                      multiple = FALSE),
            
            # Load Expression data
            fileInput(inputId = "img",
                      label = "Image",
                      multiple = FALSE),
            # Slider to determine point size of UMAP plots
            sliderInput("size", "Dot size:",
                        min = 0,
                        max = 5,
                        value = 3,
                        step = 0.05),
            # Slider to determine point size of UMAP plots
            sliderInput("alpha", "Dot transparency:",
                        min = 0,
                        max = 1,
                        value = 1,
                        step = 0.1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                     column(6,shiny::plotOutput('img_plot')),
                     column(6,plotly::plotlyOutput('sel_plot'))
            ),
            # shiny::plotOutput("img_plot", width = "50%"),
            # plotly::plotlyOutput("sel_plot", width = "50%"),
            DT::DTOutput("barcode_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    ##############################
    ##### Update input fields ####
    ##############################
    # observe({})
        
    # Image overlapping
    output$img_plot <- renderPlot({
        # Read data from reactive observed slots
        # metadata_df <- dfInput()
        
        # Convert image to grob object
        img_grob <- grid::rasterGrob(img,
                                     interpolate = FALSE,
                                     width = grid::unit(1, "npc"),
                                     height = grid::unit(1, "npc"))
        
        p <- ggplot2::ggplot(data = metadata_df) +
            ggplot2::annotation_custom(
                grob = img_grob,
                xmin = 0,
                xmax = ncol(img),
                ymin = 0,
                ymax = -nrow(img)) +
        ggplot2::geom_point(ggplot2::aes(x = coord_x,
                                         y = coord_y,
                                         # color = metadata_df[, groupby_var()],
                                         color = metadata_df[, "Spatial_snn_res.0.8"]),
                            size = 1,
                            alpha = input$alpha) +
        ggplot2::theme_classic()
        
        p
    })
    
    # Lasso selection plot
    output$sel_plot <- renderPlotly({
        
        # Read data from reactive observed slots
        # metadata_df <- dfInput()
        
        p <- ggplot2::ggplot(data = metadata_df) +
            ggplot2::geom_point(ggplot2::aes(x = coord_x,
                                             y = coord_y,
                                             # color = metadata_df[, groupby_var()],
                                             color = metadata_df[, "Spatial_snn_res.0.8"],
                                             text = barcode),
                                size = as.numeric(input$size)) +
            ggplot2::theme_classic()
        
        # Get selected barcodes
        g <- ggplotly(p, tooltip = "text") %>%
            layout(dragmode = "lasso",
                   yaxis = list(
                       scaleanchor = "x",
                       scaleratio = 1))
        
        g
        })
    
    # Lasso selection table
    output$barcode_table <- DT::renderDT(
        # It is important here to set server = FALSE so that when we save the table in CSV format i saves ALL the entries. By default it is TRUE and it only saves those entries currently shown in the table!
        # Note that your users might run into performance and memory issues using server=FALSE if your data table is very large.
        # https://stackoverflow.com/questions/50508854/button-extension-to-download-all-data-or-only-visible-data
        server = FALSE,
        {
            # Read data from reactive observed slots
            # metadata_df <- dfInput()
            metadata_df <- metadata_df %>%
                # Round and make as character so we can then join with the d dataframe representing the coordinates in the shiny app.
                dplyr::mutate(
                    coord_x = as.character(round(coord_x, 6)),
                    coord_y = as.character(round(coord_y, 6))
                )
            
            
            d <- plotly::event_data(event = "plotly_selected")
            if (!is.null(d)) {
                d <- d %>%
                    # we need to round the numbers to 6 since those are the coord that the metadata gives, shiny app returns up tp 14 decimals so the left join doesn't work
                    dplyr::mutate(
                        x = as.character(round(x, 6)),
                        y = as.character(round(y, 6))
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
                    options = list(
                        dom = 'Bfrtip',
                        buttons =  list(list(extend = "csv", filename = filename))
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
shinyApp(ui = ui, server = server)
