# This is a Shiny web application. You can run the application by clicking

library(shiny)
library(shinyjs)
library(shinythemes)
library(grDevices)
library(jpeg)
library(png)
library(plotly)
library(dplyr)
library(DT)

# Source functions used
source("functions.R")

# Metadata dataframe set as NULL at the beginning to avoid showing error
if (! base::exists("metadata_df")) metadata_df <- NULL
if (! base::exists("img")) img_obj <- NULL

# Define UI for application that draws a histogram
ui <- shiny::fluidPage(
    # Add a theme
    theme = shinythemes::shinytheme("flatly"),

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
            shiny::fileInput(inputId = "expr",
                             label = "Expression matrix",
                             multiple = FALSE),
            
            # Load Image Jpeg
            shiny::fileInput(inputId = "img_fn",
                      label = "Image",
                      multiple = FALSE),
            
            # Enter value to group by
            shiny::selectizeInput("groupby", "Coloring feature:",
                           selected = NULL,
                           choices = NULL,
                           options = base::list(create = TRUE),
                           # multiple = TRUE,
                           multiple = FALSE),
            shiny::actionButton(inputId = "apply_groupby",
                                label = "Update Grouping"),
            
            # Select genes to use as marker
            shiny::selectizeInput("gene_mrkr", "Select marker gene:",
                           selected = NULL,
                           choices = NULL,
                           # Don't allow to create new gene names
                           options = base::list(create = FALSE),
                           multiple = FALSE),
            shiny::actionButton(inputId = "apply_marker",
                                label = "Update marker"),
            
            # Radio button to choose to show gene or metadata feature
            shiny::radioButtons("radio", "Choose coloring:",
                                c("Group" = "grp",
                                  "Gene" = "gene"))
            
            ),
        # Show a plot of the generated distribution
        shiny::mainPanel(
            plotly::plotlyOutput("sel_plot", width = "100%", height = 800),
            DT::DTOutput("barcode_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # Setting maximum file size to 8GB
    base::options(shiny.maxRequestSize = 8000 * 1024 ^ 2)
    
    shiny::observeEvent(input$apply_checkbox, {
        shinyjs::toggle("interactive_filter")
    })
    
    ##############################
    ##### Update input fields ####
    ##############################
    shiny::observe({
        
        # Load metadata
        file1 <- input$metadata
        if (base::is.null(file1)) return()
        metadata_df <<- base::readRDS(file1$datapath)
        
        # Load jpeg image
        file2 <- input$img_fn
        if (base::is.null(file2)) return()
        # Load image as jpg
        img_obj <<- jpeg::readJPEG(file2$datapath)
        
        # Load expression matrix
        file3 <- input$expr
        if(base::is.null(file3)) return()
        # Keep the count matrix sparse
        expr_mtrx <<- base::readRDS(file3$datapath) %>%
            Matrix::Matrix(data = ., sparse = TRUE)
        
        # Update filter gene
        shiny::updateSelectizeInput(
            session,
            inputId = "gene_mrkr",
            choices = base::rownames(expr_mtrx),
            selected = base::rownames(expr_mtrx)[1],
            # consider using server-side selectize for massively improved performance. See the Details section of the ?selectizeInput help topic.
            # The server-side selectize input uses R to process searching, and R will return the filtered data to selectize.
            server = TRUE)
        
        # Subset character/factor columns with <= 50 unique values
        feat_ch1 <- base::sapply(metadata_df, function(x) {
            base::length(base::unique(x)) <= 50
            })
        # Keep numeric columns
        feat_num <- base::sapply(base::colnames(metadata_df), function(x) {
            base::is.numeric(metadata_df[, x])
            })
        # Remove coord_x, coord_y
        feat_num <- feat_num[!feat_num %in% c("coord_x", "coord_y")]
        
        # Get features to show for the metadata either character with less than 50 categories or numeric
        feats <- base::colnames(metadata_df)[feat_ch1 | feat_num]
        
        # Update groupby selection
        shiny::updateSelectizeInput(
            session,
            inputId = "groupby",
            choices = c("", 
                        feats),
            selected = feats[1],
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
    
    gene_viz <- shiny::eventReactive(input$apply_marker, {
        input$gene_mrkr
    })
    
    df_input <- shiny::reactive({
        ## modifying metadata to see if this helps
        ## Need to add this groupby_var() ordering so the app can read the metata file, need to check why
        metadata_df[base::order(metadata_df[, groupby_var()]), ]
        metadata_df[base::colnames(expr_mtrx), ]
    })
    
    # Lasso selection plot
    output$sel_plot <- plotly::renderPlotly({
        
        # Read data from reactive observed slots
        metadata_df <- df_input()
        
        # Set the color vector
        color_vr <- if (input$radio == "grp") {
            col_vec <- metadata_df[, groupby_var()]
        } else if (input$radio == "gene") {
            col_vec <- expr_mtrx[gene_viz() , ]
            # Below doesn't work well with selectize
            # col_vec <- expr_mtrx[input$gene_mrkr, ]
        }
        
        pp <- plotly::plot_ly(
            data = metadata_df,
            x = metadata_df[, "coord_x"],
            y = metadata_df[, "coord_y"],
            # color = metadata_df[, groupby_var()],
            color = col_vec,
            size = 0.4,
            alpha = 0.9
        ) %>%
        # https://plotly-r.com/embedding-images.html
        plotly::layout(
            dragmode = "lasso",
            yaxis = base::list(
                scaleanchor = "x",
                scaleratio = 1),
            images = base::list(
                source = plotly::raster2uri(grDevices::as.raster(img_obj)),
                x = 0,
                y = -base::nrow(img_obj),
                sizex = base::ncol(img_obj),
                sizey = -base::nrow(img_obj),
                xref = "x",
                yref = "y",
                xanchor = "left",
                yanchor = "bottom",
                # https://plotly.com/r/reference/#Layout_and_layout_style_objects
                sizing = "contain",
                layer = "below"
            )
        )
            pp
        
        })
    
    # Lasso selection table
    output$barcode_table <- DT::renderDT(
        # It is important here to set server = FALSE so that when we save the table in CSV format i saves ALL the entries. By default it is TRUE and it only saves those entries currently shown in the table!
        # Note that your users might run into performance and memory issues using server=FALSE if your data table is very large.
        # https://stackoverflow.com/questions/50508854/button-extension-to-download-all-data-or-only-visible-data
        server = FALSE, {
            # Read data from reactive observed slots
            metadata_df <- df_input()
            
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
                    dplyr::left_join(metadata_df,
                                     by = c("y" = "coord_y",
                                            "x" = "coord_x")) %>%
                    dplyr::select("pointNumber", "x", "y", "barcode")
                
                filename <- "selection-shinyapp"
                # Return datatable with the csv option to save the table directly
                # Download options following this -> https://rstudio.github.io/DT/003-tabletools-buttons.html
                # Also this -> https://github.com/rstudio/DT/issues/409
                DT::datatable(
                    data = d,
                    extensions = c("Buttons"),
                    options = base::list(
                        dom = "Bfrtip",
                        buttons =  base::list(base::list(extend = "csv",
                                                         filename = filename))
                    )
                )
            }
        })}

# Run the application 
shiny::shinyApp(ui = ui, server = server)
