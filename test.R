library(ggplot2)

metadata_df <- readRDS("../tonsil_atlas/spatial_transcriptomics/misc/metadata_esvq52_nluss5.rds")
img <- jpeg::readJPEG("../tonsil_atlas/spatial_transcriptomics/01-spaceranger/img/esvq52_nluss5_V19S23-039_C1.jpg")



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
                                 color = Spatial_snn_res.0.8,
                                 text = barcode),
                    size = as.numeric(1),
                    alpha = 1) +
  cowplot::theme_half_open(11, rel_small = 1) +
  ggplot2::theme_void() +
  ggplot2::coord_fixed(ratio = 1,
                       xlim = NULL,
                       ylim = NULL,
                       expand = TRUE,
                       clip = "on")

# Get selected barcodes
ggplotly(p, tooltip = "text") %>%
  layout(dragmode = "lasso")

library(ggplot2)
library(plotly)

df <- data.frame(
  id = 1:100,
  value = rnorm(100)
)

# image_destination <- "https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/tonsil_atlas/blob/spatial/spatial_transcriptomics/01-spaceranger/img/c28w2r_7jne4i_V19S23-039_B1.jpg"
# image_destination <- "https://raw.githubusercontent.com/MarcElosua/test/main/esvq52_nluss5_V19S23-039_C1.jpg"
image_destination <- "https://raw.githubusercontent.com/MarcElosua/test/main/esvq52_nluss5_V19S23-039_C1.jpg?raw=true"
# image_destination <- "https://github.com/MarcElosua/test/blob/main/esvq52_nluss5_V19S23-039_C1.jpg"
# gg <- ggplot(metadata_df, aes(x = coord_x, y = coord_y)) +
#   geom_point()
library(plotly)
plot_ly(
  data = metadata_df,
  x = metadata_df[, "coord_x"],
  y = metadata_df[, "coord_y"],
  color = metadata_df[, "Spatial_snn_res.0.8"],
  size = 1,
  alpha = 0.8
  ) %>%
  layout(
    images = list(
      list(
        source =  image_destination,
        xref = "x",
        yref = "y",
        x = 1,
        y = 3,
        sizex = 2,
        sizey = 2,
        sizing = "stretch",
        opacity = 0.4,
        layer = "below"
      )
    )
)

EBImage::readImage('https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Crab_Nebula.jpg/240px-Crab_Nebula.jpg')
img <- EBImage::readImage("https://raw.githubusercontent.com/MarcElosua/test/main/esvq52_nluss5_V19S23-039_C1.jpg")
fig <- plot_ly(type="image", z=img*250)

# Required packages
install.packages("biOps")
install.packages("raster")
library(biOps)
library(raster)
img <- EBImage::readImage(files = "../tonsil_atlas/spatial_transcriptomics/01-spaceranger/img/esvq52_nluss5_V19S23-039_C1.jpg")
plot_ly(type="image", z=img)
# img = EBImage::readImage('https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Crab_Nebula.jpg/240px-Crab_Nebula.jpg')
# img = EBImage::readImage("https://raw.githubusercontent.com/Single-Cell-Genomics-Group-CNAG-CRG/shiny-annotation/main/esvq52_nluss5_V19S23-039_C1.jpg?token=AD42NAYYD5BR4CK6IPYQLUDBMXK4A")
plotly::layout(pp, images = list(
  list(
    source =  image_destination,
    xref = "x",
    yref = "y",
    x = 0,
    y = 2,
    sizex = 100,
    sizey = 4,
    sizing = "stretch",
    opacity = 0.8,
    layer = "below"
  )
))
  

library(plotly)
p <- plot_ly(x=c(1, 2, 3), y=c(1, 2, 3), type='scatter')
p2 <-
  p %>%
  layout(images = list(
    list(source = "https://github.com/MarcElosua/test/blob/main/esvq52_nluss5_V19S23-039_C1.jpg",
         xref = "paper",
         yref = "paper",
         x= 0,
         y= 1,
         sizex = 0.2,
         sizey = 0.2,
         opacity = 0.8)
  )
  )
p2
