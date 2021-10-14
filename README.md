# shiny-pathology
This repo contains a very simple shiny app to annotate spatial transcriptomics slides.

The user just needs to upload 2 files: \
- metadata file from the Seurat object using the function seurat2shiny that can be found [here](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/shiny-pathology/blob/main/seurat2shiny.R) \
- jpeg image from the visium slide, note that the heavier the image the longer it will take to render the plot \

After uploading the necessary files one just needs to select the coloring feature, wait a couple of seconds and the spatial transcriptomic spots with the underlying image will appear!

Then just use the lasso selection tool as shown below to select regions of interest, the barcode of the spots falling within selected are will appear in a dynamic table below which can be downloaded pressing the *CSV* button

![](img/shiny-pathology.gif)
