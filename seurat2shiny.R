#' This function takes in a Seurat 3.0 object and returns a named list with 2
#' objects formated to be loaded in the ShinyApp:
#'
#'      1. Metadata + coordinates of the 2D embedding.
#'      2. Expression matrix of the selected slot.
#'
#' ShinyApp: https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/shiny-annotation/blob/main/seurat2shiny.R
#'
#' @param object Object of class Seurat. Mandatory.
#' @param tech Character string. Specify the technology 
#' @param assay Character string. Assay within the Seurat object from which to extract the expression matrix. Default: active assay.
#' @param slot Character string. Slot containing the expression matrix. Default: data.
#' @param reduction Character string. Dimensionality reduction from which to extract the 2D coordinates. Default: umap.
#' @param image Character string or NULL. When tech is sp, name of the image from which we want to extract the coordinates as found in names(object/@images), by default NULL.
#' @param asfactors Character vector. Metadata columns to be converted to factors. Default: NULL.
#' @param save Logical value. Save metadata and expression matrix objects as RDS files. Default: FALSE
#' @param path Character string. Path to save output files if 'save = TRUE'. Default: working directory.
#'
#' @return Named list containing the joint metadata + 2D embedding and the expression matrix.
#'
#' @examples
#' seurat2shiny( object = seurat_object, asfactors = c("plate", "replicate") )
#'
#' shiny_list = seurat2shiny(object = seurat_object)
#'
#' @export

seurat2shiny = function(
  object                         ,
  tech      = c("sc", "sp")      ,
  assay     = object@active.assay,
  slot      = "data"             ,
  reduction = "umap"             ,
  image     = NULL               ,
  asfactors = NULL               ,
  save      = FALSE               ,
  path      = "."                  # path = getwd()
) {
  suppressMessages( library(Seurat) );
  
  # Input check.
  if ( ! is(object, "Seurat") )
    stop("'object' is not a Seurat object.");
  
  if ( ! assay %in% Seurat::Assays(object) )
    stop("'assay' not in the Seurat object's available assays.");
  
  if ( tech == "sc" & ! (reduction %in% names(object@reductions)) )
    stop("'reduction' not in the Seurat object's available reductions.");
  
  if ( ! slot %in% c("counts", "data", "scale.data") )
    stop("'slot' not in the Seurat object's available slots.");
  
  if ( ! tech %in% c("sc", "sp") )
    stop("tech must be sc or sp.");
  
  
  # Check Which technology it is processing
  if (tech == "sc") {
    # Extract 2D coordinates. Keep only first 2 dimensions, remove the rest if any.
    embeds <- as.data.frame(object@reductions[[reduction]]@cell.embeddings)[1:2];
    names(embeds) <- c("coord_x", "coord_y");
  } else if (tech == "sp") {
    # If the image is null select the first one
    if (is.null(image)) {
      image <- names(object@images)[1]
      warning(sprintf("image is not set, we will use %s", image))
    } 
    
    embeds <- data.frame(object@images[[image]]@coordinates[, c("imagerow", "imagecol")])
    colnames(embeds) <- c("coord_y", "coord_x");
    
    # Inverse coord_y
    embeds$coord_y <- - embeds$coord_y
  }
  
  
  # Join metadata with coordinates.
  metadata <- object@meta.data;
  
  for (col in asfactors) {
    metadata[[col]] <- as.factor(metadata[[col]]);
  };
  
  metadata <- merge(x = metadata, y = embeds, by = "row.names");
  names(metadata)[1] <-  "barcode"; # names(metadata)[names(metadata) == "Row.names"] = "barcode";
  rownames(metadata) <- metadata$barcode
  metadata$barcode <- as.character(metadata$barcode)
  # Reset barcode order which is shuffled in merge
  metadata <- metadata[rownames(object@meta.data), ]
  
  # Extract expression data.
  # expression = as.matrix( Seurat::GetAssayData(object = object, slot = slot, assay = assay) );
  expression = Seurat::GetAssayData(object = object, slot = slot, assay = assay);
  
  if ( ! identical( as.character(metadata$barcode), colnames(expression) ) )
    warning("Cells in metadata and expression matrix do not match.");
  
  if (save) {
    saveRDS( object = metadata  , file = paste0(path, "/metadata.rds"  ) );
    saveRDS( object = expression, file = paste0(path, "/expression.rds") );
  };
  
  invisible(
    list(metadata = metadata, expression = expression)
  );
}