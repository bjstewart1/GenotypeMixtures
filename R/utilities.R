#########################
# Miscellaneous functions
#########################

#' @title cluster_color_ramp
#' @param cols vector of colours.
#' @import grDevices
#' 
cluster_color_ramp <- function(cols = c("aquamarine", "orange", "blue", "grey", "red", "pink", "deepskyblue", "yellow", "coral", "darkseagreen1"), ...) {
    pal <- colorRampPalette(colors = cols, ...)
    return(pal)
}
