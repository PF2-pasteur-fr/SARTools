#' Cairo Size Wrapper
#'
#' Limits the pixel size to either to the minimum of the given and the max allowed pixel size 32767x32767 pixels if
#' Cairo graphical has been choosen
#'
#' @param requested pixel size
#' @return the minimum between the requested and maximal allowed pixel size if Cairo is selected otherwise the requested pixel size
#' @author vipul patel




cairoSizeWrapper <- function(pixel_re){
  if(options("bitmapType") == "cairo") return(min(pixel_re,32767)) else return(pixel_re)

}