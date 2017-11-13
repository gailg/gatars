#' @title Display a correlation type matrix
#' 
#' @description  
#' \code{matrix_image_fn} uses \code{geom_tile} from \code{ggplot2} to display a 
#' correlation type image in the following orientation: x increases from left to right
#' and y increases from top to bottom.
#' expect.
#' @param the_matrix The matrix you would like to display.
#' @param main The title if you would like one.
#' @param legend_title The title of the legend. A character string.
#' @param horizontal_line_after If you would like to delineate the rows with a drawn
#' line, you can specify an integer (or integers) an a line (or lines) will be drawn 
#' below that row
#' @param vertical_line_after Analogous to horizontal_line_after but for rows
#' instead of columns
#' @param size The size of the \code{vertical_line_after} and \code{horizontal_line_after}.
#' A positive real numer.
#' @param linetype A \code{ggplot2} linetype.  Examples are integers \code{1}, \code{2}, \code{3}, or 
#' \code{"solid"}, \code{"dashed"}, \code{"dotted"}.
#' @param xlab A character string specifying a label for the x-axis.
#' @param ylab Analogous to xlab but for the y-axis.
#' @param rowticks If you would like a different set of x-axis ticks to be labeled besides the 
#' default ticks of \code{1:nrow(the_matrix)} specify a vector of integers.
#' @param column_ticks Anagloous to \code{rowticks} but for columns.
#' @param legend_position Specifies the position for the legend.  This is the \code{legend_position}
#' from \code{ggplot2}.  Examples are \code{"none"}, \code{"left"}, \code{"right"}, \code{"bottom"},
#' \code{"top"}
#' @examples
#' set.seed(1)
#' NNN = 30
#' epsilon = .8
#' x = rnorm(NNN, 0, 1) 
#' x1 = x + rnorm(NNN, 0, epsilon)
#' x2 = x + rnorm(NNN, 0, epsilon)
#' x3 = x + rnorm(NNN, 0, epsilon)
#' independent_columns = matrix(rnorm(5 * NNN, 0, 1), ncol = 5)
#' data = cbind(x1, x2, x3, independent_columns)
#' matrix_image_fn(cor(data), "The correlation of some multivariate normal data")
#' matrix_image_fn(cor(data), 
#' main = "Showing other options",
#' legend_title = "larry's legend",
#' horizontal_line_after = 3,
#' vertical_line_after = c(3, 5),
#' size = 1, linetype = "dotdash",
#' xlab = "George", ylab = "Harry",
#' row_ticks = 1:3, column_ticks = 1:3,
#' legend_position = "top")
#' @import ggplot2
#' @export
matrix_image_fn = function(
  the_matrix,
  main = "", legend_title = "",
  horizontal_line_after = NULL, 
  vertical_line_after = NULL, 
  size = 3, linetype = 3,
  xlab = "", ylab = "", 
  row_ticks = 1:nrow(the_matrix), 
  column_ticks = 1:ncol(the_matrix),
  legend_position = "right"
){
  III = nrow(the_matrix)
  JJJ = ncol(the_matrix)
  zzz = the_matrix[III:1, ]
  x = 1:III
  y = 1:JJJ
  df = expand.grid(x = x, y = y)
  z = apply(df, 1, function(row){
    zzz[row["x"], row["y"]]
  })
  df$z = z
  image = ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradient(name = legend_title, low = "white", high = "black") +
    xlab(xlab) +
    ylab(ylab) +
    coord_flip() +
    theme(legend.position = legend_position) + 
    ggtitle(main)    
  if(!is.null(row_ticks)){
    image = image + scale_x_continuous(breaks = (III:1)[row_ticks], labels = row_ticks)
  }
  if(!is.null(column_ticks)){
    image = image + scale_y_continuous(breaks = column_ticks, labels = column_ticks)
  }
  if(!is.null(vertical_line_after)){
    image = image + geom_hline(
      yintercept = vertical_line_after + 0.5, 
      color = "red",
      linetype = linetype, 
      size = size)
  }
  if(!is.null(horizontal_line_after)){
    image = image + geom_vline(
      xintercept = (III:1)[horizontal_line_after] - 0.5, 
      color = "red",
      linetype = linetype, 
      size = size)
  }
  image
}

  
