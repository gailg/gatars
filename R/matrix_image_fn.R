#' @title Display a correlation type matrix
#' 
#' @description  
#' \code{matrix_image_fn} uses \code{geom_tile} from \code{ggplot2} to display a 
#' correlation type image in the orientation that I think a mathematician would
#' expect.
#' @param the_matrix The matrix you would like to display.
#' @param main The title if you would like one.
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
#' @import ggplot2
#' @export
matrix_image_fn = function(
  the_matrix,
  main = "", legend_name = "",
  vertical_line_after = NULL, horizontal_line_after = NULL, 
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
    scale_fill_gradient(name = legend_name, low = "white", high = "black") +
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

  
