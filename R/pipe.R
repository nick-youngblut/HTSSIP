#' Pipe
#'
#' Like dplyr, HTSSIP also uses the pipe function, \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @examples
#' data(physeq_S2D2_l)
#' physeq = physeq_S2D2_l[[1]]
#' # Instead of:
#' df_l2fc = HRSIP(physeq, design=~Substrate)
#' head(df_l2fc)
#' # You can write
#' df_l2fc = physeq %>% HRSIP(design=~Substrate)
#' head(df_l2fc)
NULL
