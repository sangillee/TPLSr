#' Sample participant data from a left-right button press task
#'
#' A dataset containing five sample participant's binary button presses inside the scanner (left/right).
#'
#' @format A data frame with following variables
#' \describe{
#'   \item{X}{Brain image single trial coefficients. N-by-v matrix}
#'   \item{Y}{Left = 0, Right = 1, binary indicator of participant choice}
#'   \item{subj}{Subject number (i.e., 1, 2, 3)}
#'   \item{run}{Run number (i.e., 1, 2, 3, 4, 5, 6, 7, 8)}
#'   \item{mask}{Binary 3D brain image that indexes where the variables in X came from.}
#' }
#' @source Kable, J. W., Caulfield, M. K., Falcone, M., McConnell, M., Bernardo, L., Parthasarathi, T., ... & Diefenbach, P. (2017). No effect of commercial cognitive training on brain activity, choice behavior, or cognitive performance. Journal of Neuroscience, 37(31), 7390-7402.
"TPLSdat"
