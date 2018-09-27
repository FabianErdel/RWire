#'CorMatrix class that is used to store correlation matrices

setClass(
  "CorMatrix",
  slots = list(
    coord1 = "data.frame",
    coord2 = "data.frame",
    cormat = "data.frame"
  )
)