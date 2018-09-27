#'AccMatrix class that is used to store accessibility matrices

setClass(
  "AccMatrix",
  slots = list(
    coord = "data.frame",
    accmat = "data.frame"
  )
)
