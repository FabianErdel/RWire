#'PEPairs class that is used to store sets of promoter-enhancer pairs

setClass(
  "PEPairs",
  slots = list(
    promoters = "data.frame",
    enhancers = "data.frame"
  )
)
