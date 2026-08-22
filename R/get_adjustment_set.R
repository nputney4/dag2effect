get_adjustment_set <- function(cg, exposure, outcome) {
  list(
    backdoor = caugi::adjustment_set(cg, X = exposure, Y = outcome, type = "backdoor"),
    parent = caugi::adjustment_set(cg, X = exposure, Y = outcome, type = "parents"),
    optimal = caugi::adjustment_set(cg, X = exposure, Y = outcome, type = "optimal"),
    ancestor = caugi::ancestors(cg, exposure)
  )
}

