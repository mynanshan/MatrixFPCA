library(rTensor)

source("matfdata.R")
source(file.path("experiment", "gendata-simple.R"))

check_visit_range <- function(ratio, expected) {
  generated <- gen_data_simple(
    n = 100L, p = 2L, q = 2L, npc_row = 2L, npc_col = 2L,
    npc_time = 3L, nt = 51L, sampling_ratio = ratio, sig = 0,
    spatial_model = "ar1"
  )
  observed <- lengths(get_t_list(generated$data))
  stopifnot(
    identical(generated$visit_counts, observed),
    min(observed) >= expected[1L],
    max(observed) <= expected[2L]
  )
}

set.seed(20260816)
check_visit_range(c(0.10, 0.20), c(6L, 10L))
check_visit_range(c(0.20, 0.40), c(11L, 20L))

cat("visit-density tests passed\n")
