# ============================================================================
# Installed-package runtime bridge.
# This file intentionally overrides development-only sourceCpp loaders.
# ============================================================================

bs52_package_cpp_loaded52 <- function(functions) {
  ns <- environment(bs52_package_cpp_loaded52)
  all(vapply(
    functions,
    function(name) {
      exists(name, envir = ns, mode = "function", inherits = TRUE)
    },
    logical(1)
  ))
}

bs52_init_rcpp <- function(cache_dir = NULL) {
  required_functions <- c(
    "bs52_grow_live_block_cpp",
    "bs52_grow_live_block_dead_cpp",
    "bs52_grow_live_block_dead_pools_cpp",
    "bs52_ward_failure_log_cpp",
    "bs52_recruit_scan_cpp",
    "bs52_shade_light_cpp",
    "bs52_append_recruits_cpp",
    "bs52_recruit_records_cpp",
    "bs52_record_landscape_cpp"
  )

  if (!bs52_package_cpp_loaded52(required_functions)) {
    ns <- environment(bs52_init_rcpp)
    missing <- required_functions[!vapply(
      required_functions,
      function(name) {
        exists(name, envir = ns, mode = "function", inherits = TRUE)
      },
      logical(1)
    )]
    stop(
      "Exports C++ de Biomass Succession absents : ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

bs52_large_init_bulk52 <- function() {
  required_functions <- c(
    "bs52_grow_states_bulk_cpp",
    "bs52_grow_states_bulk_indexed_cpp",
    "bs52_grow_states_bulk_indexed_features_cpp",
    "bs52_state_features_bulk_cpp",
    "bs52_rebuild_states_bulk_cpp",
    "bs52_state_snapshot_bulk_cpp",
    "bs52_grow_csr_indexed_features_cpp",
    "bs52_rebuild_csr_cpp",
    "bs52_ward_probability_fused_bitset_cpp",
    "bs52_ward_probability_fused_bitset_v310_reference_cpp",
    "bs52_growth_kernel_version_cpp"
  )

  ns <- environment(bs52_large_init_bulk52)
  loaded <- vapply(
    required_functions,
    function(name) {
      exists(name, envir = ns, mode = "function", inherits = TRUE)
    },
    logical(1)
  )

  if (!all(loaded)) {
    bs52_large_stop(
      "Exports C++ du package manquants : ",
      paste(required_functions[!loaded], collapse = ", ")
    )
  }

  kernel_version <- bs52_growth_kernel_version_cpp()

  if (!identical(kernel_version, "v3.11-ward-interior-fastpath")) {
    bs52_large_stop(
      "Le kernel C++ charge est ",
      kernel_version,
      " au lieu de v3.11-ward-interior-fastpath."
    )
  }

  invisible(TRUE)
}

bs52_large_code_fingerprint52 <- function() {
  c(
    succession_version = as.character(utils::packageVersion("landisR.succession.biomass")),
    core_version = as.character(utils::packageVersion("landisR")),
    kernel = bs52_growth_kernel_version_cpp()
  )
}

bs52_dead_unsupported52 <- function(...) {
  stop(
    "Dead biomass pools are outside the validated landisR beta.2 scope.",
    call. = FALSE
  )
}

bs52_dead_death <- bs52_dead_unsupported52
bs52_dead_end_year <- bs52_dead_unsupported52
bs52_dead_growth <- bs52_dead_unsupported52
bs52_dead_landscape_snapshot <- bs52_dead_unsupported52
bs52_dead_landscape_spinup <- bs52_dead_unsupported52
bs52_dead_landscape_zero <- bs52_dead_unsupported52
bs52_dead_snapshot <- bs52_dead_unsupported52
bs52_dead_state <- bs52_dead_unsupported52
