# Apply the scientific timestep to the shared forest. The core advances the clock.
biomass_step <- function(simulation, end_year) {
  timestep <- as.integer(simulation$project$canonical$timestep)
  deferred <- landisR::landis_deferred_regeneration(simulation)
  if (deferred && is.null(simulation$ctx$ward_v7_groups))
    simulation$ctx$ward_v7_groups <- bs52_ward531_groups(simulation$ctx$config, native_order = TRUE)

  checkpoint <-
    .biomass_core$services$checkpoint(
      simulation$year,
      simulation$states,
      simulation$state_eco,
      simulation$cell_state,
      simulation$rng_state,
      simulation$session_signature
    )

  result <-
    bs52_run_large_landscape52(
      ctx =
        simulation$ctx,
      duration =
        end_year,
      seed =
        simulation$seed,
      engine =
        "Rcpp",
      dispersal =
        simulation$dispersal,
      ward_variant =
        simulation$ward_variant,
      output_years =
        sort(
          unique(
            c(
              0L,
              end_year
            )
          )
        ),
      checkpoint_every =
        timestep,
      output_directory =
        NULL,
      run_signature =
        simulation$session_signature,
      resume_checkpoint =
        checkpoint,
      collect_snapshots =
        TRUE,
      rng_chunk_cells =
        simulation$rng_chunk_cells,
      snapshot_compression =
        "gzip1",
      checkpoint_compression =
        "gzip1",
      establishment = simulation$establishment,
      landis_rng = simulation$landis_rng,
      deferred = deferred,
      pending_regeneration = simulation$pending_regeneration
    )

  if (deferred) simulation$pending_regeneration <- NULL
  if (!is.null(simulation$establishment$planting)) {
    simulation$establishment$planting <- NULL
    if (!length(simulation$establishment$excluded_cells))
      simulation$establishment <- NULL
  }
  simulation$states <-
    result$states

  simulation$state_eco <-
    as.integer(
      result$state_eco
    )

  simulation$cell_state <-
    as.integer(
      result$cell_state
    )

  simulation$rng_state <-
    result$final_rng_state
  simulation$landis_rng <- result$landis_rng

  if (nrow(
        result$profile
      )) {

    simulation$profile <-
      if (!nrow(
            simulation$profile
          )) {
        result$profile
      } else {
        rbind(
          simulation$profile,
          result$profile
        )
      }
  }

  simulation
}
