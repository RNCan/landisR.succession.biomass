.biomass_core <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .biomass_core$services <- landisR::landis_succession_services()
}

biomass_context <- function(canonical, v7 = FALSE) {
  ctx <- bs52_large_context52(canonical)
  if (v7) ctx <- bs52_ward531_context(ctx)
  ctx
}

landis_biomass_backend <- function() {
  list(
    api_version = 2L,
    context = biomass_context,
    step = biomass_step,
    fire = list(capture_sources = landisR_serotiny_prefire_sources,
                recruit = landisR_serotiny_recruit_after_fire),
    normalize_fire_traits = landisR_normalize_serotiny_traits,
    read_canonical = bs52_read_canonical52,
    import_landis = bs52_import_landis_to_r,
    run_production = bs52_run_large_production52
  )
}
