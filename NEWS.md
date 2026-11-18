# landisR.succession.biomass 0.1.0.9008

- Read the LANDIS-II Biomass Succession configuration independently of the landscape.

# landisR.succession.biomass 0.1.0.9005

- Report growth, regeneration and output writing in the console progress bar
  for `landisR::landis_run()`.

# landisR.succession.biomass 0.1.0.9004

- Index explicit initial cohorts by MapCode and reuse their ordered groups.
- Provide a versioned binary canonical cache with input checksums.

# landisR.succession.biomass 0.1.0.9002

- Keep new age-one cohorts with their existing species group in the CSR store.
  When a disturbance leaves tied-age cohorts of different species, this preserves
  the species encounter order used by LANDIS-II during biomass competition.
