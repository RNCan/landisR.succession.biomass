\name{landis_biomass_backend}
\alias{landis_biomass_backend}
\title{Biomass Succession provider for landisR}
\usage{landis_biomass_backend()}
\description{
Provides the engine used by biomass projects. Use \code{landisR::landis_build()},
\code{landisR::landis_sim()}, or \code{landisR::landis_run()} to run the model.
}
\details{
The provider exposes its fire response through \code{fire$capture_sources}
and \code{fire$recruit}. The first callback is called before the core
applies cohort mortality. It receives the current simulation and affected
cell IDs and returns the sources to carry across the mortality edit.
The second callback receives the post-mortality state, those sources, the
stochastic setting, and the event year. It returns the updated simulation
and source and recruitment counts for the fire journal. The core keeps
both calls in the fire event, before the next extension runs.
}
\value{
A versioned list of functions used by the landisR core to initialize and
advance biomass projects, process fire-driven recruitment, and run
production simulations.
}
\keyword{internal}
