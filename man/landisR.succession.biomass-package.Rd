\name{landisR.succession.biomass-package}
\alias{landisR.succession.biomass-package}
\alias{landisR.succession.biomass}
\title{landisR.succession.biomass: Living Biomass Succession for LANDIS-R}
\description{
Simulate cohort growth, mortality, seed dispersal and establishment.
The engine also provides initial community preparation and post-fire
recruitment for biomass projects.
}
\details{
The \pkg{landisR} core owns the shared forest state and simulation clock.
Build or import a biomass project with \code{landisR::landis_build()} or
\code{biomass_succession_read()}, then combine them with \code{landisR::landis_initialize()}. Run the initialized simulation with \code{landisR::landis_sim()}
or \code{landisR::landis_run()}. The core loads the biomass engine
and calls its provider automatically.
}
\seealso{\code{\link{landis_biomass_backend}}}
\docType{package}
