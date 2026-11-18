\name{biomass_succession_read}
\alias{biomass_succession_read}
\title{Read LANDIS-II Biomass Succession inputs}
\usage{
biomass_succession_read(source_directory, scenario = "scenario.txt",
  prefer_cells_csv = TRUE)
}
\arguments{
  \item{source_directory}{Base directory for the LANDIS-II scenario.}
  \item{scenario}{Scenario path, relative to \code{source_directory} or absolute.}
  \item{prefer_cells_csv}{Use \file{cells.csv} from the scenario directory
  when available. Its \code{cell_id} column identifies cells and
  \code{community} gives the initial community code. With \code{FALSE},
  use the initial communities raster.}
}
\description{
Read the Biomass Succession configuration named in the scenario. It includes
the timestep, species and ecoregion parameters, dynamic inputs, light and
shade tables, and initial communities and their cell assignments.
}
\value{
A \code{landisR_biomass_succession} configuration to pass as
\code{succession} to \code{landisR::landis_initialize()} with a shared
landscape. The initializer checks the species catalog and the community
assignments against the landscape, then creates the initial cohorts.
}
\examples{
\dontrun{
landscape <- landisR::landis_read("inputs")
succession <- biomass_succession_read("inputs")
simulation <- landisR::landis_initialize(landscape, succession = succession)
}
}
