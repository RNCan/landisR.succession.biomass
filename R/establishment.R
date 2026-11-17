landisR_establishment_mask <- function(mask, establishment, ids, cell_state,
                                       state_eco, species_names, sep) {
  mask[match(establishment$excluded_cells, ids)] <- 0L
  planting <- establishment$planting
  if (!NROW(planting)) return(mask)

  position <- match(planting$cell_id, ids)
  species <- match(planting$species, species_names)
  eligible <- sep[cbind(state_eco[cell_state[position]], species)] > 0
  if (!any(eligible)) return(mask)

  # A successful planting takes the place of natural seeding for this step.
  # The base engine has already consumed its usual cell-major random draws.
  planted <- rowsum(bitwShiftL(1L, species[eligible] - 1L),
                    position[eligible], reorder = FALSE)
  mask[as.integer(rownames(planted))] <- planted[, 1L]
  mask
}
