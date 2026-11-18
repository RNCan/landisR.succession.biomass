#include <Rcpp.h>
#include <landisR/mt19937.hpp>
#include <landisR/progress.hpp>
#include <landisR/array_sort.hpp>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

struct Offset {
  int row, col;
};

struct NeighborGroup {
  double probability;
  Offset offsets[4];
};

bool is_mature(const LogicalMatrix& mature, const IntegerVector& cell_state,
               const IntegerVector& positions, int row, int col,
               int nrow, int ncol, int species) {
  if (row < 0 || col < 0 || row >= nrow || col >= ncol) return false;
  const int cell = positions[row * ncol + col] - 1;
  return cell >= 0 && mature(cell_state[cell] - 1, species);
}

} // namespace

// [[Rcpp::export(rng = false)]]
IntegerVector bs52_v7_distance_order_cpp(const NumericVector& distance) {
  IntegerVector order(distance.size());
  std::iota(order.begin(), order.end(), 1);
  landis_array_sort::landis_sort(order.begin(), order.size(), [&](int a, int b) {
    return distance[a - 1] < distance[b - 1];
  });
  return order;
}

// [[Rcpp::export(rng = false)]]
List bs52_ward531_draw_cpp(LogicalMatrix mature, IntegerVector cell_state,
                          IntegerVector state_eco, IntegerVector shade_class,
                          IntegerMatrix targets, int nrow, int ncol,
                          IntegerVector positions, IntegerVector site_order,
                          NumericMatrix sep, NumericMatrix light, IntegerVector tolerance,
                          List groups, LogicalVector strict, IntegerVector mask,
                          IntegerVector excluded, std::string rng_state, double rng_offset,
                          bool deferred, IntegerVector pending_cell, IntegerVector pending_species,
                          LogicalVector pending_serotiny, NumericVector vegetative_probability) {
  const int n_sites = cell_state.size(), n_species = mature.ncol();
  if (n_species < 1 || n_species > 30 || mask.size() != n_sites ||
      targets.nrow() != n_sites || targets.ncol() != 2 ||
      mature.nrow() != state_eco.size() || shade_class.size() != state_eco.size() ||
      sep.ncol() != n_species || tolerance.size() != n_species ||
      groups.size() != n_species || strict.size() != n_species || light.ncol() != 6 ||
      nrow < 1 || ncol < 1 || static_cast<std::int64_t>(nrow) * ncol != positions.size() ||
      site_order.size() != n_sites)
    stop("Invalid Ward 5.3.1 input dimensions");

  landis_progress::Reporter progress;
  progress.phase(0, "Preparing Ward neighborhoods");
  landis_rng::Stream rng(0);
  rng.restore(rng_state);
  rng.draws = landis_rng::Stream::offset(rng_offset);
  for (int cell = 0; cell < n_sites; ++cell) {
    const int row = targets(cell, 0), col = targets(cell, 1), state = cell_state[cell];
    if (row < 1 || row > nrow || col < 1 || col > ncol || state < 1 || state > mature.nrow())
      stop("Invalid active position for Ward 5.3.1");
    const int eco = state_eco[state - 1], shade = shade_class[state - 1];
    if (eco < 1 || eco > sep.nrow() || shade < 0 || shade >= light.ncol())
      stop("Invalid state ecology for Ward 5.3.1");
  }
  for (int cell : excluded) {
    if (cell < 1 || cell > n_sites) stop("Invalid excluded position for Ward 5.3.1");
    mask[cell - 1] = -1;
  }
  std::vector<std::vector<NeighborGroup>> neighborhoods(n_species);
  std::vector<bool> universal(n_species, false);
  for (int sp = 0; sp < n_species; ++sp) {
    if (tolerance[sp] < 1 || tolerance[sp] > light.nrow())
      stop("Invalid shade tolerance for Ward 5.3.1");
    SEXP entry = groups[sp];
    if (Rf_isNull(entry)) {
      universal[sp] = true;
      continue;
    }
    List group(entry);
    IntegerVector dr = group["dr"], dc = group["dc"];
    NumericVector probability = group["p"];
    if (dr.size() != dc.size() || dr.size() != probability.size())
      stop("Invalid neighborhood for Ward 5.3.1");
    for (R_xlen_t k = 0; k < dr.size(); ++k) {
      const int r = dr[k], c = dc[k];
      const bool axis = c == 0;
      neighborhoods[sp].push_back({probability[k],
        {{r, c}, {axis ? 0 : -r, axis ? r : c}, {-r, -c}, {axis ? 0 : r, axis ? -r : -c}}});
    }
  }

  std::vector<unsigned int> serotiny(deferred ? n_sites : 0, 0);
  std::vector<unsigned int> resprout(deferred ? n_sites : 0, 0);
  if (deferred) {
    if (pending_cell.size() != pending_species.size() || pending_cell.size() != pending_serotiny.size() ||
        vegetative_probability.size() != n_species)
      stop("Invalid pending regeneration dimensions");
    for (R_xlen_t i = 0; i < pending_cell.size(); ++i) {
      const int cell = pending_cell[i], sp = pending_species[i];
      if (cell < 1 || cell > n_sites || sp < 1 || sp > n_species ||
          pending_serotiny[i] == NA_LOGICAL)
        stop("Invalid pending regeneration signal");
      auto& flags = pending_serotiny[i] ? serotiny : resprout;
      flags[cell - 1] |= 1u << (sp - 1);
    }
  }
  progress.start(0.14, 1, "Ward dispersal & recruitment", n_sites);
  for (int order = 0; order < n_sites; ++order) {
    progress.advance(order);
    if ((order & 16383) == 0) checkUserInterrupt();
    const int cell = site_order[order] - 1;
    if (mask[cell] == -1) {
      mask[cell] = 0;
      continue;
    }
    if (mask[cell] != 0) continue;  // Successful planting takes priority.
    const int state = cell_state[cell] - 1;
    const int eco = state_eco[state] - 1, shade = shade_class[state];
    const int row = targets(cell, 0) - 1, col = targets(cell, 1) - 1;
    if (deferred && (serotiny[cell] || resprout[cell])) {
      for (int sp = 0; sp < n_species; ++sp) {
        if (!(serotiny[cell] & (1u << sp))) continue;
        if (rng.next_double() >= light(tolerance[sp] - 1, shade)) continue;
        if (rng.next_double() < sep(eco, sp)) mask[cell] |= 1 << sp;
      }
      if (mask[cell] == 0) {
        for (int sp = 0; sp < n_species; ++sp) {
          if (!(resprout[cell] & (1u << sp))) continue;
          if (rng.next_double() >= light(tolerance[sp] - 1, shade)) continue;
          const double probability = static_cast<float>(vegetative_probability[sp]);
          if (rng.next_double() < probability) mask[cell] |= 1 << sp;
        }
      }
      if (mask[cell] != 0) continue;
    }
    for (int sp = 0; sp < n_species; ++sp) {
      if (rng.next_double() >= light(tolerance[sp] - 1, shade)) continue;
      if (rng.next_double() >= sep(eco, sp) || strict[sp]) continue;
      bool recruited = universal[sp] || mature(state, sp);
      if (!recruited) {
        for (const NeighborGroup& group : neighborhoods[sp]) {
          for (const Offset& offset : group.offsets) {
            // A draw is consumed even when the neighbor is outside the grid.
            if (rng.next_double() >= group.probability) continue;
            if (is_mature(mature, cell_state, positions, row + offset.row,
                          col + offset.col, nrow, ncol, sp)) {
              recruited = true;
              break;
            }
          }
          if (recruited) break;
        }
      }
      if (recruited) mask[cell] |= 1 << sp;
    }
  }
  return List::create(_["mask"] = mask,
                      _["rng_state"] = rng.encode(),
                      _["rng_offset"] = static_cast<double>(rng.draws));
}
