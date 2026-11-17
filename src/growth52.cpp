#include <Rcpp.h>
#include <landisR/csr_store.hpp>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Rcpp;
using landis_csr::make_csr_store;
using landis_csr::numeric_column;
using landis_csr::validate_csr_store;



// [[Rcpp::export]]
std::string bs52_growth_kernel_version_cpp() {
  return "v3.11-ward-interior-fastpath";
}


// [[Rcpp::plugins(cpp11)]]

namespace {

inline void bs52_progress_update(const std::string& label,
                                 int current,
                                 int total,
                                 int& last_percent,
                                 bool show_progress) {
  if (!show_progress || total <= 0) return;

  int percent = static_cast<int>(
    std::floor(
      100.0 * static_cast<double>(current) /
      static_cast<double>(total)
    )
  );

  if (current >= total) percent = 100;

  // Redraw only every 2 percentage points.
  if (percent < 100 &&
      last_percent >= 0 &&
      percent < last_percent + 2) {
    return;
  }

  if (percent == last_percent && current < total) return;

  const int width = 30;
  const int filled = static_cast<int>(
    std::floor(
      static_cast<double>(width) *
      static_cast<double>(percent) / 100.0
    )
  );

  Rcpp::Rcout << "\r  " << label << " [";

  for (int i = 0; i < width; ++i) {
    Rcpp::Rcout << (i < filled ? "=" : " ");
  }

  Rcpp::Rcout << "] ";

  if (percent < 10) {
    Rcpp::Rcout << "  ";
  } else if (percent < 100) {
    Rcpp::Rcout << " ";
  }

  Rcpp::Rcout << percent << "%";
  Rcpp::Rcout.flush();

  last_percent = percent;

  if (current >= total) {
    Rcpp::Rcout << "\n";
    Rcpp::Rcout.flush();
  }
}


struct CohortState {
  std::vector<std::string> species;
  std::vector<double> age;
  std::vector<double> biomass;
};

// The living-biomass update is also needed by the dead-pool bookkeeping.
// Keep all intermediate values in one structure so the optional dead-pool
// path follows the same annual calculation as the validated living path.
struct CohortChange {
  double new_biomass;
  double raw_anpp;
  double actual_anpp;
  double mortality_age;
  double total_mortality;
};

struct DeadStaticParameters {
  double leaf_longevity;
  double woody_decay;
  double leaf_lignin;
};

struct DeadInputEvent {
  int year;
  std::string species;
  double wood;
  double litter;
};

// A group of sites can share an identical living-cohort state while retaining
// different dead-pool stocks from earlier stochastic recruitment.  Keep the
// pool vectors in C++ during an annual block so that every cohort event is
// applied in its original order without creating an R data.frame per event.
struct DeadPoolVectors {
  std::vector<double> wood;
  std::vector<double> wood_decay;
  std::vector<double> litter;
  std::vector<double> litter_decay;

  DeadPoolVectors(const NumericVector& wood_in,
                  const NumericVector& wood_decay_in,
                  const NumericVector& litter_in,
                  const NumericVector& litter_decay_in) :
    wood(as<std::vector<double> >(wood_in)),
    wood_decay(as<std::vector<double> >(wood_decay_in)),
    litter(as<std::vector<double> >(litter_in)),
    litter_decay(as<std::vector<double> >(litter_decay_in)) {
    const std::size_t n = wood.size();
    if (wood_decay.size() != n || litter.size() != n ||
        litter_decay.size() != n)
      stop("Incompatible dead-pool vectors");
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(wood[i]) || !std::isfinite(wood_decay[i]) ||
          !std::isfinite(litter[i]) || !std::isfinite(litter_decay[i]) ||
          wood[i] < 0.0 || litter[i] < 0.0)
        stop("Invalid dead-pool vector state");
    }
  }

  void add(double wood_input, double litter_input,
           const DeadStaticParameters& parameters, double aet) {
    if (!std::isfinite(wood_input) || !std::isfinite(litter_input) ||
        wood_input < 0.0 || litter_input < 0.0)
      stop("Invalid dead-pool input");
    if (wood_input > 65535.0)
      stop("Woody input outside supported C# UInt16 range");
    if (wood_input == 0.0 && litter_input == 0.0) return;
    const double litter_k = (-0.5365 + 0.00241 * aet) -
      ((-0.01586 + 0.000056 * aet) * parameters.leaf_lignin * 100.0);
    for (std::size_t i = 0; i < wood.size(); ++i) {
      if (wood_input > 0.0) {
        const double total = wood[i] + wood_input;
        wood_decay[i] = (wood[i] * wood_decay[i] +
          wood_input * parameters.woody_decay) / total;
        wood[i] = total;
      }
      if (litter_input > 0.0) {
        const double total = litter[i] + litter_input;
        litter_decay[i] = (litter[i] * litter_decay[i] +
          litter_input * litter_k) / total;
        litter[i] = total;
      }
    }
  }

  void decompose() {
    for (std::size_t i = 0; i < wood.size(); ++i) {
      const double next_wood = std::trunc(wood[i] * std::exp(-wood_decay[i]));
      const double next_litter = std::trunc(
        litter[i] * std::exp(-litter_decay[i]));
      if (!std::isfinite(next_wood) || !std::isfinite(next_litter) ||
          next_wood > 4294967295.0 || next_litter > 4294967295.0)
        stop("Dead pool outside supported C# UInt32 range");
      wood[i] = next_wood;
      litter[i] = next_litter;
    }
  }
};

std::vector<std::string> unique_species(const std::vector<std::string>& x) {
  std::vector<std::string> out;
  std::unordered_map<std::string, bool> seen;
  for (const auto& value : x) {
    if (!seen[value]) {
      seen[value] = true;
      out.push_back(value);
    }
  }
  return out;
}


// Equivalent to CohortBiomass.ComputeChange for the living-biomass fields
// used in this port.  The source truncates only the returned biomass delta.
CohortChange cohort_change(double biomass, double age, double longevity,
                            double mortality_shape, double growth_shape,
                            double max_anpp, double max_biomass,
                            double site_biomass, double competition,
                            bool spinup,
                            double spinup_mortality_fraction) {
  double mortality_age = std::min(
    biomass * std::exp(age / longevity * mortality_shape) /
      std::exp(mortality_shape), biomass);
  if (spinup && spinup_mortality_fraction > 0.0) {
    mortality_age = std::min(
      mortality_age + biomass * spinup_mortality_fraction, biomass);
  }
  const double potential_biomass = std::max(
    0.0, std::max(1.0, max_biomass - site_biomass + biomass));
  const double b_ap = biomass / potential_biomass;
  const double power = R_pow(b_ap, growth_shape);
  double raw_anpp = max_anpp * std::exp(1.0) * power *
    std::exp(-power) * competition;
  raw_anpp = std::min(max_anpp * competition, raw_anpp);
  const double actual_anpp = std::max(1.0, raw_anpp - mortality_age);
  double mortality_growth = b_ap > 1.0 ? max_anpp * competition :
    max_anpp * (2.0 * b_ap) / (1.0 + b_ap) * competition;
  mortality_growth = std::min(
    biomass, std::min(mortality_growth, max_anpp * competition));
  mortality_growth = std::min(
    std::max(0.0, mortality_growth - mortality_age), actual_anpp);
  CohortChange out;
  out.new_biomass = biomass + std::trunc(
    actual_anpp - mortality_age - mortality_growth);
  out.raw_anpp = raw_anpp;
  out.actual_anpp = actual_anpp;
  out.mortality_age = mortality_age;
  out.total_mortality = mortality_age + mortality_growth;
  return out;
}

double new_biomass(double biomass, double age, double longevity,
                   double mortality_shape, double growth_shape,
                   double max_anpp, double max_biomass,
                   double site_biomass, double competition, bool spinup,
                   double spinup_mortality_fraction) {
  return cohort_change(biomass, age, longevity, mortality_shape,
                        growth_shape, max_anpp, max_biomass, site_biomass,
                        competition, spinup,
                        spinup_mortality_fraction).new_biomass;
}

std::unordered_map<std::string, DeadStaticParameters>
dead_static_parameters(const DataFrame& parameters) {
  CharacterVector species = parameters["species"];
  NumericVector leaf_longevity = numeric_column(parameters, "leaf_longevity");
  NumericVector woody_decay = numeric_column(parameters, "woody_decay");
  NumericVector leaf_lignin = numeric_column(parameters, "leaf_lignin");
  if (species.size() != leaf_longevity.size() ||
      species.size() != woody_decay.size() ||
      species.size() != leaf_lignin.size())
    stop("Incompatible static dead-pool parameters");
  std::unordered_map<std::string, DeadStaticParameters> out;
  for (int i = 0; i < species.size(); ++i) {
    const std::string name = as<std::string>(species[i]);
    if (!std::isfinite(leaf_longevity[i]) || leaf_longevity[i] < 1.0 ||
        !std::isfinite(woody_decay[i]) || !std::isfinite(leaf_lignin[i]))
      stop("Invalid static dead-pool parameters");
    if (out.find(name) != out.end())
      stop("Duplicate static dead-pool species");
    out[name] = DeadStaticParameters{
      leaf_longevity[i], woody_decay[i], leaf_lignin[i]
    };
  }
  return out;
}

void combine_young(CohortState& state,
                   const std::vector<std::string>& order,
                   int timestep) {
  if (timestep <= 1) return;
  std::vector<bool> keep(state.age.size(), true);
  std::vector<std::string> add_species;
  std::vector<double> add_age;
  std::vector<double> add_biomass;
  for (const auto& s : order) {
    bool found = false;
    double total = 0.0;
    for (std::size_t i = 0; i < state.age.size(); ++i) {
      if (state.species[i] == s && state.age[i] <= timestep) {
        keep[i] = false;
        found = true;
        total += state.biomass[i];
      }
    }
    if (found) {
      add_species.push_back(s);
      add_age.push_back(timestep - 1.0);
      add_biomass.push_back(total);
    }
  }
  CohortState next;
  next.species.reserve(state.species.size() + add_species.size());
  next.age.reserve(state.age.size() + add_age.size());
  next.biomass.reserve(state.biomass.size() + add_biomass.size());
  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (keep[i]) {
      next.species.push_back(state.species[i]);
      next.age.push_back(state.age[i]);
      next.biomass.push_back(state.biomass[i]);
    }
  }
  next.species.insert(next.species.end(), add_species.begin(), add_species.end());
  next.age.insert(next.age.end(), add_age.begin(), add_age.end());
  next.biomass.insert(next.biomass.end(), add_biomass.begin(), add_biomass.end());
  state = std::move(next);
}

void grow_year(CohortState& state, const DataFrame& params, int timestep,
               bool do_combine, double spinup_fraction, bool spinup,
               int year = 0,
               const std::unordered_map<std::string, DeadStaticParameters>*
                 dead_static = NULL,
               std::vector<DeadInputEvent>* dead_events = NULL,
               DeadPoolVectors* dead_pools = NULL,
               double dead_aet = 0.0) {
  if (state.age.empty()) return;
  const std::vector<std::string> spp_order = unique_species(state.species);
  double site_biomass = 0.0;
  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (state.age[i] >= timestep) site_biomass += state.biomass[i];
  }
  if (do_combine) combine_young(state, spp_order, timestep);
  if (state.age.empty()) return;

  CohortState grouped;
  grouped.species.reserve(state.species.size());
  grouped.age.reserve(state.age.size());
  grouped.biomass.reserve(state.biomass.size());
  std::vector<std::vector<int> > queues(spp_order.size());
  for (std::size_t s = 0; s < spp_order.size(); ++s) {
    for (std::size_t i = 0; i < state.species.size(); ++i) {
      if (state.species[i] == spp_order[s]) {
        const int id = static_cast<int>(grouped.age.size());
        grouped.species.push_back(state.species[i]);
        grouped.age.push_back(state.age[i]);
        grouped.biomass.push_back(state.biomass[i]);
        queues[s].push_back(id);
      }
    }
  }
  state = std::move(grouped);

  CharacterVector p_species = params["species"];
  NumericVector longevity = numeric_column(params, "longevity");
  NumericVector mortality_shape = numeric_column(params, "mortality_shape");
  NumericVector growth_shape = numeric_column(params, "growth_shape");
  NumericVector max_anpp = numeric_column(params, "max_anpp");
  NumericVector max_biomass = numeric_column(params, "max_biomass");
  std::unordered_map<std::string, int> param_index;
  for (int i = 0; i < p_species.size(); ++i)
    param_index[as<std::string>(p_species[i])] = i;
  std::vector<int> pi(spp_order.size());
  for (std::size_t s = 0; s < spp_order.size(); ++s) {
    const auto found = param_index.find(spp_order[s]);
    if (found == param_index.end()) stop("Missing species parameters");
    pi[s] = found->second;
  }

  std::vector<std::size_t> queue_position(spp_order.size(), 0);
  std::vector<bool> alive(state.age.size(), true);
  auto head_age = [&](int s) {
    return state.age[queues[s][queue_position[s]]];
  };
  std::vector<int> iterators;
  for (std::size_t s = 0; s < spp_order.size(); ++s) {
    auto at = iterators.end();
    for (auto it = iterators.begin(); it != iterators.end(); ++it) {
      if (head_age(*it) < head_age(static_cast<int>(s))) {
        at = it;
        break;
      }
    }
    iterators.insert(at, static_cast<int>(s));
  }

  while (!iterators.empty()) {
    const int s = iterators.front();
    const int id = queues[s][queue_position[s]];
    if (state.age[id] >= longevity[pi[s]]) {
      if (dead_events != NULL || dead_pools != NULL) {
        const auto static_found = dead_static->find(spp_order[s]);
        if (static_found == dead_static->end())
          stop("Missing static dead-pool species parameters");
        const double focal_weight = std::max(
          R_pow(state.biomass[id], 0.95), 1.0);
        long double other_weight_accumulator = 0.0L;
        for (std::size_t j = 0; j < state.age.size(); ++j) {
          if (!alive[j]) continue;
          // Same age argument as the R reference's bs52_dead_death path.
          const bool include = state.age[j] + 1.0 != state.age[id] ||
            state.species[j] != spp_order[s];
          if (include)
            other_weight_accumulator +=
              std::max(R_pow(state.biomass[j], 0.95), 1.0);
        }
        const double competition = focal_weight /
          (focal_weight + static_cast<double>(other_weight_accumulator));
        const CohortChange death_change = cohort_change(
          state.biomass[id], state.age[id], longevity[pi[s]],
          mortality_shape[pi[s]], growth_shape[pi[s]], max_anpp[pi[s]],
          max_biomass[pi[s]], site_biomass, competition, spinup,
          spinup_fraction);
        const double death_anpp = std::max(
          0.0, death_change.raw_anpp - death_change.mortality_age);
        const double foliage = std::trunc(std::min(
          std::max(death_anpp * 0.35 * static_found->second.leaf_longevity,
                   state.biomass[id] * 0.025),
          state.biomass[id] * 0.35));
        const double wood_input = state.biomass[id] - foliage;
        if (dead_events != NULL) {
          dead_events->push_back(DeadInputEvent{
            year, spp_order[s], wood_input, foliage
          });
        }
        if (dead_pools != NULL)
          dead_pools->add(wood_input, foliage, static_found->second, dead_aet);
      }
      alive[id] = false;
    } else {
      const double updated_age = state.age[id] + 1.0;
      const double focal_weight = std::max(
        R_pow(state.biomass[id], 0.95), 1.0);
      long double other_weight_accumulator = 0.0L;
      for (std::size_t j = 0; j < state.age.size(); ++j) {
        if (!alive[j]) continue;
        const bool include = state.age[j] + 1.0 != updated_age ||
          state.species[j] != spp_order[s];
        if (include)
          other_weight_accumulator +=
            std::max(R_pow(state.biomass[j], 0.95), 1.0);
      }
      const double other_weight =
        static_cast<double>(other_weight_accumulator);
      const double competition = focal_weight /
        (focal_weight + other_weight);
      const CohortChange change = cohort_change(
        state.biomass[id], updated_age, longevity[pi[s]],
        mortality_shape[pi[s]], growth_shape[pi[s]], max_anpp[pi[s]],
        max_biomass[pi[s]], site_biomass, competition, spinup,
        spinup_fraction);
      if (dead_events != NULL || dead_pools != NULL) {
        const auto static_found = dead_static->find(spp_order[s]);
        if (static_found == dead_static->end())
          stop("Missing static dead-pool species parameters");
        if (change.new_biomass <= 0.0)
          stop("Zero post-growth biomass in dead-pool bookkeeping");
        const double annual_leaf = change.actual_anpp * 0.35;
        const double standing_nonwood = annual_leaf *
          (static_found->second.leaf_longevity - 1.0);
        const double remaining_mortality = change.total_mortality - annual_leaf;
        const double mortality_nonwood = std::max(
          0.0, remaining_mortality * standing_nonwood / change.new_biomass);
        const double mortality_wood = std::max(
          0.0, remaining_mortality - mortality_nonwood);
        if (mortality_wood > 65535.0)
          stop("Woody input outside supported C# UInt16 range");
        const double wood_input = std::trunc(mortality_wood);
        const double litter_input = annual_leaf + mortality_nonwood;
        if (dead_events != NULL) {
          dead_events->push_back(DeadInputEvent{
            year, spp_order[s], wood_input, litter_input
          });
        }
        if (dead_pools != NULL)
          dead_pools->add(wood_input, litter_input, static_found->second,
                          dead_aet);
      }
      state.age[id] = updated_age;
      state.biomass[id] = std::max(0.0, change.new_biomass);
      if (state.biomass[id] == 0.0) alive[id] = false;
    }

    ++queue_position[s];
    if (queue_position[s] >= queues[s].size()) {
      iterators.erase(iterators.begin());
    } else if (iterators.size() > 1 &&
               head_age(s) < head_age(iterators[1])) {
      iterators.erase(iterators.begin());
      auto at = iterators.end();
      for (auto it = iterators.begin(); it != iterators.end(); ++it) {
        if (head_age(*it) < head_age(s)) {
          at = it;
          break;
        }
      }
      iterators.insert(at, s);
    }
  }

  CohortState kept;
  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (alive[i]) {
      kept.species.push_back(state.species[i]);
      kept.age.push_back(state.age[i]);
      kept.biomass.push_back(state.biomass[i]);
    }
  }
  state = std::move(kept);
}


// ============================================================================
// Indexed living-growth path for large52 v2.1.
//
// The validated reference path keeps species as strings because that mirrors
// LANDIS input tables closely. For very large landscapes this becomes costly:
// grow_year() rebuilds species-name maps for every state and every simulated
// year. The indexed path below keeps the exact cohort order and floating-point
// operations, but translates species names to stable integer IDs once per
// state/block and compiles the annual parameter tables once per eco/block.
// ============================================================================

struct CohortStateIndexed {
  std::vector<int> species;
  std::vector<double> age;
  std::vector<double> biomass;
};

struct GrowthParametersIndexed {
  std::vector<double> longevity;
  std::vector<double> mortality_shape;
  std::vector<double> growth_shape;
  std::vector<double> max_anpp;
  std::vector<double> max_biomass;
};

std::vector<int> unique_species_indexed(const std::vector<int>& x,
                                        int n_species) {
  std::vector<int> out;
  std::vector<unsigned char> seen(
    static_cast<std::size_t>(n_species),
    static_cast<unsigned char>(0)
  );
  out.reserve(static_cast<std::size_t>(n_species));

  for (const int value : x) {
    if (value < 0 || value >= n_species)
      stop("Invalid indexed cohort species");
    if (!seen[static_cast<std::size_t>(value)]) {
      seen[static_cast<std::size_t>(value)] = 1;
      out.push_back(value);
    }
  }

  return out;
}

void combine_young_indexed(CohortStateIndexed& state,
                           const std::vector<int>& order,
                           int timestep) {
  if (timestep <= 1) return;

  std::vector<bool> keep(state.age.size(), true);
  std::vector<int> add_species;
  std::vector<double> add_age;
  std::vector<double> add_biomass;

  add_species.reserve(order.size());
  add_age.reserve(order.size());
  add_biomass.reserve(order.size());

  for (const int s : order) {
    bool found = false;
    double total = 0.0;

    for (std::size_t i = 0; i < state.age.size(); ++i) {
      if (state.species[i] == s && state.age[i] <= timestep) {
        keep[i] = false;
        found = true;
        total += state.biomass[i];
      }
    }

    if (found) {
      add_species.push_back(s);
      add_age.push_back(timestep - 1.0);
      add_biomass.push_back(total);
    }
  }

  CohortStateIndexed next;
  next.species.reserve(state.species.size() + add_species.size());
  next.age.reserve(state.age.size() + add_age.size());
  next.biomass.reserve(state.biomass.size() + add_biomass.size());

  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (keep[i]) {
      next.species.push_back(state.species[i]);
      next.age.push_back(state.age[i]);
      next.biomass.push_back(state.biomass[i]);
    }
  }

  next.species.insert(
    next.species.end(),
    add_species.begin(),
    add_species.end()
  );
  next.age.insert(
    next.age.end(),
    add_age.begin(),
    add_age.end()
  );
  next.biomass.insert(
    next.biomass.end(),
    add_biomass.begin(),
    add_biomass.end()
  );

  state = std::move(next);
}

void grow_year_indexed(CohortStateIndexed& state,
                       const GrowthParametersIndexed& params,
                       int timestep,
                       bool do_combine,
                       double spinup_fraction,
                       bool spinup) {
  if (state.age.empty()) return;

  const int n_species = static_cast<int>(params.longevity.size());

  if (params.mortality_shape.size() != params.longevity.size() ||
      params.growth_shape.size() != params.longevity.size() ||
      params.max_anpp.size() != params.longevity.size() ||
      params.max_biomass.size() != params.longevity.size())
    stop("Incompatible indexed growth parameters");

  const std::vector<int> spp_order = unique_species_indexed(
    state.species,
    n_species
  );

  double site_biomass = 0.0;
  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (state.age[i] >= timestep)
      site_biomass += state.biomass[i];
  }

  if (do_combine)
    combine_young_indexed(state, spp_order, timestep);

  if (state.age.empty()) return;

  CohortStateIndexed grouped;
  grouped.species.reserve(state.species.size());
  grouped.age.reserve(state.age.size());
  grouped.biomass.reserve(state.biomass.size());

  std::vector<std::vector<int> > queues(spp_order.size());

  for (std::size_t s = 0; s < spp_order.size(); ++s) {
    for (std::size_t i = 0; i < state.species.size(); ++i) {
      if (state.species[i] == spp_order[s]) {
        const int id = static_cast<int>(grouped.age.size());
        grouped.species.push_back(state.species[i]);
        grouped.age.push_back(state.age[i]);
        grouped.biomass.push_back(state.biomass[i]);
        queues[s].push_back(id);
      }
    }
  }

  state = std::move(grouped);

  std::vector<std::size_t> queue_position(
    spp_order.size(),
    static_cast<std::size_t>(0)
  );
  std::vector<bool> alive(state.age.size(), true);

  // v3.2: cache the competition weight for each live cohort.
  //
  // The previous validated path recomputed
  //   max(R_pow(biomass, 0.95), 1.0)
  // for every cohort inside every focal-cohort competition scan.
  // The biomass of a cohort changes only when that cohort itself is processed,
  // so we can cache this pure function and refresh exactly once after an
  // update. The subsequent long-double summation still visits cohorts in the
  // identical index order, preserving the floating-point accumulation path.
  std::vector<double> competition_weight(
    state.biomass.size(),
    0.0
  );

  for (std::size_t i = 0; i < state.biomass.size(); ++i) {
    competition_weight[i] = std::max(
      R_pow(
        state.biomass[i],
        0.95
      ),
      1.0
    );
  }

  auto head_age = [&](int s) {
    return state.age[
      queues[static_cast<std::size_t>(s)][
        queue_position[static_cast<std::size_t>(s)]
      ]
    ];
  };

  std::vector<int> iterators;
  iterators.reserve(spp_order.size());

  for (std::size_t s = 0; s < spp_order.size(); ++s) {
    auto at = iterators.end();

    for (auto it = iterators.begin(); it != iterators.end(); ++it) {
      if (head_age(*it) < head_age(static_cast<int>(s))) {
        at = it;
        break;
      }
    }

    iterators.insert(
      at,
      static_cast<int>(s)
    );
  }

  while (!iterators.empty()) {
    const int s = iterators.front();
    const std::size_t su = static_cast<std::size_t>(s);
    const int id = queues[su][queue_position[su]];
    const int species_id = spp_order[su];

    if (species_id < 0 || species_id >= n_species)
      stop("Invalid indexed species in growth");

    if (state.age[id] >= params.longevity[species_id]) {
      alive[id] = false;
    } else {
      const double updated_age = state.age[id] + 1.0;

      const double focal_weight =
        competition_weight[
          static_cast<std::size_t>(id)
        ];

      long double other_weight_accumulator = 0.0L;

      for (std::size_t j = 0; j < state.age.size(); ++j) {
        if (!alive[j]) continue;

        const bool include =
          state.age[j] + 1.0 != updated_age ||
          state.species[j] != species_id;

        if (include) {
          other_weight_accumulator +=
            competition_weight[j];
        }
      }

      const double other_weight =
        static_cast<double>(other_weight_accumulator);

      const double competition = focal_weight /
        (focal_weight + other_weight);

      const CohortChange change = cohort_change(
        state.biomass[id],
        updated_age,
        params.longevity[species_id],
        params.mortality_shape[species_id],
        params.growth_shape[species_id],
        params.max_anpp[species_id],
        params.max_biomass[species_id],
        site_biomass,
        competition,
        spinup,
        spinup_fraction
      );

      state.age[id] = updated_age;
      state.biomass[id] = std::max(
        0.0,
        change.new_biomass
      );

      if (state.biomass[id] == 0.0) {
        alive[id] = false;
        competition_weight[
          static_cast<std::size_t>(id)
        ] = 0.0;
      } else {
        competition_weight[
          static_cast<std::size_t>(id)
        ] = std::max(
          R_pow(
            state.biomass[id],
            0.95
          ),
          1.0
        );
      }
    }

    ++queue_position[su];

    if (queue_position[su] >= queues[su].size()) {
      iterators.erase(iterators.begin());
    } else if (
      iterators.size() > 1 &&
      head_age(s) < head_age(iterators[1])
    ) {
      iterators.erase(iterators.begin());

      auto at = iterators.end();

      for (auto it = iterators.begin(); it != iterators.end(); ++it) {
        if (head_age(*it) < head_age(s)) {
          at = it;
          break;
        }
      }

      iterators.insert(at, s);
    }
  }

  CohortStateIndexed kept;
  kept.species.reserve(state.species.size());
  kept.age.reserve(state.age.size());
  kept.biomass.reserve(state.biomass.size());

  for (std::size_t i = 0; i < state.age.size(); ++i) {
    if (alive[i]) {
      kept.species.push_back(state.species[i]);
      kept.age.push_back(state.age[i]);
      kept.biomass.push_back(state.biomass[i]);
    }
  }

  state = std::move(kept);
}

GrowthParametersIndexed compile_growth_parameters_indexed(
    const DataFrame& params,
    const std::unordered_map<std::string, int>& species_index,
    int n_species) {

  CharacterVector p_species = params["species"];
  NumericVector longevity = numeric_column(params, "longevity");
  NumericVector mortality_shape = numeric_column(params, "mortality_shape");
  NumericVector growth_shape = numeric_column(params, "growth_shape");
  NumericVector max_anpp = numeric_column(params, "max_anpp");
  NumericVector max_biomass = numeric_column(params, "max_biomass");

  const int n = p_species.size();

  if (longevity.size() != n ||
      mortality_shape.size() != n ||
      growth_shape.size() != n ||
      max_anpp.size() != n ||
      max_biomass.size() != n)
    stop("Incompatible indexed growth parameter columns");

  GrowthParametersIndexed out;
  out.longevity.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );
  out.mortality_shape.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );
  out.growth_shape.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );
  out.max_anpp.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );
  out.max_biomass.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );

  std::vector<unsigned char> seen(
    static_cast<std::size_t>(n_species),
    static_cast<unsigned char>(0)
  );

  for (int i = 0; i < n; ++i) {
    const std::string name = as<std::string>(p_species[i]);
    const auto found = species_index.find(name);

    if (found == species_index.end())
      stop("Unknown species in indexed growth parameters");

    const int id = found->second;

    if (seen[static_cast<std::size_t>(id)])
      stop("Duplicate species in indexed growth parameters");

    seen[static_cast<std::size_t>(id)] = 1;

    out.longevity[id] = longevity[i];
    out.mortality_shape[id] = mortality_shape[i];
    out.growth_shape[id] = growth_shape[i];
    out.max_anpp[id] = max_anpp[i];
    out.max_biomass[id] = max_biomass[i];
  }

  for (int j = 0; j < n_species; ++j) {
    if (!seen[static_cast<std::size_t>(j)])
      stop("Missing species in indexed growth parameters");
  }

  return out;
}

DataFrame as_data_frame_indexed(
    const CohortStateIndexed& state,
    const CharacterVector& species_names) {

  const int n = static_cast<int>(state.species.size());

  CharacterVector species(n);
  NumericVector age(n);
  NumericVector biomass(n);

  for (int i = 0; i < n; ++i) {
    const int id = state.species[static_cast<std::size_t>(i)];

    if (id < 0 || id >= species_names.size())
      stop("Invalid indexed species in output state");

    species[i] = species_names[id];
    age[i] = state.age[static_cast<std::size_t>(i)];
    biomass[i] = state.biomass[static_cast<std::size_t>(i)];
  }

  return DataFrame::create(
    _["species"] = species,
    _["age"] = age,
    _["biomass"] = biomass
  );
}

DataFrame dead_events_as_data_frame(const std::vector<DeadInputEvent>& events) {
  const int n = static_cast<int>(events.size());
  IntegerVector year(n);
  CharacterVector species(n);
  NumericVector wood(n), litter(n);
  for (int i = 0; i < n; ++i) {
    year[i] = events[i].year;
    species[i] = events[i].species;
    wood[i] = events[i].wood;
    litter[i] = events[i].litter;
  }
  return DataFrame::create(
    _["year"] = year, _["species"] = species,
    _["wood"] = wood, _["litter"] = litter
  );
}

DataFrame as_data_frame(const CohortState& state) {
  return DataFrame::create(
    _["species"] = wrap(state.species),
    _["age"] = wrap(state.age),
    _["biomass"] = wrap(state.biomass)
  );
}

DataFrame record_landscape_impl(List cohorts, IntegerVector ids, CharacterVector species_names, double year) {
  int ns=species_names.size(), nsite=ids.size(), n=nsite*(ns+1), z=0;
  std::unordered_map<std::string,int> ix; for(int j=0;j<ns;++j) ix[as<std::string>(species_names[j])]=j;
  NumericVector years(n,year), biomass(n); IntegerVector cell_id(n), max_age(n); CharacterVector species(n);
  for(int k=0;k<nsite;++k) { NumericVector total(ns); int oldest=0; DataFrame d=cohorts[ids[k]-1]; CharacterVector sp=d["species"]; NumericVector age=numeric_column(d,"age"), b=numeric_column(d,"biomass");
    for(int i=0;i<age.size();++i) { auto it=ix.find(as<std::string>(sp[i])); if(it!=ix.end()) total[it->second]+=b[i]; if(age[i]>oldest) oldest=age[i]; }
    for(int j=0;j<ns;++j,++z){cell_id[z]=ids[k];species[z]=species_names[j];biomass[z]=total[j];max_age[z]=oldest;}
    cell_id[z]=ids[k];species[z]="TotalBiomass";biomass[z]=std::accumulate(total.begin(),total.end(),0.0);max_age[z]=oldest;++z;
  }
  return DataFrame::create(_["year"]=years,_["cell_id"]=cell_id,_["species"]=species,_["biomass"]=biomass,_["max_age"]=max_age);
}

} // namespace

// [[Rcpp::export]]
DataFrame bs52_record_landscape_cpp(List cohorts, IntegerVector ids,
                                    CharacterVector species_names, double year) {
  return record_landscape_impl(cohorts, ids, species_names, year);
}

// [[Rcpp::export]]
DataFrame bs52_grow_live_block_cpp(DataFrame cohorts, List params_by_year,
                                   int timestep, bool combine_last = true,
                                   double spinup_fraction = 0.0,
                                   bool spinup = false) {
  CohortState state;
  state.species = as<std::vector<std::string> >(cohorts["species"]);
  state.age = as<std::vector<double> >(cohorts["age"]);
  state.biomass = as<std::vector<double> >(cohorts["biomass"]);
  const int n_years = params_by_year.size();
  for (int i = 0; i < n_years; ++i) {
    grow_year(state, as<DataFrame>(params_by_year[i]), timestep,
              combine_last && i == n_years - 1,
              spinup_fraction, spinup);
  }
  return as_data_frame(state);
}

// [[Rcpp::export]]
List bs52_grow_live_block_dead_cpp(DataFrame cohorts, List params_by_year,
                                   DataFrame static_parameters,
                                   IntegerVector years, int timestep,
                                   bool combine_last = true,
                                   double spinup_fraction = 0.0,
                                   bool spinup = false) {
  if (years.size() != params_by_year.size())
    stop("One absolute year is required for each parameter table");
  CohortState state;
  state.species = as<std::vector<std::string> >(cohorts["species"]);
  state.age = as<std::vector<double> >(cohorts["age"]);
  state.biomass = as<std::vector<double> >(cohorts["biomass"]);
  const std::unordered_map<std::string, DeadStaticParameters> dead_static =
    dead_static_parameters(static_parameters);
  std::vector<DeadInputEvent> events;
  const int n_years = params_by_year.size();
  for (int i = 0; i < n_years; ++i) {
    grow_year(state, as<DataFrame>(params_by_year[i]), timestep,
              combine_last && i == n_years - 1,
              spinup_fraction, spinup, years[i], &dead_static, &events);
  }
  return List::create(
    _["cohorts"] = as_data_frame(state),
    _["events"] = dead_events_as_data_frame(events)
  );
}

// [[Rcpp::export]]
List bs52_grow_live_block_dead_pools_cpp(
    DataFrame cohorts, List params_by_year, DataFrame static_parameters,
    double aet, IntegerVector years, int timestep,
    NumericVector wood, NumericVector wood_decay,
    NumericVector litter, NumericVector litter_decay,
    bool combine_last = true, double spinup_fraction = 0.0,
    bool spinup = false) {
  if (years.size() != params_by_year.size())
    stop("One absolute year is required for each parameter table");
  if (!std::isfinite(aet)) stop("Invalid AET for dead pools");
  CohortState state;
  state.species = as<std::vector<std::string> >(cohorts["species"]);
  state.age = as<std::vector<double> >(cohorts["age"]);
  state.biomass = as<std::vector<double> >(cohorts["biomass"]);
  const std::unordered_map<std::string, DeadStaticParameters> dead_static =
    dead_static_parameters(static_parameters);
  DeadPoolVectors pools(wood, wood_decay, litter, litter_decay);
  const int n_years = params_by_year.size();
  // Fire mortality happens after the succession extension's last annual
  // update.  At that moment SiteVars.TotalBiomass still holds the non-young
  // stock calculated immediately before that final annual growth calculation.
  // Return it so the prescribed-fire bridge can reproduce
  // Cohort.ComputeNonWoodyBiomass without re-running the block.
  double pre_final_site_biomass = NA_REAL;
  for (int i = 0; i < n_years; ++i) {
    if (i == n_years - 1) {
      pre_final_site_biomass = 0.0;
      for (std::size_t j = 0; j < state.age.size(); ++j) {
        if (state.age[j] >= timestep)
          pre_final_site_biomass += state.biomass[j];
      }
    }
    grow_year(state, as<DataFrame>(params_by_year[i]), timestep,
              combine_last && i == n_years - 1,
              spinup_fraction, spinup, years[i], &dead_static, NULL,
              &pools, aet);
    // Pool.Decompose runs after all cohort inputs for a calendar year, even
    // when the representative cohort state is empty.
    pools.decompose();
  }
  return List::create(
    _["cohorts"] = as_data_frame(state),
    _["wood"] = wrap(pools.wood),
    _["wood_decay"] = wrap(pools.wood_decay),
    _["litter"] = wrap(pools.litter),
    _["litter_decay"] = wrap(pools.litter_decay),
    _["pre_final_site_biomass"] = pre_final_site_biomass
  );
}

// Spatial accumulation for Ward dispersal. Kernel probabilities and their
// logarithms are calculated in R so the validated floating-point pathway is
// preserved. Targets and offsets are one-based, as in the R implementation.
// [[Rcpp::export]]
NumericVector bs52_ward_failure_log_cpp(LogicalMatrix mature,
                                        LogicalMatrix active,
                                        IntegerMatrix targets,
                                        IntegerVector neighbor_row,
                                        IntegerVector neighbor_col,
                                        NumericVector log_failure,
                                        bool ward2019 = true) {
  const int n_targets = targets.nrow();
  const int n_neighbors = neighbor_row.size();
  if (neighbor_col.size() != n_neighbors ||
      log_failure.size() != n_neighbors)
    stop("Incompatible Ward neighborhood vectors");
  if (mature.nrow() != active.nrow() || mature.ncol() != active.ncol())
    stop("Incompatible Ward landscape matrices");

  NumericVector result(n_targets);
  const int nr = active.nrow();
  const int nc = active.ncol();
  for (int j = 0; j < n_neighbors; ++j) {
    if (!(log_failure[j] < 0.0)) continue;
    const int dr = neighbor_row[j];
    const int dc = neighbor_col[j];
    int row_offset[4];
    int col_offset[4];
    if (dc == 0) {
      row_offset[0] = dr;  col_offset[0] = 0;
      row_offset[1] = 0;   col_offset[1] = dr;
      row_offset[2] = -dr; col_offset[2] = 0;
      row_offset[3] = 0;   col_offset[3] = -dr;
    } else {
      row_offset[0] = dr;  col_offset[0] = dc;
      row_offset[1] = -dr; col_offset[1] = dc;
      row_offset[2] = -dr; col_offset[2] = -dc;
      row_offset[3] = dr;  col_offset[3] = -dc;
    }
    for (int k = 0; k < n_targets; ++k) {
      const int target_row = targets(k, 0) - 1;
      const int target_col = targets(k, 1) - 1;
      int count = 0;
      for (int quadrant = 0; quadrant < 4; ++quadrant) {
        const int r = target_row + row_offset[quadrant];
        const int c = target_col + col_offset[quadrant];
        if (r >= 0 && r < nr && c >= 0 && c < nc &&
            mature(r, c) == TRUE && active(r, c) == TRUE)
          ++count;
      }
      const int power = ward2019 ? count : (count > 0 ? 1 : 0);
      if (power > 0) result[k] += power * log_failure[j];
    }
  }
  return result;
}

// Cohort properties needed by one synchronous recruitment step.
// [[Rcpp::export]]
List bs52_recruit_scan_cpp(List cohorts, IntegerVector ids,
                           CharacterVector species_names,
                           IntegerVector maturity, int timestep) {
  const int n_sites = ids.size();
  const int n_species = species_names.size();
  if (maturity.size() != n_species) stop("Incompatible maturity vector");
  std::unordered_map<std::string, int> species_index;
  for (int j = 0; j < n_species; ++j)
    species_index[as<std::string>(species_names[j])] = j;
  LogicalMatrix mature(n_sites, n_species);
  NumericVector shade_biomass(n_sites);
  NumericVector establishment_biomass(n_sites);
  for (int k = 0; k < n_sites; ++k) {
    DataFrame cohort = cohorts[ids[k] - 1];
    CharacterVector species = cohort["species"];
    NumericVector age = numeric_column(cohort, "age");
    NumericVector biomass = numeric_column(cohort, "biomass");
    for (int i = 0; i < age.size(); ++i) {
      const double b = biomass[i];
      if (age[i] > 5.0) shade_biomass[k] += b;
      if (age[i] >= timestep && age[i] > 1.0)
        establishment_biomass[k] += b;
      if (b > 0.0) {
        const auto found = species_index.find(as<std::string>(species[i]));
        if (found != species_index.end() && age[i] >= maturity[found->second])
          mature(k, found->second) = TRUE;
      }
    }
  }
  return List::create(
    _["mature"] = mature,
    _["shade_biomass"] = shade_biomass,
    _["establishment_biomass"] = establishment_biomass
  );
}

// Shade class and sufficient-light probability, before seed arrival.
// [[Rcpp::export]]
List bs52_shade_light_cpp(NumericVector shade_biomass,
                          NumericVector eco_max,
                          NumericMatrix thresholds,
                          NumericMatrix sep,
                          IntegerVector shade_tolerance,
                          NumericMatrix light) {
  const int n_sites = shade_biomass.size();
  const int n_species = sep.ncol();
  if (eco_max.size() != n_sites || thresholds.nrow() != n_sites ||
      sep.nrow() != n_sites || shade_tolerance.size() != n_species)
    stop("Incompatible shade inputs");
  IntegerVector shade(n_sites);
  NumericMatrix probability(n_sites, n_species);
  for (int k = 0; k < n_sites; ++k) {
    const double actual = std::min(eco_max[k], shade_biomass[k]);
    const double relative = actual / eco_max[k];
    int cls = 0;
    for (int h = 0; h < thresholds.ncol(); ++h)
      if (relative >= thresholds(k, h)) cls = h + 1;
    shade[k] = cls;
    for (int j = 0; j < n_species; ++j) {
      const int tolerance = shade_tolerance[j];
      if (tolerance < 1 || tolerance > light.nrow() || cls >= light.ncol())
        stop("Invalid shade-tolerance lookup");
      const double p = sep(k, j) * light(tolerance - 1, cls);
      probability(k, j) = std::min(1.0, std::max(0.0, p));
    }
  }
  return List::create(_["shade"] = shade, _["probability"] = probability);
}

// Adds sexual recruits in the exact species order supplied by R.
// [[Rcpp::export]]
List bs52_append_recruits_cpp(List cohorts, IntegerVector ids,
                              LogicalMatrix recruited,
                              CharacterVector species_names,
                              NumericMatrix max_anpp,
                              NumericVector eco_max,
                              NumericVector establishment_biomass) {
  const int n_sites = ids.size();
  const int n_species = species_names.size();
  if (recruited.nrow() != n_sites || recruited.ncol() != n_species ||
      max_anpp.nrow() != n_sites || max_anpp.ncol() != n_species ||
      eco_max.size() != n_sites || establishment_biomass.size() != n_sites)
    stop("Incompatible recruitment inputs");
  List out = clone(cohorts);
  for (int k = 0; k < n_sites; ++k) {
    DataFrame cohort = cohorts[ids[k] - 1];
    CharacterVector old_species = cohort["species"];
    NumericVector old_age = numeric_column(cohort, "age");
    NumericVector old_biomass = numeric_column(cohort, "biomass");
    int additions = 0;
    for (int j = 0; j < n_species; ++j)
      if (recruited(k, j) == TRUE) ++additions;
    if (!additions) continue;
    CharacterVector species(old_species.size() + additions);
    NumericVector age(old_age.size() + additions);
    NumericVector biomass(old_biomass.size() + additions);
    for (int i = 0; i < old_age.size(); ++i) {
      species[i] = old_species[i];
      age[i] = old_age[i];
      biomass[i] = old_biomass[i];
    }
    int next = old_age.size();
    for (int j = 0; j < n_species; ++j) {
      if (recruited(k, j) != TRUE) continue;
      const double anpp = max_anpp(k, j);
      const double candidate = anpp * std::exp(
        -1.6 * establishment_biomass[k] / eco_max[k]);
      species[next] = species_names[j];
      age[next] = 1.0;
      biomass[next] = std::trunc(std::max(2.0, std::min(anpp, candidate)));
      ++next;
    }
    out[ids[k] - 1] = DataFrame::create(
      _["species"] = species, _["age"] = age, _["biomass"] = biomass
    );
  }
  return out;
}

// Output table in cell-major then species-major order, matching do.call(rbind,
// records) in the R implementation.
// [[Rcpp::export]]
DataFrame bs52_recruit_records_cpp(CharacterVector species_names,
                                   LogicalMatrix mature,
                                   NumericMatrix arrival,
                                   NumericMatrix probability,
                                   LogicalMatrix recruited,
                                   IntegerVector ids,
                                   IntegerVector rows,
                                   IntegerVector cols, int year) {
  const int n_sites = ids.size();
  const int n_species = species_names.size();
  const int n = n_sites * n_species;
  if (mature.nrow() != n_sites || mature.ncol() != n_species ||
      arrival.nrow() != n_sites || arrival.ncol() != n_species ||
      probability.nrow() != n_sites || probability.ncol() != n_species ||
      recruited.nrow() != n_sites || recruited.ncol() != n_species ||
      rows.size() != n_sites || cols.size() != n_sites)
    stop("Incompatible recruitment-output inputs");
  CharacterVector species(n);
  LogicalVector mature_out(n), recruited_out(n);
  NumericVector arrival_out(n), probability_out(n);
  IntegerVector cell_id(n), row(n), col(n), year_out(n, year);
  int index = 0;
  for (int k = 0; k < n_sites; ++k) {
    for (int j = 0; j < n_species; ++j, ++index) {
      species[index] = species_names[j];
      mature_out[index] = mature(k, j);
      arrival_out[index] = arrival(k, j);
      probability_out[index] = probability(k, j);
      recruited_out[index] = recruited(k, j);
      cell_id[index] = ids[k];
      row[index] = rows[k];
      col[index] = cols[k];
    }
  }
  return DataFrame::create(
    _["species"] = species, _["mature"] = mature_out,
    _["seed_probability"] = arrival_out, _["probability"] = probability_out,
    _["recruited"] = recruited_out, _["cell_id"] = cell_id,
    _["row"] = row, _["col"] = col, _["year"] = year_out
  );
}


// ============================================================================
// Large-landscape bulk kernels.
//
// These functions preserve the validated single-state equations above while
// amortizing the R <-> C++ transition over all unique landscape states.
// They are used only by large52 v2 and do not change the reference API.
// ============================================================================

// [[Rcpp::export]]
List bs52_grow_states_bulk_cpp(List states, IntegerVector state_eco,
                               List params_by_eco, int timestep,
                               bool combine_last = true,
                               double spinup_fraction = 0.0,
                               bool spinup = false) {
  const int n_states = states.size();
  if (state_eco.size() != n_states)
    stop("Incompatible state/ecoregion vectors");

  List out(n_states);
  for (int s = 0; s < n_states; ++s) {
    if ((s % 10000) == 0) checkUserInterrupt();

    const int eco = state_eco[s];
    if (eco < 1 || eco > params_by_eco.size())
      stop("Invalid state ecoregion");

    SEXP block_sexp = params_by_eco[eco - 1];
    if (Rf_isNull(block_sexp))
      stop("Missing block parameters for active ecoregion");
    List block(block_sexp);

    DataFrame cohorts = states[s];
    CohortState state;
    state.species = as<std::vector<std::string> >(cohorts["species"]);
    state.age = as<std::vector<double> >(numeric_column(cohorts, "age"));
    state.biomass = as<std::vector<double> >(numeric_column(cohorts, "biomass"));

    const int n_years = block.size();
    for (int i = 0; i < n_years; ++i) {
      grow_year(state, as<DataFrame>(block[i]), timestep,
                combine_last && i == n_years - 1,
                spinup_fraction, spinup);
    }

    out[s] = as_data_frame(state);
  }

  return out;
}


// [[Rcpp::export]]
List bs52_grow_states_bulk_indexed_cpp(
    List states,
    IntegerVector state_eco,
    List params_by_eco,
    CharacterVector species_names,
    int timestep,
    bool combine_last = true,
    double spinup_fraction = 0.0,
    bool spinup = false,
    bool show_progress = false,
    std::string progress_label = "Croissance") {

  const int n_states = states.size();
  const int n_species = species_names.size();
  const int n_eco = params_by_eco.size();

  if (state_eco.size() != n_states)
    stop("Incompatible state/ecoregion vectors");
  if (n_species < 1)
    stop("No species supplied to indexed bulk growth");

  std::unordered_map<std::string, int> species_index;
  species_index.reserve(
    static_cast<std::size_t>(n_species * 2)
  );

  for (int j = 0; j < n_species; ++j) {
    const std::string name = as<std::string>(species_names[j]);

    if (species_index.find(name) != species_index.end())
      stop("Duplicate species name in indexed bulk growth");

    species_index[name] = j;
  }

  // Compile all annual eco/species parameter tables once for this block.
  std::vector<std::vector<GrowthParametersIndexed> > compiled(
    static_cast<std::size_t>(n_eco)
  );

  int n_years = -1;

  for (int eco = 0; eco < n_eco; ++eco) {
    SEXP block_sexp = params_by_eco[eco];

    if (Rf_isNull(block_sexp))
      continue;

    List block(block_sexp);

    if (n_years < 0) {
      n_years = block.size();
    } else if (block.size() != n_years) {
      stop("Inconsistent indexed growth block lengths");
    }

    compiled[static_cast<std::size_t>(eco)].reserve(
      static_cast<std::size_t>(block.size())
    );

    for (int i = 0; i < block.size(); ++i) {
      compiled[static_cast<std::size_t>(eco)].push_back(
        compile_growth_parameters_indexed(
          as<DataFrame>(block[i]),
          species_index,
          n_species
        )
      );
    }
  }

  if (n_years < 0)
    stop("No active parameter block for indexed growth");

  List out(n_states);
  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_states,
    progress_percent,
    show_progress
  );

  for (int s = 0; s < n_states; ++s) {
    if ((s % 10000) == 0)
      checkUserInterrupt();

    const int eco_1 = state_eco[s];

    if (eco_1 < 1 || eco_1 > n_eco)
      stop("Invalid state ecoregion in indexed bulk growth");

    const int eco = eco_1 - 1;

    if (compiled[static_cast<std::size_t>(eco)].size() !=
        static_cast<std::size_t>(n_years))
      stop("Missing indexed growth parameters for active ecoregion");

    DataFrame cohorts = states[s];

    CharacterVector cohort_species = cohorts["species"];
    NumericVector cohort_age = numeric_column(
      cohorts,
      "age"
    );
    NumericVector cohort_biomass = numeric_column(
      cohorts,
      "biomass"
    );

    if (cohort_species.size() != cohort_age.size() ||
        cohort_species.size() != cohort_biomass.size())
      stop("Incompatible cohort columns in indexed bulk growth");

    CohortStateIndexed state;
    state.species.reserve(
      static_cast<std::size_t>(cohort_species.size())
    );
    state.age.reserve(
      static_cast<std::size_t>(cohort_age.size())
    );
    state.biomass.reserve(
      static_cast<std::size_t>(cohort_biomass.size())
    );

    for (int i = 0; i < cohort_species.size(); ++i) {
      const std::string name = as<std::string>(
        cohort_species[i]
      );

      const auto found = species_index.find(name);

      if (found == species_index.end())
        stop("Unknown cohort species in indexed bulk growth");

      state.species.push_back(found->second);
      state.age.push_back(cohort_age[i]);
      state.biomass.push_back(cohort_biomass[i]);
    }

    for (int i = 0; i < n_years; ++i) {
      grow_year_indexed(
        state,
        compiled[static_cast<std::size_t>(eco)][
          static_cast<std::size_t>(i)
        ],
        timestep,
        combine_last && i == n_years - 1,
        spinup_fraction,
        spinup
      );
    }

    out[s] = as_data_frame_indexed(
      state,
      species_names
    );

    bs52_progress_update(
      progress_label,
      s + 1,
      n_states,
      progress_percent,
      show_progress
    );
  }

  return out;
}


// [[Rcpp::export]]
List bs52_grow_states_bulk_indexed_features_cpp(
    List states,
    IntegerVector state_eco,
    List params_by_eco,
    CharacterVector species_names,
    NumericVector maturity,
    int timestep,
    NumericMatrix eco_thresholds,
    NumericMatrix eco_sep,
    NumericVector eco_max,
    IntegerVector shade_tolerance,
    NumericMatrix light,
    bool combine_last = true,
    double spinup_fraction = 0.0,
    bool spinup = false,
    bool show_progress = false,
    std::string progress_label = "Croissance + features") {

  const int n_states = states.size();
  const int n_species = species_names.size();
  const int n_eco = params_by_eco.size();

  if (state_eco.size() != n_states)
    stop("Incompatible state/ecoregion vectors");

  if (n_species < 1)
    stop("No species supplied to fused indexed growth");

  if (maturity.size() != n_species ||
      eco_sep.ncol() != n_species ||
      eco_thresholds.nrow() != n_eco ||
      eco_sep.nrow() != n_eco ||
      eco_max.size() != n_eco ||
      shade_tolerance.size() != n_species)
    stop("Incompatible fused growth/feature inputs");

  // One name -> integer registry for the whole 10-y block.
  std::unordered_map<std::string, int> species_index;
  species_index.reserve(
    static_cast<std::size_t>(n_species * 2)
  );

  for (int j = 0; j < n_species; ++j) {
    const std::string name = as<std::string>(species_names[j]);

    if (species_index.find(name) != species_index.end())
      stop("Duplicate species name in fused indexed growth");

    species_index[name] = j;
  }

  // Compile all annual eco/species growth tables exactly once.
  std::vector<std::vector<GrowthParametersIndexed> > compiled(
    static_cast<std::size_t>(n_eco)
  );

  int n_years = -1;

  for (int eco = 0; eco < n_eco; ++eco) {
    SEXP block_sexp = params_by_eco[eco];

    if (Rf_isNull(block_sexp))
      continue;

    List block(block_sexp);

    if (n_years < 0) {
      n_years = block.size();
    } else if (block.size() != n_years) {
      stop("Inconsistent fused indexed growth block lengths");
    }

    compiled[static_cast<std::size_t>(eco)].reserve(
      static_cast<std::size_t>(block.size())
    );

    for (int i = 0; i < block.size(); ++i) {
      compiled[static_cast<std::size_t>(eco)].push_back(
        compile_growth_parameters_indexed(
          as<DataFrame>(block[i]),
          species_index,
          n_species
        )
      );
    }
  }

  if (n_years < 0)
    stop("No active parameter block for fused indexed growth");

  // Precompute SEP x light for every eco x shade class x species.
  const int n_shade_classes = light.ncol();

  std::vector<std::vector<double> > base_probability_lookup(
    static_cast<std::size_t>(n_eco * n_shade_classes),
    std::vector<double>(
      static_cast<std::size_t>(n_species),
      0.0
    )
  );

  for (int eco = 0; eco < n_eco; ++eco) {
    if (!std::isfinite(eco_max[eco]) || eco_max[eco] <= 0.0)
      continue;

    for (int cls = 0; cls < n_shade_classes; ++cls) {
      std::vector<double>& row =
        base_probability_lookup[
          static_cast<std::size_t>(
            eco * n_shade_classes + cls
          )
        ];

      for (int j = 0; j < n_species; ++j) {
        const int tolerance = shade_tolerance[j];

        if (tolerance < 1 || tolerance > light.nrow())
          stop("Invalid shade tolerance in fused growth");

        const double p =
          eco_sep(eco, j) *
          light(tolerance - 1, cls);

        row[static_cast<std::size_t>(j)] =
          std::min(
            1.0,
            std::max(
              0.0,
              p
            )
          );
      }
    }
  }

  List out_states(n_states);

  LogicalMatrix mature(
    n_states,
    n_species
  );

  NumericVector establishment_biomass(
    n_states
  );

  NumericMatrix base_probability(
    n_states,
    n_species
  );

  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_states,
    progress_percent,
    show_progress
  );

  for (int s = 0; s < n_states; ++s) {
    if ((s % 10000) == 0)
      checkUserInterrupt();

    const int eco_1 = state_eco[s];

    if (eco_1 < 1 || eco_1 > n_eco)
      stop("Invalid state ecoregion in fused indexed growth");

    const int eco = eco_1 - 1;

    if (compiled[static_cast<std::size_t>(eco)].size() !=
        static_cast<std::size_t>(n_years))
      stop("Missing fused growth parameters for active ecoregion");

    if (!std::isfinite(eco_max[eco]) ||
        eco_max[eco] <= 0.0)
      stop("Invalid ecoregion biomass maximum in fused growth");

    DataFrame cohorts = states[s];

    CharacterVector cohort_species =
      cohorts["species"];

    NumericVector cohort_age =
      numeric_column(
        cohorts,
        "age"
      );

    NumericVector cohort_biomass =
      numeric_column(
        cohorts,
        "biomass"
      );

    if (cohort_species.size() != cohort_age.size() ||
        cohort_species.size() != cohort_biomass.size())
      stop("Incompatible cohort columns in fused indexed growth");

    CohortStateIndexed state;

    state.species.reserve(
      static_cast<std::size_t>(
        cohort_species.size()
      )
    );

    state.age.reserve(
      static_cast<std::size_t>(
        cohort_age.size()
      )
    );

    state.biomass.reserve(
      static_cast<std::size_t>(
        cohort_biomass.size()
      )
    );

    for (int i = 0; i < cohort_species.size(); ++i) {
      const std::string name =
        as<std::string>(
          cohort_species[i]
        );

      const auto found =
        species_index.find(name);

      if (found == species_index.end())
        stop("Unknown cohort species in fused indexed growth");

      state.species.push_back(
        found->second
      );

      state.age.push_back(
        cohort_age[i]
      );

      state.biomass.push_back(
        cohort_biomass[i]
      );
    }

    // Exact same annual biological growth as v2.1.
    for (int i = 0; i < n_years; ++i) {
      grow_year_indexed(
        state,
        compiled[
          static_cast<std::size_t>(eco)
        ][
          static_cast<std::size_t>(i)
        ],
        timestep,
        combine_last &&
          i == n_years - 1,
        spinup_fraction,
        spinup
      );
    }

    // Compute all recruitment features immediately while the final
    // indexed cohort vectors are still hot in cache.
    double shade_biomass = 0.0;
    double establishment = 0.0;

    for (std::size_t i = 0;
         i < state.age.size();
         ++i) {

      const int species_id =
        state.species[i];

      const double age =
        state.age[i];

      const double biomass =
        state.biomass[i];

      if (species_id < 0 ||
          species_id >= n_species)
        stop("Invalid indexed species after fused growth");

      if (age > 5.0)
        shade_biomass += biomass;

      if (age >= timestep &&
          age > 1.0)
        establishment += biomass;

      if (biomass > 0.0 &&
          age >= maturity[species_id]) {
        mature(
          s,
          species_id
        ) = true;
      }
    }

    establishment_biomass[s] =
      establishment;

    const double actual =
      std::min(
        eco_max[eco],
        shade_biomass
      );

    const double relative =
      actual /
      eco_max[eco];

    int cls = 0;

    for (int h = 0;
         h < eco_thresholds.ncol();
         ++h) {

      const double threshold =
        eco_thresholds(
          eco,
          h
        );

      if (!NumericVector::is_na(threshold) &&
          relative >= threshold) {
        cls = h + 1;
      }
    }

    if (cls < 0 ||
        cls >= n_shade_classes)
      stop("Invalid shade class in fused growth");

    const std::vector<double>& probability_row =
      base_probability_lookup[
        static_cast<std::size_t>(
          eco * n_shade_classes + cls
        )
      ];

    for (int j = 0;
         j < n_species;
         ++j) {

      base_probability(
        s,
        j
      ) = probability_row[
        static_cast<std::size_t>(j)
      ];
    }

    // Convert back to the existing state representation only once.
    out_states[s] =
      as_data_frame_indexed(
        state,
        species_names
      );

    bs52_progress_update(
      progress_label,
      s + 1,
      n_states,
      progress_percent,
      show_progress
    );
  }

  return List::create(
    _["states"] = out_states,
    _["mature"] = mature,
    _["establishment_biomass"] =
      establishment_biomass,
    _["base_probability"] =
      base_probability
  );
}

// [[Rcpp::export]]
List bs52_state_features_bulk_cpp(
    List states, IntegerVector state_eco,
    CharacterVector species_names, NumericVector maturity, int timestep,
    NumericMatrix eco_thresholds, NumericMatrix eco_sep,
    NumericVector eco_max, IntegerVector shade_tolerance,
    NumericMatrix light,
    bool show_progress = false,
    std::string progress_label = "Features") {

  const int n_states = states.size();
  const int n_species = species_names.size();
  const int n_eco = eco_sep.nrow();

  if (state_eco.size() != n_states || maturity.size() != n_species ||
      eco_sep.ncol() != n_species || eco_thresholds.nrow() != n_eco ||
      eco_max.size() != n_eco || shade_tolerance.size() != n_species)
    stop("Incompatible bulk feature inputs");

  std::unordered_map<std::string, int> species_index;
  species_index.reserve(n_species);
  for (int j = 0; j < n_species; ++j)
    species_index[as<std::string>(species_names[j])] = j;

  LogicalMatrix mature(n_states, n_species);
  NumericVector shade_biomass(n_states);
  NumericVector establishment_biomass(n_states);
  NumericVector max_age(n_states);
  IntegerVector shade_class(n_states);
  NumericMatrix base_probability(n_states, n_species);
  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_states,
    progress_percent,
    show_progress
  );

  for (int s = 0; s < n_states; ++s) {
    if ((s % 10000) == 0) checkUserInterrupt();

    const int eco = state_eco[s] - 1;
    if (eco < 0 || eco >= n_eco || !std::isfinite(eco_max[eco]) ||
        eco_max[eco] <= 0.0)
      stop("Invalid state ecoregion in bulk features");

    DataFrame cohort = states[s];
    CharacterVector species = cohort["species"];
    NumericVector age = numeric_column(cohort, "age");
    NumericVector biomass = numeric_column(cohort, "biomass");

    double oldest = 0.0;
    for (int i = 0; i < age.size(); ++i) {
      const double a = age[i];
      const double b = biomass[i];
      if (a > oldest) oldest = a;
      if (a > 5.0) shade_biomass[s] += b;
      if (a >= timestep && a > 1.0)
        establishment_biomass[s] += b;

      if (b > 0.0) {
        const auto found = species_index.find(as<std::string>(species[i]));
        if (found == species_index.end())
          stop("Unknown cohort species in bulk features");
        if (a >= maturity[found->second])
          mature(s, found->second) = TRUE;
      }
    }
    max_age[s] = oldest;

    const double actual = std::min(eco_max[eco], shade_biomass[s]);
    const double relative = actual / eco_max[eco];
    int cls = 0;
    for (int h = 0; h < eco_thresholds.ncol(); ++h) {
      const double threshold = eco_thresholds(eco, h);
      if (!NumericVector::is_na(threshold) && relative >= threshold)
        cls = h + 1;
    }
    shade_class[s] = cls;

    if (cls < 0 || cls >= light.ncol())
      stop("Invalid shade class in bulk features");

    for (int j = 0; j < n_species; ++j) {
      const int tolerance = shade_tolerance[j];
      if (tolerance < 1 || tolerance > light.nrow())
        stop("Invalid shade tolerance in bulk features");
      const double p = eco_sep(eco, j) * light(tolerance - 1, cls);
      base_probability(s, j) = std::min(1.0, std::max(0.0, p));
    }

    bs52_progress_update(
      progress_label,
      s + 1,
      n_states,
      progress_percent,
      show_progress
    );
  }

  return List::create(
    _["mature"] = mature,
    _["shade_biomass"] = shade_biomass,
    _["establishment_biomass"] = establishment_biomass,
    _["max_age"] = max_age,
    _["shade_class"] = shade_class,
    _["base_probability"] = base_probability
  );
}

// [[Rcpp::export]]
List bs52_rebuild_states_bulk_cpp(
    List states, IntegerVector state_eco, IntegerVector cell_state,
    IntegerVector recruit_mask, CharacterVector species_names,
    NumericMatrix max_anpp_by_eco, NumericVector eco_max_by_eco,
    NumericVector establishment_biomass,
    bool show_progress = false,
    std::string progress_label = "Rebuild") {

  const int n_states = states.size();
  const int n_sites = cell_state.size();
  const int n_species = species_names.size();
  const int n_eco = max_anpp_by_eco.nrow();

  if (state_eco.size() != n_states || recruit_mask.size() != n_sites ||
      max_anpp_by_eco.ncol() != n_species || eco_max_by_eco.size() != n_eco ||
      establishment_biomass.size() != n_states)
    stop("Incompatible bulk rebuild inputs");
  if (n_species < 1 || n_species > 30)
    stop("Bulk rebuild supports 1 to 30 species");

  List out_states(n_sites);
  IntegerVector out_eco(n_sites);
  IntegerVector next_cell_state(n_sites);

  std::unordered_map<std::uint64_t, int> key_to_state;
  key_to_state.reserve(static_cast<std::size_t>(n_sites * 1.05));

  int n_unique = 0;
  int n_recruited_cells = 0;
  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_sites,
    progress_percent,
    show_progress
  );

  for (int k = 0; k < n_sites; ++k) {
    if ((k % 20000) == 0) checkUserInterrupt();

    const int old_state_1 = cell_state[k];
    if (old_state_1 < 1 || old_state_1 > n_states)
      stop("Invalid cell state in bulk rebuild");
    const int old_state = old_state_1 - 1;

    const int mask = recruit_mask[k];
    if (mask < 0)
      stop("Invalid recruitment mask");
    if (mask != 0) ++n_recruited_cells;

    const std::uint64_t key =
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(old_state_1)) << 32) |
      static_cast<std::uint32_t>(mask);

    auto found = key_to_state.find(key);
    int next_id;

    if (found == key_to_state.end()) {
      next_id = n_unique;
      key_to_state.emplace(key, next_id);
      ++n_unique;

      const int eco = state_eco[old_state] - 1;
      if (eco < 0 || eco >= n_eco || !std::isfinite(eco_max_by_eco[eco]) ||
          eco_max_by_eco[eco] <= 0.0)
        stop("Invalid ecoregion in bulk rebuild");
      out_eco[next_id] = eco + 1;

      if (mask == 0) {
        out_states[next_id] = states[old_state];
      } else {
        DataFrame cohort = states[old_state];
        CharacterVector old_species = cohort["species"];
        NumericVector old_age = numeric_column(cohort, "age");
        NumericVector old_biomass = numeric_column(cohort, "biomass");

        int additions = 0;
        for (int j = 0; j < n_species; ++j)
          if ((static_cast<unsigned int>(mask) & (1u << j)) != 0u)
            ++additions;

        CharacterVector species(old_species.size() + additions);
        NumericVector age(old_age.size() + additions);
        NumericVector biomass(old_biomass.size() + additions);

        for (int i = 0; i < old_age.size(); ++i) {
          species[i] = old_species[i];
          age[i] = old_age[i];
          biomass[i] = old_biomass[i];
        }

        int next = old_age.size();
        const double site_biomass = establishment_biomass[old_state];
        const double eco_max = eco_max_by_eco[eco];

        for (int j = 0; j < n_species; ++j) {
          if ((static_cast<unsigned int>(mask) & (1u << j)) == 0u)
            continue;
          const double anpp = max_anpp_by_eco(eco, j);
          const double candidate = anpp * std::exp(
            -1.6 * site_biomass / eco_max);
          species[next] = species_names[j];
          age[next] = 1.0;
          biomass[next] = std::trunc(
            std::max(2.0, std::min(anpp, candidate)));
          ++next;
        }

        out_states[next_id] = DataFrame::create(
          _["species"] = species,
          _["age"] = age,
          _["biomass"] = biomass
        );
      }
    } else {
      next_id = found->second;
    }

    next_cell_state[k] = next_id + 1;

    bs52_progress_update(
      progress_label,
      k + 1,
      n_sites,
      progress_percent,
      show_progress
    );
  }

  List trimmed_states(n_unique);
  IntegerVector trimmed_eco(n_unique);
  for (int i = 0; i < n_unique; ++i) {
    trimmed_states[i] = out_states[i];
    trimmed_eco[i] = out_eco[i];
  }

  return List::create(
    _["states"] = trimmed_states,
    _["state_eco"] = trimmed_eco,
    _["cell_state"] = next_cell_state,
    _["n_recruited_cells"] = n_recruited_cells
  );
}

// [[Rcpp::export]]
List bs52_state_snapshot_bulk_cpp(List states,
                                  CharacterVector species_names) {
  const int n_states = states.size();
  const int n_species = species_names.size();

  std::unordered_map<std::string, int> species_index;
  species_index.reserve(n_species);
  for (int j = 0; j < n_species; ++j)
    species_index[as<std::string>(species_names[j])] = j;

  NumericMatrix state_biomass(n_states, n_species);
  NumericVector state_total_biomass(n_states);
  NumericVector state_max_age(n_states);

  for (int s = 0; s < n_states; ++s) {
    if ((s % 10000) == 0) checkUserInterrupt();

    DataFrame cohort = states[s];
    CharacterVector species = cohort["species"];
    NumericVector age = numeric_column(cohort, "age");
    NumericVector biomass = numeric_column(cohort, "biomass");

    double oldest = 0.0;
    for (int i = 0; i < age.size(); ++i) {
      const auto found = species_index.find(as<std::string>(species[i]));
      if (found == species_index.end())
        stop("Unknown cohort species in bulk snapshot");
      state_biomass(s, found->second) += biomass[i];
      if (age[i] > oldest) oldest = age[i];
    }

    double total = 0.0;
    for (int j = 0; j < n_species; ++j)
      total += state_biomass(s, j);

    state_total_biomass[s] = total;
    state_max_age[s] = oldest;
  }

  colnames(state_biomass) = species_names;

  return List::create(
    _["state_biomass"] = state_biomass,
    _["state_total_biomass"] = state_total_biomass,
    _["state_max_age"] = state_max_age
  );
}

// ============================================================================
// large52 v3.4 exact age-mortality lookup
//
// The final mortality fraction is deliberately NOT cached because replacing
//     biomass * numerator / denominator
// with
//     biomass * (numerator / denominator)
// can change the last floating-point bits.
//
// Instead, the two exp() results are cached separately and the original
// multiply-then-divide sequence is retained at runtime.
// ============================================================================

struct MortalityAgeLookupIndexed {
  std::vector<unsigned char> enabled;
  std::vector<double> longevity;
  std::vector<double> mortality_shape;
  std::vector<double> shape_exp;
  std::vector<std::vector<double> > age_exp;
};


MortalityAgeLookupIndexed bs52_build_mortality_age_lookup(
    const std::vector<std::vector<GrowthParametersIndexed> >& compiled,
    int n_species) {

  MortalityAgeLookupIndexed out;

  out.enabled.assign(
    static_cast<std::size_t>(n_species),
    static_cast<unsigned char>(0)
  );

  out.longevity.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );

  out.mortality_shape.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );

  out.shape_exp.assign(
    static_cast<std::size_t>(n_species),
    NA_REAL
  );

  out.age_exp.resize(
    static_cast<std::size_t>(n_species)
  );

  std::vector<unsigned char> initialized(
    static_cast<std::size_t>(n_species),
    static_cast<unsigned char>(0)
  );

  // First pass: obtain one reference longevity/mortality-shape pair for each
  // species, then disable the lookup if any active eco/year differs.
  for (std::size_t eco = 0;
       eco < compiled.size();
       ++eco) {

    for (std::size_t year = 0;
         year < compiled[eco].size();
         ++year) {

      const GrowthParametersIndexed& params =
        compiled[eco][year];

      if (params.longevity.size() !=
            static_cast<std::size_t>(n_species) ||
          params.mortality_shape.size() !=
            static_cast<std::size_t>(n_species))
        stop("Incompatible parameters while building mortality lookup");

      for (int j = 0;
           j < n_species;
           ++j) {

        const std::size_t ju =
          static_cast<std::size_t>(j);

        const double longevity =
          params.longevity[ju];

        const double mortality_shape =
          params.mortality_shape[ju];

        if (!std::isfinite(longevity) ||
            longevity <= 0.0 ||
            !std::isfinite(mortality_shape)) {
          out.enabled[ju] = 0;
          initialized[ju] = 2;
          continue;
        }

        if (initialized[ju] == 0) {
          out.longevity[ju] =
            longevity;

          out.mortality_shape[ju] =
            mortality_shape;

          out.enabled[ju] = 1;
          initialized[ju] = 1;
        } else if (
          initialized[ju] == 1 &&
          (
            longevity !=
              out.longevity[ju] ||
            mortality_shape !=
              out.mortality_shape[ju]
          )
        ) {
          out.enabled[ju] = 0;
          initialized[ju] = 2;
        }
      }
    }
  }

  // Second pass: one compact age table per eligible species.
  for (int j = 0;
       j < n_species;
       ++j) {

    const std::size_t ju =
      static_cast<std::size_t>(j);

    if (!out.enabled[ju])
      continue;

    out.shape_exp[ju] =
      std::exp(
        out.mortality_shape[ju]
      );

    if (!std::isfinite(
          out.shape_exp[ju]
        ) ||
        out.shape_exp[ju] == 0.0) {
      out.enabled[ju] = 0;
      out.age_exp[ju].clear();
      continue;
    }

    // A living cohort is processed only while its old age is < longevity.
    // updated_age can therefore reach ceil(longevity). Two extra slots keep
    // the lookup robust to integer-valued edge cases.
    const double max_age_double =
      std::ceil(
        out.longevity[ju]
      ) +
      2.0;

    if (!std::isfinite(
          max_age_double
        ) ||
        max_age_double < 0.0 ||
        max_age_double >
          static_cast<double>(
            std::numeric_limits<int>::max() - 1
          )) {
      out.enabled[ju] = 0;
      out.age_exp[ju].clear();
      continue;
    }

    const int max_age =
      static_cast<int>(
        max_age_double
      );

    out.age_exp[ju].resize(
      static_cast<std::size_t>(
        max_age + 1
      )
    );

    for (int age = 0;
         age <= max_age;
         ++age) {

      // Keep the exact subexpression order used in cohort_change():
      // age / longevity * mortality_shape
      out.age_exp[ju][
        static_cast<std::size_t>(age)
      ] =
        std::exp(
          static_cast<double>(age) /
            out.longevity[ju] *
            out.mortality_shape[ju]
        );
    }
  }

  return out;
}


inline bool bs52_mortality_age_from_lookup(
    double biomass,
    double age,
    int species_id,
    const MortalityAgeLookupIndexed& lookup,
    double& mortality_age) {

  if (species_id < 0 ||
      species_id >=
        static_cast<int>(
          lookup.enabled.size()
        ))
    return false;

  const std::size_t su =
    static_cast<std::size_t>(
      species_id
    );

  if (!lookup.enabled[su])
    return false;

  if (!std::isfinite(age) ||
      age < 0.0 ||
      age >
        static_cast<double>(
          std::numeric_limits<int>::max()
        ))
    return false;

  const int age_index =
    static_cast<int>(
      age
    );

  if (static_cast<double>(
        age_index
      ) != age)
    return false;

  if (age_index < 0 ||
      static_cast<std::size_t>(
        age_index
      ) >=
        lookup.age_exp[su].size())
    return false;

  const double numerator =
    lookup.age_exp[su][
      static_cast<std::size_t>(
        age_index
      )
    ];

  const double denominator =
    lookup.shape_exp[su];

  // IMPORTANT: preserve the original arithmetic order:
  //   (biomass * exp(...)) / exp(shape)
  mortality_age =
    std::min(
      biomass *
        numerator /
        denominator,
      biomass
    );

  return true;
}


inline CohortChange cohort_change_with_age_lookup(
    double biomass,
    double age,
    int species_id,
    const GrowthParametersIndexed& params,
    double site_biomass,
    double competition,
    bool spinup,
    double spinup_mortality_fraction,
    const MortalityAgeLookupIndexed& lookup) {

  if (species_id < 0 ||
      species_id >=
        static_cast<int>(
          params.longevity.size()
        ))
    stop("Invalid species in cached cohort change");

  const std::size_t su =
    static_cast<std::size_t>(
      species_id
    );

  double mortality_age =
    0.0;

  if (!bs52_mortality_age_from_lookup(
        biomass,
        age,
        species_id,
        lookup,
        mortality_age
      )) {

    // Exact fallback to the previously validated implementation.
    return cohort_change(
      biomass,
      age,
      params.longevity[su],
      params.mortality_shape[su],
      params.growth_shape[su],
      params.max_anpp[su],
      params.max_biomass[su],
      site_biomass,
      competition,
      spinup,
      spinup_mortality_fraction
    );
  }

  if (spinup &&
      spinup_mortality_fraction > 0.0) {

    mortality_age =
      std::min(
        mortality_age +
          biomass *
          spinup_mortality_fraction,
        biomass
      );
  }

  const double potential_biomass =
    std::max(
      0.0,
      std::max(
        1.0,
        params.max_biomass[su] -
          site_biomass +
          biomass
      )
    );

  const double b_ap =
    biomass /
    potential_biomass;

  const double power =
    R_pow(
      b_ap,
      params.growth_shape[su]
    );

  double raw_anpp =
    params.max_anpp[su] *
    std::exp(1.0) *
    power *
    std::exp(-power) *
    competition;

  raw_anpp =
    std::min(
      params.max_anpp[su] *
        competition,
      raw_anpp
    );

  const double actual_anpp =
    std::max(
      1.0,
      raw_anpp -
        mortality_age
    );

  double mortality_growth =
    b_ap > 1.0 ?
      params.max_anpp[su] *
        competition :
      params.max_anpp[su] *
        (2.0 * b_ap) /
        (1.0 + b_ap) *
        competition;

  mortality_growth =
    std::min(
      biomass,
      std::min(
        mortality_growth,
        params.max_anpp[su] *
          competition
      )
    );

  mortality_growth =
    std::min(
      std::max(
        0.0,
        mortality_growth -
          mortality_age
      ),
      actual_anpp
    );

  CohortChange out;

  out.new_biomass =
    biomass +
    std::trunc(
      actual_anpp -
        mortality_age -
        mortality_growth
    );

  out.raw_anpp =
    raw_anpp;

  out.actual_anpp =
    actual_anpp;

  out.mortality_age =
    mortality_age;

  out.total_mortality =
    mortality_age +
    mortality_growth;

  return out;
}


// ============================================================================
// large52 v3.4 structural growth kernel
//
// This path is production-CSR specific. It keeps the validated biological
// equations and processing order but removes repeated small heap allocations
// and redundant regrouping between annual steps.
// ============================================================================


struct CompetitionWeightLookup52 {
  std::vector<double> value;
  mutable long long fallback_count;

  explicit CompetitionWeightLookup52(
      int max_integer_biomass = 65535) :
    value(
      static_cast<std::size_t>(
        max_integer_biomass + 1
      )
    ),
    fallback_count(0) {

    for (int i = 0;
         i <= max_integer_biomass;
         ++i) {

      value[
        static_cast<std::size_t>(i)
      ] =
        std::max(
          R_pow(
            static_cast<double>(i),
            0.95
          ),
          1.0
        );
    }
  }

  inline double get(
      double biomass) const {

    if (biomass >= 0.0 &&
        biomass <=
          static_cast<double>(
            value.size() - 1
          ) &&
        biomass ==
          std::trunc(
            biomass
          )) {

      return value[
        static_cast<std::size_t>(
          biomass
        )
      ];
    }

    ++fallback_count;

    return std::max(
      R_pow(
        biomass,
        0.95
      ),
      1.0
    );
  }
};

struct GrowthScratchIndexed {
  std::vector<unsigned char> seen_species;
  std::vector<int> species_order;

  std::vector<int> tmp_species;
  std::vector<double> tmp_age;
  std::vector<double> tmp_biomass;

  std::vector<int> run_begin;
  std::vector<int> run_end;
  std::vector<int> run_position;
  std::vector<int> iterators;

  std::vector<unsigned char> alive;
  std::vector<double> competition_weight;

  std::vector<double> young_total;
  std::vector<unsigned char> young_found;

  explicit GrowthScratchIndexed(int n_species = 0) {
    reset_species_capacity(n_species);
  }

  void reset_species_capacity(int n_species) {
    seen_species.assign(
      static_cast<std::size_t>(n_species),
      static_cast<unsigned char>(0)
    );

    young_total.assign(
      static_cast<std::size_t>(n_species),
      0.0
    );

    young_found.assign(
      static_cast<std::size_t>(n_species),
      static_cast<unsigned char>(0)
    );

    species_order.reserve(
      static_cast<std::size_t>(n_species)
    );

    run_begin.reserve(
      static_cast<std::size_t>(n_species)
    );

    run_end.reserve(
      static_cast<std::size_t>(n_species)
    );

    run_position.reserve(
      static_cast<std::size_t>(n_species)
    );

    iterators.reserve(
      static_cast<std::size_t>(n_species)
    );
  }
};


inline void bs52_species_order_reuse(
    const CohortStateIndexed& state,
    int n_species,
    GrowthScratchIndexed& scratch) {

  if (static_cast<int>(
        scratch.seen_species.size()
      ) != n_species) {
    scratch.reset_species_capacity(
      n_species
    );
  }

  std::fill(
    scratch.seen_species.begin(),
    scratch.seen_species.end(),
    static_cast<unsigned char>(0)
  );

  scratch.species_order.clear();

  for (std::size_t i = 0;
       i < state.species.size();
       ++i) {

    const int species_id =
      state.species[i];

    if (species_id < 0 ||
        species_id >= n_species)
      stop("Invalid indexed cohort species in v3.3 growth");

    if (!scratch.seen_species[
          static_cast<std::size_t>(
            species_id
          )
        ]) {

      scratch.seen_species[
        static_cast<std::size_t>(
          species_id
        )
      ] = 1;

      scratch.species_order.push_back(
        species_id
      );
    }
  }
}


inline void bs52_stable_group_reuse(
    CohortStateIndexed& state,
    const std::vector<int>& species_order,
    GrowthScratchIndexed& scratch) {

  const std::size_t n =
    state.species.size();

  scratch.tmp_species.clear();
  scratch.tmp_age.clear();
  scratch.tmp_biomass.clear();

  scratch.tmp_species.reserve(n);
  scratch.tmp_age.reserve(n);
  scratch.tmp_biomass.reserve(n);

  // Same nested-loop stable grouping as grow_year_indexed().
  for (const int species_id :
       species_order) {

    for (std::size_t i = 0;
         i < n;
         ++i) {

      if (state.species[i] ==
          species_id) {

        scratch.tmp_species.push_back(
          state.species[i]
        );

        scratch.tmp_age.push_back(
          state.age[i]
        );

        scratch.tmp_biomass.push_back(
          state.biomass[i]
        );
      }
    }
  }

  state.species.swap(
    scratch.tmp_species
  );

  state.age.swap(
    scratch.tmp_age
  );

  state.biomass.swap(
    scratch.tmp_biomass
  );
}


inline void bs52_combine_young_and_group_reuse(
    CohortStateIndexed& state,
    const std::vector<int>& species_order,
    int timestep,
    int n_species,
    GrowthScratchIndexed& scratch) {

  if (timestep <= 1) {
    bs52_stable_group_reuse(
      state,
      species_order,
      scratch
    );

    return;
  }

  if (static_cast<int>(
        scratch.young_total.size()
      ) != n_species) {
    scratch.reset_species_capacity(
      n_species
    );
  }

  std::fill(
    scratch.young_total.begin(),
    scratch.young_total.end(),
    0.0
  );

  std::fill(
    scratch.young_found.begin(),
    scratch.young_found.end(),
    static_cast<unsigned char>(0)
  );

  const std::size_t n =
    state.species.size();

  // Preserve the exact per-species biomass accumulation order from
  // combine_young_indexed(): for each species in first-occurrence order,
  // scan cohorts from first to last.
  for (const int species_id :
       species_order) {

    double total =
      0.0;

    bool found =
      false;

    for (std::size_t i = 0;
         i < n;
         ++i) {

      if (state.species[i] ==
            species_id &&
          state.age[i] <=
            timestep) {

        total +=
          state.biomass[i];

        found =
          true;
      }
    }

    if (found) {
      scratch.young_total[
        static_cast<std::size_t>(
          species_id
        )
      ] = total;

      scratch.young_found[
        static_cast<std::size_t>(
          species_id
        )
      ] = 1;
    }
  }

  scratch.tmp_species.clear();
  scratch.tmp_age.clear();
  scratch.tmp_biomass.clear();

  scratch.tmp_species.reserve(
    n +
    species_order.size()
  );

  scratch.tmp_age.reserve(
    n +
    species_order.size()
  );

  scratch.tmp_biomass.reserve(
    n +
    species_order.size()
  );

  // This directly produces the result of:
  //   combine_young_indexed(...)
  //   stable grouping by species_order
  //
  // For each species: surviving old cohorts retain their original order,
  // followed by the one combined young cohort.
  for (const int species_id :
       species_order) {

    for (std::size_t i = 0;
         i < n;
         ++i) {

      if (state.species[i] ==
            species_id &&
          state.age[i] >
            timestep) {

        scratch.tmp_species.push_back(
          state.species[i]
        );

        scratch.tmp_age.push_back(
          state.age[i]
        );

        scratch.tmp_biomass.push_back(
          state.biomass[i]
        );
      }
    }

    if (scratch.young_found[
          static_cast<std::size_t>(
            species_id
          )
        ]) {

      scratch.tmp_species.push_back(
        species_id
      );

      scratch.tmp_age.push_back(
        timestep - 1.0
      );

      scratch.tmp_biomass.push_back(
        scratch.young_total[
          static_cast<std::size_t>(
            species_id
          )
        ]
      );
    }
  }

  state.species.swap(
    scratch.tmp_species
  );

  state.age.swap(
    scratch.tmp_age
  );

  state.biomass.swap(
    scratch.tmp_biomass
  );
}


inline void bs52_build_group_runs_reuse(
    const CohortStateIndexed& state,
    GrowthScratchIndexed& scratch) {

  scratch.run_begin.clear();
  scratch.run_end.clear();
  scratch.run_position.clear();

  if (state.species.empty())
    return;

  int begin =
    0;

  int current =
    state.species[0];

  for (int i = 1;
       i < static_cast<int>(
         state.species.size()
       );
       ++i) {

    if (state.species[
          static_cast<std::size_t>(i)
        ] != current) {

      scratch.run_begin.push_back(
        begin
      );

      scratch.run_end.push_back(
        i
      );

      scratch.run_position.push_back(
        begin
      );

      begin =
        i;

      current =
        state.species[
          static_cast<std::size_t>(i)
        ];
    }
  }

  scratch.run_begin.push_back(
    begin
  );

  scratch.run_end.push_back(
    static_cast<int>(
      state.species.size()
    )
  );

  scratch.run_position.push_back(
    begin
  );
}


inline void bs52_compact_alive_in_place(
    CohortStateIndexed& state,
    const std::vector<unsigned char>& alive) {

  std::size_t write =
    0;

  for (std::size_t read = 0;
       read < state.species.size();
       ++read) {

    if (!alive[read])
      continue;

    if (write != read) {
      state.species[write] =
        state.species[read];

      state.age[write] =
        state.age[read];

      state.biomass[write] =
        state.biomass[read];
    }

    ++write;
  }

  state.species.resize(
    write
  );

  state.age.resize(
    write
  );

  state.biomass.resize(
    write
  );
}


void grow_year_indexed_structural(
    CohortStateIndexed& state,
    const GrowthParametersIndexed& params,
    int timestep,
    bool do_combine,
    double spinup_fraction,
    bool spinup,
    bool& state_is_grouped,
    GrowthScratchIndexed& scratch,
    const MortalityAgeLookupIndexed& mortality_lookup,
    const CompetitionWeightLookup52& competition_weight_lookup) {

  if (state.age.empty()) {
    state_is_grouped =
      true;

    return;
  }

  const int n_species =
    static_cast<int>(
      params.longevity.size()
    );

  if (params.mortality_shape.size() !=
        params.longevity.size() ||
      params.growth_shape.size() !=
        params.longevity.size() ||
      params.max_anpp.size() !=
        params.longevity.size() ||
      params.max_biomass.size() !=
        params.longevity.size())
    stop("Incompatible indexed growth parameters in v3.3");

  // Exactly the same site-biomass accumulation order as the validated path:
  // calculate BEFORE any grouping or young-cohort combination.
  double site_biomass =
    0.0;

  for (std::size_t i = 0;
       i < state.age.size();
       ++i) {

    if (state.age[i] >=
        timestep) {
      site_biomass +=
        state.biomass[i];
    }
  }

  // Species order is still first occurrence in the current state, matching
  // unique_species_indexed().
  bs52_species_order_reuse(
    state,
    n_species,
    scratch
  );

  if (do_combine) {
    bs52_combine_young_and_group_reuse(
      state,
      scratch.species_order,
      timestep,
      n_species,
      scratch
    );

    state_is_grouped =
      true;
  } else if (!state_is_grouped) {
    bs52_stable_group_reuse(
      state,
      scratch.species_order,
      scratch
    );

    state_is_grouped =
      true;
  }

  if (state.age.empty())
    return;

  // Because the state is species-grouped, each old vector<int> queue is now
  // exactly one contiguous [begin,end) run.
  bs52_build_group_runs_reuse(
    state,
    scratch
  );

  const int n_runs =
    static_cast<int>(
      scratch.run_begin.size()
    );

  scratch.alive.assign(
    state.age.size(),
    static_cast<unsigned char>(1)
  );

  scratch.competition_weight.resize(
    state.biomass.size()
  );

  for (std::size_t i = 0;
       i < state.biomass.size();
       ++i) {

    scratch.competition_weight[i] =
      competition_weight_lookup.get(
        state.biomass[i]
      );
  }

  auto head_age =
    [&](int run) {

      return state.age[
        static_cast<std::size_t>(
          scratch.run_position[
            static_cast<std::size_t>(
              run
            )
          ]
        )
      ];
    };

  scratch.iterators.clear();

  for (int run = 0;
       run < n_runs;
       ++run) {

    auto at =
      scratch.iterators.end();

    for (auto it =
           scratch.iterators.begin();
         it !=
           scratch.iterators.end();
         ++it) {

      if (head_age(*it) <
          head_age(run)) {

        at =
          it;

        break;
      }
    }

    scratch.iterators.insert(
      at,
      run
    );
  }

  while (!scratch.iterators.empty()) {
    const int run =
      scratch.iterators.front();

    const std::size_t ru =
      static_cast<std::size_t>(
        run
      );

    const int id =
      scratch.run_position[ru];

    const int species_id =
      state.species[
        static_cast<std::size_t>(
          id
        )
      ];

    if (species_id < 0 ||
        species_id >= n_species)
      stop("Invalid indexed species in v3.3 growth");

    if (state.age[
          static_cast<std::size_t>(
            id
          )
        ] >=
        params.longevity[
          static_cast<std::size_t>(
            species_id
          )
        ]) {

      scratch.alive[
        static_cast<std::size_t>(
          id
        )
      ] = 0;
    } else {
      const double updated_age =
        state.age[
          static_cast<std::size_t>(
            id
          )
        ] +
        1.0;

      const double focal_weight =
        scratch.competition_weight[
          static_cast<std::size_t>(
            id
          )
        ];

      long double
        other_weight_accumulator =
          0.0L;

      // Preserve the exact old index-order accumulation.
      for (std::size_t j = 0;
           j < state.age.size();
           ++j) {

        if (!scratch.alive[j])
          continue;

        const bool include =
          state.age[j] + 1.0 !=
            updated_age ||
          state.species[j] !=
            species_id;

        if (include) {
          other_weight_accumulator +=
            scratch.competition_weight[j];
        }
      }

      const double other_weight =
        static_cast<double>(
          other_weight_accumulator
        );

      const double competition =
        focal_weight /
        (
          focal_weight +
          other_weight
        );

      const CohortChange change =
        cohort_change_with_age_lookup(
          state.biomass[
            static_cast<std::size_t>(
              id
            )
          ],
          updated_age,
          species_id,
          params,
          site_biomass,
          competition,
          spinup,
          spinup_fraction,
          mortality_lookup
        );

      state.age[
        static_cast<std::size_t>(
          id
        )
      ] =
        updated_age;

      state.biomass[
        static_cast<std::size_t>(
          id
        )
      ] =
        std::max(
          0.0,
          change.new_biomass
        );

      if (state.biomass[
            static_cast<std::size_t>(
              id
            )
          ] == 0.0) {

        scratch.alive[
          static_cast<std::size_t>(
            id
          )
        ] = 0;

        scratch.competition_weight[
          static_cast<std::size_t>(
            id
          )
        ] =
          0.0;
      } else {
        scratch.competition_weight[
          static_cast<std::size_t>(
            id
          )
        ] =
          competition_weight_lookup.get(
            state.biomass[
              static_cast<std::size_t>(
                id
              )
            ]
          );
      }
    }

    ++scratch.run_position[ru];

    if (scratch.run_position[ru] >=
        scratch.run_end[ru]) {

      scratch.iterators.erase(
        scratch.iterators.begin()
      );
    } else if (
      scratch.iterators.size() > 1 &&
      head_age(run) <
        head_age(
          scratch.iterators[1]
        )
    ) {

      scratch.iterators.erase(
        scratch.iterators.begin()
      );

      auto at =
        scratch.iterators.end();

      for (auto it =
             scratch.iterators.begin();
           it !=
             scratch.iterators.end();
           ++it) {

        if (head_age(*it) <
            head_age(run)) {

          at =
            it;

          break;
        }
      }

      scratch.iterators.insert(
        at,
        run
      );
    }
  }

  bs52_compact_alive_in_place(
    state,
    scratch.alive
  );

  // In-place compaction cannot break species grouping.
  state_is_grouped =
    true;
}


// [[Rcpp::export]]
List bs52_grow_csr_indexed_features_cpp(
    List store,
    IntegerVector state_eco,
    List params_by_eco,
    CharacterVector species_names,
    NumericVector maturity,
    int timestep,
    NumericMatrix eco_thresholds,
    NumericMatrix eco_sep,
    NumericVector eco_max,
    IntegerVector shade_tolerance,
    NumericMatrix light,
    bool combine_last = true,
    double spinup_fraction = 0.0,
    bool spinup = false,
    bool show_progress = false,
    std::string progress_label = "Croissance CSR + features") {

  const int n_species =
    species_names.size();

  validate_csr_store(
    store,
    n_species
  );

  IntegerVector offset =
    store["offset"];

  IntegerVector input_species =
    store["species_id"];

  NumericVector input_age =
    store["age"];

  NumericVector input_biomass =
    store["biomass"];

  const int n_states =
    offset.size() - 1;

  const int n_eco =
    params_by_eco.size();

  if (state_eco.size() != n_states)
    stop("Incompatible CSR state/ecoregion vectors");

  if (maturity.size() != n_species ||
      eco_sep.ncol() != n_species ||
      eco_thresholds.nrow() != n_eco ||
      eco_sep.nrow() != n_eco ||
      eco_max.size() != n_eco ||
      shade_tolerance.size() != n_species)
    stop("Incompatible CSR growth/feature inputs");

  std::unordered_map<std::string, int>
    species_index;

  species_index.reserve(
    static_cast<std::size_t>(
      n_species * 2
    )
  );

  for (int j = 0; j < n_species; ++j) {
    const std::string name =
      as<std::string>(
        species_names[j]
      );

    if (species_index.find(name) !=
        species_index.end())
      stop("Duplicate species name in CSR growth");

    species_index[name] = j;
  }

  // Compile annual eco/species growth parameters once.
  std::vector<
    std::vector<GrowthParametersIndexed>
  > compiled(
    static_cast<std::size_t>(
      n_eco
    )
  );

  int n_years = -1;

  for (int eco = 0; eco < n_eco; ++eco) {
    SEXP block_sexp =
      params_by_eco[eco];

    if (Rf_isNull(block_sexp))
      continue;

    List block(
      block_sexp
    );

    if (n_years < 0) {
      n_years =
        block.size();
    } else if (
      block.size() != n_years
    ) {
      stop("Inconsistent CSR growth block lengths");
    }

    compiled[
      static_cast<std::size_t>(eco)
    ].reserve(
      static_cast<std::size_t>(
        block.size()
      )
    );

    for (int i = 0;
         i < block.size();
         ++i) {

      compiled[
        static_cast<std::size_t>(eco)
      ].push_back(
        compile_growth_parameters_indexed(
          as<DataFrame>(
            block[i]
          ),
          species_index,
          n_species
        )
      );
    }
  }

  if (n_years < 0)
    stop("No active parameter block for CSR growth");

  // v3.4: build one shared mortality-age lookup across the entire growth
  // block. Static mortality parameters are normally identical across
  // eco/year tables; any species that is not identical automatically falls
  // back to the original cohort_change() path.
  const MortalityAgeLookupIndexed mortality_lookup =
    bs52_build_mortality_age_lookup(
      compiled,
      n_species
    );

  // Exact lookup for max(R_pow(biomass, 0.95), 1.0).
  // Biomass in the validated LANDIS path is integer-valued.  Any unexpected
  // non-integer or >65535 value automatically falls back to the original
  // R_pow() calculation, preserving exact semantics.
  const CompetitionWeightLookup52 competition_weight_lookup(
    65535
  );

  const int n_shade_classes =
    light.ncol();

  std::vector<
    std::vector<double>
  > base_probability_lookup(
    static_cast<std::size_t>(
      n_eco *
      n_shade_classes
    ),
    std::vector<double>(
      static_cast<std::size_t>(
        n_species
      ),
      0.0
    )
  );

  for (int eco = 0;
       eco < n_eco;
       ++eco) {

    if (!std::isfinite(
          eco_max[eco]
        ) ||
        eco_max[eco] <= 0.0)
      continue;

    for (int cls = 0;
         cls < n_shade_classes;
         ++cls) {

      std::vector<double>& row =
        base_probability_lookup[
          static_cast<std::size_t>(
            eco *
            n_shade_classes +
            cls
          )
        ];

      for (int j = 0;
           j < n_species;
           ++j) {

        const int tolerance =
          shade_tolerance[j];

        if (tolerance < 1 ||
            tolerance > light.nrow())
          stop("Invalid shade tolerance in CSR growth");

        const double p =
          eco_sep(
            eco,
            j
          ) *
          light(
            tolerance - 1,
            cls
          );

        row[
          static_cast<std::size_t>(j)
        ] = std::min(
          1.0,
          std::max(
            0.0,
            p
          )
        );
      }
    }
  }

  // Growth never creates cohorts, so reserving the input cohort count
  // is sufficient for the output vectors.
  std::vector<int>
    out_species;

  std::vector<double>
    out_age;

  std::vector<double>
    out_biomass;

  out_species.reserve(
    static_cast<std::size_t>(
      input_species.size()
    )
  );

  out_age.reserve(
    static_cast<std::size_t>(
      input_age.size()
    )
  );

  out_biomass.reserve(
    static_cast<std::size_t>(
      input_biomass.size()
    )
  );

  IntegerVector out_offset(
    n_states + 1
  );

  out_offset[0] = 0;

  LogicalMatrix mature(
    n_states,
    n_species
  );

  NumericVector establishment_biomass(
    n_states
  );

  NumericMatrix base_probability(
    n_states,
    n_species
  );
  IntegerVector shade_class(n_states);

  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_states,
    progress_percent,
    show_progress
  );

  // v3.3: reuse these buffers across every state instead of allocating them
  // again for each state/year.
  CohortStateIndexed state;
  GrowthScratchIndexed growth_scratch(
    n_species
  );

  for (int s = 0;
       s < n_states;
       ++s) {

    if ((s % 10000) == 0)
      checkUserInterrupt();

    const int eco_1 =
      state_eco[s];

    if (eco_1 < 1 ||
        eco_1 > n_eco)
      stop("Invalid state ecoregion in CSR growth");

    const int eco =
      eco_1 - 1;

    if (compiled[
          static_cast<std::size_t>(eco)
        ].size() !=
        static_cast<std::size_t>(
          n_years
        ))
      stop("Missing CSR growth parameters for active ecoregion");

    if (!std::isfinite(
          eco_max[eco]
        ) ||
        eco_max[eco] <= 0.0)
      stop("Invalid ecoregion biomass maximum in CSR growth");

    const int begin =
      offset[s];

    const int end =
      offset[s + 1];

    state.species.clear();
    state.age.clear();
    state.biomass.clear();

    const std::size_t state_size =
      static_cast<std::size_t>(
        end - begin
      );

    if (state.species.capacity() <
        state_size)
      state.species.reserve(
        state_size
      );

    if (state.age.capacity() <
        state_size)
      state.age.reserve(
        state_size
      );

    if (state.biomass.capacity() <
        state_size)
      state.biomass.reserve(
        state_size
      );

    for (int k = begin;
         k < end;
         ++k) {

      state.species.push_back(
        input_species[k]
      );

      state.age.push_back(
        input_age[k]
      );

      state.biomass.push_back(
        input_biomass[k]
      );
    }

    bool state_is_grouped =
      false;

    for (int i = 0;
         i < n_years;
         ++i) {

      grow_year_indexed_structural(
        state,
        compiled[
          static_cast<std::size_t>(eco)
        ][
          static_cast<std::size_t>(i)
        ],
        timestep,
        combine_last &&
          i == n_years - 1,
        spinup_fraction,
        spinup,
        state_is_grouped,
        growth_scratch,
        mortality_lookup,
        competition_weight_lookup
      );
    }

    double shade_biomass =
      0.0;

    double establishment =
      0.0;

    for (std::size_t i = 0;
         i < state.age.size();
         ++i) {

      const int species_id =
        state.species[i];

      const double age =
        state.age[i];

      const double biomass =
        state.biomass[i];

      if (species_id < 0 ||
          species_id >= n_species)
        stop("Invalid species after CSR growth");

      if (age > 5.0)
        shade_biomass +=
          biomass;

      if (age >= timestep &&
          age > 1.0)
        establishment +=
          biomass;

      if (biomass > 0.0 &&
          age >= maturity[species_id]) {
        mature(
          s,
          species_id
        ) = true;
      }
    }

    establishment_biomass[s] =
      establishment;

    const double actual =
      std::min(
        eco_max[eco],
        shade_biomass
      );

    const double relative =
      actual /
      eco_max[eco];

    int cls = 0;

    for (int h = 0;
         h < eco_thresholds.ncol();
         ++h) {

      const double threshold =
        eco_thresholds(
          eco,
          h
        );

      if (!NumericVector::is_na(
            threshold
          ) &&
          relative >= threshold) {
        cls = h + 1;
      }
    }

    if (cls < 0 ||
        cls >= n_shade_classes)
      stop("Invalid shade class in CSR growth");
    shade_class[s] = cls;

    const std::vector<double>& probability_row =
      base_probability_lookup[
        static_cast<std::size_t>(
          eco *
          n_shade_classes +
          cls
        )
      ];

    for (int j = 0;
         j < n_species;
         ++j) {

      base_probability(
        s,
        j
      ) = probability_row[
        static_cast<std::size_t>(j)
      ];
    }

    for (std::size_t i = 0;
         i < state.age.size();
         ++i) {

      out_species.push_back(
        state.species[i]
      );

      out_age.push_back(
        state.age[i]
      );

      out_biomass.push_back(
        state.biomass[i]
      );
    }

    if (out_species.size() >
        static_cast<std::size_t>(
          std::numeric_limits<int>::max()
        ))
      stop("CSR cohort count exceeds integer offset capacity");

    out_offset[s + 1] =
      static_cast<int>(
        out_species.size()
      );

    bs52_progress_update(
      progress_label,
      s + 1,
      n_states,
      progress_percent,
      show_progress
    );
  }

  IntegerVector species_out =
    wrap(
      out_species
    );

  NumericVector age_out =
    wrap(
      out_age
    );

  NumericVector biomass_out =
    wrap(
      out_biomass
    );

  List out_store =
    make_csr_store(
      out_offset,
      species_out,
      age_out,
      biomass_out
    );

  return List::create(
    _["states"] = out_store,
    _["mature"] = mature,
    _["establishment_biomass"] =
      establishment_biomass,
    _["base_probability"] =
      base_probability,
    _["shade_class"] = shade_class,
    _["weight_lookup_fallbacks"] =
      static_cast<double>(
        competition_weight_lookup.fallback_count
      )
  );
}


// [[Rcpp::export]]
List bs52_rebuild_csr_cpp(
    List store,
    IntegerVector state_eco,
    IntegerVector cell_state,
    IntegerVector recruit_mask,
    int n_species,
    NumericMatrix max_anpp_by_eco,
    NumericVector eco_max_by_eco,
    NumericVector establishment_biomass,
    bool show_progress = false,
    std::string progress_label = "Rebuild CSR") {

  validate_csr_store(
    store,
    n_species
  );

  IntegerVector offset =
    store["offset"];

  IntegerVector input_species =
    store["species_id"];

  NumericVector input_age =
    store["age"];

  NumericVector input_biomass =
    store["biomass"];

  const int n_states =
    offset.size() - 1;

  const int n_sites =
    cell_state.size();

  const int n_eco =
    max_anpp_by_eco.nrow();

  if (state_eco.size() != n_states ||
      recruit_mask.size() != n_sites ||
      max_anpp_by_eco.ncol() != n_species ||
      eco_max_by_eco.size() != n_eco ||
      establishment_biomass.size() != n_states)
    stop("Incompatible CSR rebuild inputs");

  if (n_species < 1 ||
      n_species > 30)
    stop("CSR rebuild supports 1 to 30 species");

  IntegerVector next_cell_state(
    n_sites
  );

  std::unordered_map<
    std::uint64_t,
    int
  > key_to_state;

  key_to_state.reserve(
    static_cast<std::size_t>(
      n_sites * 1.05
    )
  );

  std::vector<int>
    unique_old_state;

  std::vector<int>
    unique_mask;

  std::vector<int>
    unique_eco;

  unique_old_state.reserve(
    static_cast<std::size_t>(
      std::min(
        n_sites,
        n_states * 2
      )
    )
  );

  unique_mask.reserve(
    unique_old_state.capacity()
  );

  unique_eco.reserve(
    unique_old_state.capacity()
  );

  int n_recruited_cells =
    0;

  int progress_percent =
    -1;

  bs52_progress_update(
    progress_label,
    0,
    n_sites,
    progress_percent,
    show_progress
  );

  // Phase 1: preserve the exact v2 ordering of unique
  // (old state, recruitment mask) pairs.
  for (int k = 0;
       k < n_sites;
       ++k) {

    if ((k % 20000) == 0)
      checkUserInterrupt();

    const int old_state_1 =
      cell_state[k];

    if (old_state_1 < 1 ||
        old_state_1 > n_states)
      stop("Invalid cell state in CSR rebuild");

    const int old_state =
      old_state_1 - 1;

    const int mask =
      recruit_mask[k];

    if (mask < 0)
      stop("Invalid recruitment mask in CSR rebuild");

    if (mask != 0)
      ++n_recruited_cells;

    const std::uint64_t key =
      (
        static_cast<std::uint64_t>(
          static_cast<std::uint32_t>(
            old_state_1
          )
        ) << 32
      ) |
      static_cast<std::uint32_t>(
        mask
      );

    const auto found =
      key_to_state.find(
        key
      );

    int next_id;

    if (found ==
        key_to_state.end()) {

      next_id =
        static_cast<int>(
          unique_old_state.size()
        );

      key_to_state.emplace(
        key,
        next_id
      );

      const int eco =
        state_eco[old_state] - 1;

      if (eco < 0 ||
          eco >= n_eco ||
          !std::isfinite(
            eco_max_by_eco[eco]
          ) ||
          eco_max_by_eco[eco] <= 0.0)
        stop("Invalid ecoregion in CSR rebuild");

      unique_old_state.push_back(
        old_state
      );

      unique_mask.push_back(
        mask
      );

      unique_eco.push_back(
        eco + 1
      );

    } else {

      next_id =
        found->second;
    }

    next_cell_state[k] =
      next_id + 1;

    bs52_progress_update(
      progress_label,
      k + 1,
      n_sites,
      progress_percent,
      show_progress
    );
  }

  const int n_unique =
    static_cast<int>(
      unique_old_state.size()
    );

  IntegerVector out_eco(
    n_unique
  );

  IntegerVector out_offset(
    n_unique + 1
  );

  out_offset[0] = 0;

  std::size_t total_cohorts =
    0;

  for (int u = 0;
       u < n_unique;
       ++u) {

    const int old_state =
      unique_old_state[
        static_cast<std::size_t>(u)
      ];

    const int mask =
      unique_mask[
        static_cast<std::size_t>(u)
      ];

    const int old_n =
      offset[old_state + 1] -
      offset[old_state];

    int additions = 0;

    for (int j = 0;
         j < n_species;
         ++j) {

      if (
        (
          static_cast<unsigned int>(
            mask
          ) &
          (1u << j)
        ) != 0u
      ) {
        ++additions;
      }
    }

    total_cohorts +=
      static_cast<std::size_t>(
        old_n + additions
      );

    if (total_cohorts >
        static_cast<std::size_t>(
          std::numeric_limits<int>::max()
        ))
      stop("CSR rebuild exceeds integer offset capacity");

    out_offset[u + 1] =
      static_cast<int>(
        total_cohorts
      );

    out_eco[u] =
      unique_eco[
        static_cast<std::size_t>(u)
      ];
  }

  IntegerVector out_species(
    static_cast<int>(
      total_cohorts
    )
  );

  NumericVector out_age(
    static_cast<int>(
      total_cohorts
    )
  );

  NumericVector out_biomass(
    static_cast<int>(
      total_cohorts
    )
  );

  for (int u = 0;
       u < n_unique;
       ++u) {

    if ((u % 20000) == 0)
      checkUserInterrupt();

    const int old_state =
      unique_old_state[
        static_cast<std::size_t>(u)
      ];

    const int mask =
      unique_mask[
        static_cast<std::size_t>(u)
      ];

    const int eco =
      unique_eco[
        static_cast<std::size_t>(u)
      ] - 1;

    int cursor =
      out_offset[u];

    const int begin =
      offset[old_state];

    const int end =
      offset[old_state + 1];

    for (int k = begin;
         k < end;
         ++k) {

      out_species[cursor] =
        input_species[k];

      out_age[cursor] =
        input_age[k];

      out_biomass[cursor] =
        input_biomass[k];

      ++cursor;
    }

    if (mask != 0) {
      const double site_biomass =
        establishment_biomass[
          old_state
        ];

      const double eco_max =
        eco_max_by_eco[
          eco
        ];

      for (int j = 0;
           j < n_species;
           ++j) {

        if (
          (
            static_cast<unsigned int>(
              mask
            ) &
            (1u << j)
          ) == 0u
        ) {
          continue;
        }

        const double anpp =
          max_anpp_by_eco(
            eco,
            j
          );

        const double candidate =
          anpp *
          std::exp(
            -1.6 *
            site_biomass /
            eco_max
          );

        out_species[cursor] =
          j;

        out_age[cursor] =
          1.0;

        out_biomass[cursor] =
          std::trunc(
            std::max(
              2.0,
              std::min(
                anpp,
                candidate
              )
            )
          );

        ++cursor;
      }
    }

    if (cursor !=
        out_offset[u + 1])
      stop("Internal CSR rebuild offset mismatch");
  }

  List out_store =
    make_csr_store(
      out_offset,
      out_species,
      out_age,
      out_biomass
    );

  return List::create(
    _["states"] = out_store,
    _["state_eco"] = out_eco,
    _["cell_state"] = next_cell_state,
    _["n_recruited_cells"] =
      n_recruited_cells
  );
}


// ============================================================================
// large52 v3.1 Ward
//
// Fuses the 17 species-level R -> C++ Ward calls into one C++ call and stores
// mature source cells as packed bits. The accumulation order for each target
// remains neighborhood-order, exactly as in bs52_ward_failure_log_cpp.
// ============================================================================

// [[Rcpp::export]]
NumericMatrix bs52_ward_probability_fused_bitset_cpp(
    LogicalMatrix state_mature,
    IntegerVector cell_state,
    IntegerMatrix targets,
    int nr,
    int nc,
    IntegerVector neighbor_row,
    IntegerVector neighbor_col,
    NumericMatrix log_failure,
    LogicalVector omnidisperse,
    NumericMatrix base_probability,
    bool local_only = false,
    bool ward2019 = true,
    bool show_progress = false,
    std::string progress_label = "Ward bitset") {

  const int n_sites = cell_state.size();
  const int n_states = state_mature.nrow();
  const int n_species = state_mature.ncol();
  const int n_neighbors = neighbor_row.size();

  if (nr < 1 || nc < 1)
    stop("Invalid Ward landscape dimensions");

  if (targets.nrow() != n_sites ||
      targets.ncol() != 2)
    stop("Incompatible Ward targets");

  if (neighbor_col.size() != n_neighbors ||
      log_failure.nrow() != n_neighbors ||
      log_failure.ncol() != n_species)
    stop("Incompatible Ward neighborhood tables");

  if (omnidisperse.size() != n_species)
    stop("Incompatible Ward omnidisperse vector");

  if (base_probability.nrow() != n_states ||
      base_probability.ncol() != n_species)
    stop("Incompatible Ward base-probability matrix");

  NumericMatrix probability(
    n_sites,
    n_species
  );

  // Local-only dispersal does not need a spatial source grid.
  if (local_only) {
    int progress_percent = -1;

    bs52_progress_update(
      progress_label,
      0,
      n_species,
      progress_percent,
      show_progress
    );

    for (int j = 0; j < n_species; ++j) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        if (state < 0 ||
            state >= n_states)
          stop("Invalid cell state in Ward local kernel");

        const double arrival =
          state_mature(
            state,
            j
          ) == TRUE
            ? 1.0
            : 0.0;

        probability(
          k,
          j
        ) =
          arrival *
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );
    }

    return probability;
  }

  const std::size_t n_cells =
    static_cast<std::size_t>(nr) *
    static_cast<std::size_t>(nc);

  const std::size_t words_per_species =
    (n_cells + 63u) / 64u;

  if (words_per_species == 0u)
    stop("Invalid Ward bitset size");

  std::vector<std::uint64_t> mature_bits(
    words_per_species *
      static_cast<std::size_t>(n_species),
    0u
  );

  std::vector<unsigned char> species_has_mature_source(
    static_cast<std::size_t>(n_species),
    0u
  );

  // Full-grid zero-based index for every active target. This is also the
  // source-grid location because targets are the active cells themselves.
  std::vector<std::size_t> target_index(
    static_cast<std::size_t>(n_sites)
  );

  for (int k = 0; k < n_sites; ++k) {
    const int row =
      targets(
        k,
        0
      ) - 1;

    const int col =
      targets(
        k,
        1
      ) - 1;

    if (row < 0 ||
        row >= nr ||
        col < 0 ||
        col >= nc)
      stop("Ward target out of bounds");

    const int state =
      cell_state[k] - 1;

    if (state < 0 ||
        state >= n_states)
      stop("Invalid cell state in Ward bitset kernel");

    target_index[
      static_cast<std::size_t>(k)
    ] =
      static_cast<std::size_t>(row) +
      static_cast<std::size_t>(col) *
      static_cast<std::size_t>(nr);
  }

  // Pack all mature-source layers once. Inactive cells are never inserted,
  // so a bit test replaces mature(r,c) && active(r,c).
  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      if (state_mature(
            state,
            j
          ) != TRUE)
        continue;

      const std::size_t index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      mature_bits[
        species_base +
        (index >> 6)
      ] |=
        (
          static_cast<std::uint64_t>(1u) <<
          (index & 63u)
        );

      species_has_mature_source[
        static_cast<std::size_t>(j)
      ] = 1u;
    }
  }

  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_species,
    progress_percent,
    show_progress
  );

  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    // Preserve the old active_failure <- probability > 0 behavior,
    // but precompile the four symmetric offsets once per active Ward radius.
    struct WardPreparedNeighbor52 {
      int h;
      int row_offset[4];
      int col_offset[4];
      long long linear_offset[4];
    };

    std::vector<WardPreparedNeighbor52> active_neighbors;

    active_neighbors.reserve(
      static_cast<std::size_t>(
        n_neighbors
      )
    );

    int max_active_row_offset =
      0;

    int max_active_col_offset =
      0;

    for (int h = 0; h < n_neighbors; ++h) {
      if (log_failure(
            h,
            j
          ) >= 0.0) {
        continue;
      }

      const int dr =
        neighbor_row[h];

      const int dc =
        neighbor_col[h];

      WardPreparedNeighbor52 prepared;

      prepared.h =
        h;

      if (dc == 0) {
        prepared.row_offset[0] = dr;
        prepared.col_offset[0] = 0;

        prepared.row_offset[1] = 0;
        prepared.col_offset[1] = dr;

        prepared.row_offset[2] = -dr;
        prepared.col_offset[2] = 0;

        prepared.row_offset[3] = 0;
        prepared.col_offset[3] = -dr;
      } else {
        prepared.row_offset[0] = dr;
        prepared.col_offset[0] = dc;

        prepared.row_offset[1] = -dr;
        prepared.col_offset[1] = dc;

        prepared.row_offset[2] = -dr;
        prepared.col_offset[2] = -dc;

        prepared.row_offset[3] = dr;
        prepared.col_offset[3] = -dc;
      }

      for (int quadrant = 0;
           quadrant < 4;
           ++quadrant) {

        prepared.linear_offset[
          quadrant
        ] =
          static_cast<long long>(
            prepared.row_offset[
              quadrant
            ]
          ) +
          static_cast<long long>(
            prepared.col_offset[
              quadrant
            ]
          ) *
          static_cast<long long>(
            nr
          );

        max_active_row_offset =
          std::max(
            max_active_row_offset,
            std::abs(
              prepared.row_offset[
                quadrant
              ]
            )
          );

        max_active_col_offset =
          std::max(
            max_active_col_offset,
            std::abs(
              prepared.col_offset[
                quadrant
              ]
            )
          );
      }

      active_neighbors.push_back(
        prepared
      );
    }

    if (omnidisperse[j] == TRUE) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        probability(
          k,
          j
        ) =
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );

      continue;
    }

    if (species_has_mature_source[
          static_cast<std::size_t>(j)
        ] == 0u) {

      const double zero_arrival =
        -std::expm1(
          0.0
        );

      for (int k = 0;
           k < n_sites;
           ++k) {

        const int state =
          cell_state[k] - 1;

        probability(
          k,
          j
        ) =
          zero_arrival *
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );

      continue;
    }

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      const int target_row =
        targets(
          k,
          0
        ) - 1;

      const int target_col =
        targets(
          k,
          1
        ) - 1;

      const std::size_t local_index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      const bool local_mature =
        (
          mature_bits[
            species_base +
            (local_index >> 6)
          ] &
          (
            static_cast<std::uint64_t>(1u) <<
            (local_index & 63u)
          )
        ) != 0u;

      double arrival;

      if (local_mature) {
        arrival = 1.0;
      } else {
        double failure = 0.0;

        const bool interior_target =
          target_row >=
            max_active_row_offset &&
          target_row <
            nr -
            max_active_row_offset &&
          target_col >=
            max_active_col_offset &&
          target_col <
            nc -
            max_active_col_offset;

        // Each target receives additions in exactly the same neighborhood
        // order as the validated Ward kernel.
        if (interior_target) {
          const long long target_linear =
            static_cast<long long>(
              local_index
            );

          for (std::size_t hh = 0;
               hh < active_neighbors.size();
               ++hh) {

            const WardPreparedNeighbor52& prepared =
              active_neighbors[hh];

            int count =
              0;

            for (int quadrant = 0;
                 quadrant < 4;
                 ++quadrant) {

              const std::size_t source_index =
                static_cast<std::size_t>(
                  target_linear +
                  prepared.linear_offset[
                    quadrant
                  ]
                );

              const bool source_mature =
                (
                  mature_bits[
                    species_base +
                    (source_index >> 6)
                  ] &
                  (
                    static_cast<std::uint64_t>(1u) <<
                    (source_index & 63u)
                  )
                ) != 0u;

              if (source_mature)
                ++count;
            }

            const int power =
              ward2019
                ? count
                : (
                    count > 0
                      ? 1
                      : 0
                  );

            if (power > 0) {
              failure +=
                static_cast<double>(
                  power
                ) *
                log_failure(
                  prepared.h,
                  j
                );
            }
          }
        } else {
          // Rare edge-cell fallback: preserve the validated bounds logic.
          for (std::size_t hh = 0;
               hh < active_neighbors.size();
               ++hh) {

            const WardPreparedNeighbor52& prepared =
              active_neighbors[hh];

            int count =
              0;

            for (int quadrant = 0;
                 quadrant < 4;
                 ++quadrant) {

              const int row =
                target_row +
                prepared.row_offset[
                  quadrant
                ];

              const int col =
                target_col +
                prepared.col_offset[
                  quadrant
                ];

              if (row < 0 ||
                  row >= nr ||
                  col < 0 ||
                  col >= nc)
                continue;

              const std::size_t source_index =
                static_cast<std::size_t>(row) +
                static_cast<std::size_t>(col) *
                static_cast<std::size_t>(nr);

              const bool source_mature =
                (
                  mature_bits[
                    species_base +
                    (source_index >> 6)
                  ] &
                  (
                    static_cast<std::uint64_t>(1u) <<
                    (source_index & 63u)
                  )
                ) != 0u;

              if (source_mature)
                ++count;
            }

            const int power =
              ward2019
                ? count
                : (
                    count > 0
                      ? 1
                      : 0
                  );

            if (power > 0) {
              failure +=
                static_cast<double>(
                  power
                ) *
                log_failure(
                  prepared.h,
                  j
                );
            }
          }
        }

        arrival =
          -std::expm1(
            failure
          );
      }

      probability(
        k,
        j
      ) =
        arrival *
        base_probability(
          state,
          j
        );
    }

    bs52_progress_update(
      progress_label,
      j + 1,
      n_species,
      progress_percent,
      show_progress
    );
  }

  return probability;
}



// [[Rcpp::export]]
NumericMatrix bs52_ward_probability_fused_bitset_v38_reference_cpp(
    LogicalMatrix state_mature,
    IntegerVector cell_state,
    IntegerMatrix targets,
    int nr,
    int nc,
    IntegerVector neighbor_row,
    IntegerVector neighbor_col,
    NumericMatrix log_failure,
    LogicalVector omnidisperse,
    NumericMatrix base_probability,
    bool local_only = false,
    bool ward2019 = true,
    bool show_progress = false,
    std::string progress_label = "Ward bitset") {

  const int n_sites = cell_state.size();
  const int n_states = state_mature.nrow();
  const int n_species = state_mature.ncol();
  const int n_neighbors = neighbor_row.size();

  if (nr < 1 || nc < 1)
    stop("Invalid Ward landscape dimensions");

  if (targets.nrow() != n_sites ||
      targets.ncol() != 2)
    stop("Incompatible Ward targets");

  if (neighbor_col.size() != n_neighbors ||
      log_failure.nrow() != n_neighbors ||
      log_failure.ncol() != n_species)
    stop("Incompatible Ward neighborhood tables");

  if (omnidisperse.size() != n_species)
    stop("Incompatible Ward omnidisperse vector");

  if (base_probability.nrow() != n_states ||
      base_probability.ncol() != n_species)
    stop("Incompatible Ward base-probability matrix");

  NumericMatrix probability(
    n_sites,
    n_species
  );

  // Local-only dispersal does not need a spatial source grid.
  if (local_only) {
    int progress_percent = -1;

    bs52_progress_update(
      progress_label,
      0,
      n_species,
      progress_percent,
      show_progress
    );

    for (int j = 0; j < n_species; ++j) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        if (state < 0 ||
            state >= n_states)
          stop("Invalid cell state in Ward local kernel");

        const double arrival =
          state_mature(
            state,
            j
          ) == TRUE
            ? 1.0
            : 0.0;

        probability(
          k,
          j
        ) =
          arrival *
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );
    }

    return probability;
  }

  const std::size_t n_cells =
    static_cast<std::size_t>(nr) *
    static_cast<std::size_t>(nc);

  const std::size_t words_per_species =
    (n_cells + 63u) / 64u;

  if (words_per_species == 0u)
    stop("Invalid Ward bitset size");

  std::vector<std::uint64_t> mature_bits(
    words_per_species *
      static_cast<std::size_t>(n_species),
    0u
  );

  // Full-grid zero-based index for every active target. This is also the
  // source-grid location because targets are the active cells themselves.
  std::vector<std::size_t> target_index(
    static_cast<std::size_t>(n_sites)
  );

  for (int k = 0; k < n_sites; ++k) {
    const int row =
      targets(
        k,
        0
      ) - 1;

    const int col =
      targets(
        k,
        1
      ) - 1;

    if (row < 0 ||
        row >= nr ||
        col < 0 ||
        col >= nc)
      stop("Ward target out of bounds");

    const int state =
      cell_state[k] - 1;

    if (state < 0 ||
        state >= n_states)
      stop("Invalid cell state in Ward bitset kernel");

    target_index[
      static_cast<std::size_t>(k)
    ] =
      static_cast<std::size_t>(row) +
      static_cast<std::size_t>(col) *
      static_cast<std::size_t>(nr);
  }

  // Pack all mature-source layers once. Inactive cells are never inserted,
  // so a bit test replaces mature(r,c) && active(r,c).
  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      if (state_mature(
            state,
            j
          ) != TRUE)
        continue;

      const std::size_t index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      mature_bits[
        species_base +
        (index >> 6)
      ] |=
        (
          static_cast<std::uint64_t>(1u) <<
          (index & 63u)
        );
    }
  }

  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_species,
    progress_percent,
    show_progress
  );

  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    // Preserve the old active_failure <- probability > 0 behavior.
    std::vector<int> active_neighbors;
    active_neighbors.reserve(
      static_cast<std::size_t>(
        n_neighbors
      )
    );

    for (int h = 0; h < n_neighbors; ++h) {
      if (log_failure(
            h,
            j
          ) < 0.0) {
        active_neighbors.push_back(
          h
        );
      }
    }

    if (omnidisperse[j] == TRUE) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        probability(
          k,
          j
        ) =
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );

      continue;
    }

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      const int target_row =
        targets(
          k,
          0
        ) - 1;

      const int target_col =
        targets(
          k,
          1
        ) - 1;

      const std::size_t local_index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      const bool local_mature =
        (
          mature_bits[
            species_base +
            (local_index >> 6)
          ] &
          (
            static_cast<std::uint64_t>(1u) <<
            (local_index & 63u)
          )
        ) != 0u;

      double arrival;

      if (local_mature) {
        arrival = 1.0;
      } else {
        double failure = 0.0;

        // Each target receives additions in exactly the same neighborhood
        // order as the validated Ward kernel.
        for (std::size_t hh = 0;
             hh < active_neighbors.size();
             ++hh) {

          const int h =
            active_neighbors[hh];

          const int dr =
            neighbor_row[h];

          const int dc =
            neighbor_col[h];

          int row_offset[4];
          int col_offset[4];

          if (dc == 0) {
            row_offset[0] = dr;
            col_offset[0] = 0;

            row_offset[1] = 0;
            col_offset[1] = dr;

            row_offset[2] = -dr;
            col_offset[2] = 0;

            row_offset[3] = 0;
            col_offset[3] = -dr;
          } else {
            row_offset[0] = dr;
            col_offset[0] = dc;

            row_offset[1] = -dr;
            col_offset[1] = dc;

            row_offset[2] = -dr;
            col_offset[2] = -dc;

            row_offset[3] = dr;
            col_offset[3] = -dc;
          }

          int count = 0;

          for (int quadrant = 0;
               quadrant < 4;
               ++quadrant) {

            const int row =
              target_row +
              row_offset[
                quadrant
              ];

            const int col =
              target_col +
              col_offset[
                quadrant
              ];

            if (row < 0 ||
                row >= nr ||
                col < 0 ||
                col >= nc)
              continue;

            const std::size_t source_index =
              static_cast<std::size_t>(row) +
              static_cast<std::size_t>(col) *
              static_cast<std::size_t>(nr);

            const bool source_mature =
              (
                mature_bits[
                  species_base +
                  (source_index >> 6)
                ] &
                (
                  static_cast<std::uint64_t>(1u) <<
                  (source_index & 63u)
                )
              ) != 0u;

            if (source_mature)
              ++count;
          }

          const int power =
            ward2019
              ? count
              : (
                  count > 0
                    ? 1
                    : 0
                );

          if (power > 0) {
            failure +=
              static_cast<double>(
                power
              ) *
              log_failure(
                h,
                j
              );
          }
        }

        arrival =
          -std::expm1(
            failure
          );
      }

      probability(
        k,
        j
      ) =
        arrival *
        base_probability(
          state,
          j
        );
    }

    bs52_progress_update(
      progress_label,
      j + 1,
      n_species,
      progress_percent,
      show_progress
    );
  }

  return probability;
}


// [[Rcpp::export]]
NumericMatrix bs52_ward_probability_fused_bitset_v310_reference_cpp(
    LogicalMatrix state_mature,
    IntegerVector cell_state,
    IntegerMatrix targets,
    int nr,
    int nc,
    IntegerVector neighbor_row,
    IntegerVector neighbor_col,
    NumericMatrix log_failure,
    LogicalVector omnidisperse,
    NumericMatrix base_probability,
    bool local_only = false,
    bool ward2019 = true,
    bool show_progress = false,
    std::string progress_label = "Ward bitset") {

  const int n_sites = cell_state.size();
  const int n_states = state_mature.nrow();
  const int n_species = state_mature.ncol();
  const int n_neighbors = neighbor_row.size();

  if (nr < 1 || nc < 1)
    stop("Invalid Ward landscape dimensions");

  if (targets.nrow() != n_sites ||
      targets.ncol() != 2)
    stop("Incompatible Ward targets");

  if (neighbor_col.size() != n_neighbors ||
      log_failure.nrow() != n_neighbors ||
      log_failure.ncol() != n_species)
    stop("Incompatible Ward neighborhood tables");

  if (omnidisperse.size() != n_species)
    stop("Incompatible Ward omnidisperse vector");

  if (base_probability.nrow() != n_states ||
      base_probability.ncol() != n_species)
    stop("Incompatible Ward base-probability matrix");

  NumericMatrix probability(
    n_sites,
    n_species
  );

  // Local-only dispersal does not need a spatial source grid.
  if (local_only) {
    int progress_percent = -1;

    bs52_progress_update(
      progress_label,
      0,
      n_species,
      progress_percent,
      show_progress
    );

    for (int j = 0; j < n_species; ++j) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        if (state < 0 ||
            state >= n_states)
          stop("Invalid cell state in Ward local kernel");

        const double arrival =
          state_mature(
            state,
            j
          ) == TRUE
            ? 1.0
            : 0.0;

        probability(
          k,
          j
        ) =
          arrival *
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );
    }

    return probability;
  }

  const std::size_t n_cells =
    static_cast<std::size_t>(nr) *
    static_cast<std::size_t>(nc);

  const std::size_t words_per_species =
    (n_cells + 63u) / 64u;

  if (words_per_species == 0u)
    stop("Invalid Ward bitset size");

  std::vector<std::uint64_t> mature_bits(
    words_per_species *
      static_cast<std::size_t>(n_species),
    0u
  );

  std::vector<unsigned char> species_has_mature_source(
    static_cast<std::size_t>(n_species),
    0u
  );

  // Full-grid zero-based index for every active target. This is also the
  // source-grid location because targets are the active cells themselves.
  std::vector<std::size_t> target_index(
    static_cast<std::size_t>(n_sites)
  );

  for (int k = 0; k < n_sites; ++k) {
    const int row =
      targets(
        k,
        0
      ) - 1;

    const int col =
      targets(
        k,
        1
      ) - 1;

    if (row < 0 ||
        row >= nr ||
        col < 0 ||
        col >= nc)
      stop("Ward target out of bounds");

    const int state =
      cell_state[k] - 1;

    if (state < 0 ||
        state >= n_states)
      stop("Invalid cell state in Ward bitset kernel");

    target_index[
      static_cast<std::size_t>(k)
    ] =
      static_cast<std::size_t>(row) +
      static_cast<std::size_t>(col) *
      static_cast<std::size_t>(nr);
  }

  // Pack all mature-source layers once. Inactive cells are never inserted,
  // so a bit test replaces mature(r,c) && active(r,c).
  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      if (state_mature(
            state,
            j
          ) != TRUE)
        continue;

      const std::size_t index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      mature_bits[
        species_base +
        (index >> 6)
      ] |=
        (
          static_cast<std::uint64_t>(1u) <<
          (index & 63u)
        );

      species_has_mature_source[
        static_cast<std::size_t>(j)
      ] = 1u;
    }
  }

  int progress_percent = -1;

  bs52_progress_update(
    progress_label,
    0,
    n_species,
    progress_percent,
    show_progress
  );

  for (int j = 0; j < n_species; ++j) {
    const std::size_t species_base =
      static_cast<std::size_t>(j) *
      words_per_species;

    // Preserve the old active_failure <- probability > 0 behavior.
    std::vector<int> active_neighbors;
    active_neighbors.reserve(
      static_cast<std::size_t>(
        n_neighbors
      )
    );

    for (int h = 0; h < n_neighbors; ++h) {
      if (log_failure(
            h,
            j
          ) < 0.0) {
        active_neighbors.push_back(
          h
        );
      }
    }

    if (omnidisperse[j] == TRUE) {
      for (int k = 0; k < n_sites; ++k) {
        const int state =
          cell_state[k] - 1;

        probability(
          k,
          j
        ) =
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );

      continue;
    }

    if (species_has_mature_source[
          static_cast<std::size_t>(j)
        ] == 0u) {

      const double zero_arrival =
        -std::expm1(
          0.0
        );

      for (int k = 0;
           k < n_sites;
           ++k) {

        const int state =
          cell_state[k] - 1;

        probability(
          k,
          j
        ) =
          zero_arrival *
          base_probability(
            state,
            j
          );
      }

      bs52_progress_update(
        progress_label,
        j + 1,
        n_species,
        progress_percent,
        show_progress
      );

      continue;
    }

    for (int k = 0; k < n_sites; ++k) {
      const int state =
        cell_state[k] - 1;

      const int target_row =
        targets(
          k,
          0
        ) - 1;

      const int target_col =
        targets(
          k,
          1
        ) - 1;

      const std::size_t local_index =
        target_index[
          static_cast<std::size_t>(k)
        ];

      const bool local_mature =
        (
          mature_bits[
            species_base +
            (local_index >> 6)
          ] &
          (
            static_cast<std::uint64_t>(1u) <<
            (local_index & 63u)
          )
        ) != 0u;

      double arrival;

      if (local_mature) {
        arrival = 1.0;
      } else {
        double failure = 0.0;

        // Each target receives additions in exactly the same neighborhood
        // order as the validated Ward kernel.
        for (std::size_t hh = 0;
             hh < active_neighbors.size();
             ++hh) {

          const int h =
            active_neighbors[hh];

          const int dr =
            neighbor_row[h];

          const int dc =
            neighbor_col[h];

          int row_offset[4];
          int col_offset[4];

          if (dc == 0) {
            row_offset[0] = dr;
            col_offset[0] = 0;

            row_offset[1] = 0;
            col_offset[1] = dr;

            row_offset[2] = -dr;
            col_offset[2] = 0;

            row_offset[3] = 0;
            col_offset[3] = -dr;
          } else {
            row_offset[0] = dr;
            col_offset[0] = dc;

            row_offset[1] = -dr;
            col_offset[1] = dc;

            row_offset[2] = -dr;
            col_offset[2] = -dc;

            row_offset[3] = dr;
            col_offset[3] = -dc;
          }

          int count = 0;

          for (int quadrant = 0;
               quadrant < 4;
               ++quadrant) {

            const int row =
              target_row +
              row_offset[
                quadrant
              ];

            const int col =
              target_col +
              col_offset[
                quadrant
              ];

            if (row < 0 ||
                row >= nr ||
                col < 0 ||
                col >= nc)
              continue;

            const std::size_t source_index =
              static_cast<std::size_t>(row) +
              static_cast<std::size_t>(col) *
              static_cast<std::size_t>(nr);

            const bool source_mature =
              (
                mature_bits[
                  species_base +
                  (source_index >> 6)
                ] &
                (
                  static_cast<std::uint64_t>(1u) <<
                  (source_index & 63u)
                )
              ) != 0u;

            if (source_mature)
              ++count;
          }

          const int power =
            ward2019
              ? count
              : (
                  count > 0
                    ? 1
                    : 0
                );

          if (power > 0) {
            failure +=
              static_cast<double>(
                power
              ) *
              log_failure(
                h,
                j
              );
          }
        }

        arrival =
          -std::expm1(
            failure
          );
      }

      probability(
        k,
        j
      ) =
        arrival *
        base_probability(
          state,
          j
        );
    }

    bs52_progress_update(
      progress_label,
      j + 1,
      n_species,
      progress_percent,
      show_progress
    );
  }

  return probability;
}
