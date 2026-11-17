// 252: detect mature pre-fire serotinous sources in compact CSR state.
// Strictly deterministic: no R RNG use and no change to simulator objects.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
Rcpp::List basefire_prefire_csr_252(
    Rcpp::List state_store,
    Rcpp::IntegerVector cell_state,
    Rcpp::IntegerVector active_positions,
    Rcpp::IntegerVector cell_ids,
    Rcpp::IntegerVector serotiny_species_1based,
    Rcpp::NumericVector maturity,
    int species_id_base) {
  if (active_positions.size()!=cell_ids.size())
    Rcpp::stop("252: incoherent cells/positions");
  if (serotiny_species_1based.size()!=maturity.size())
    Rcpp::stop("252: incoherent species/maturity");
  if (species_id_base != 0 && species_id_base != 1)
    Rcpp::stop("252: unknown species_id_base");
  Rcpp::IntegerVector off = Rcpp::as<Rcpp::IntegerVector>(state_store["offset"]);
  Rcpp::IntegerVector species = Rcpp::as<Rcpp::IntegerVector>(state_store["species_id"]);
  Rcpp::NumericVector ages = Rcpp::as<Rcpp::NumericVector>(state_store["age"]);
  Rcpp::NumericVector biom = Rcpp::as<Rcpp::NumericVector>(state_store["biomass"]);
  if(off.size()<2 || species.size()!=ages.size() || species.size()!=biom.size())
    Rcpp::stop("252: malformed CSR state store");
  // Accept either zero- or one-based offsets, detected from the store itself.
  int offset_base = off[0];
  if((offset_base!=0 && offset_base!=1) ||
     off[off.size()-1]-offset_base != species.size())
    Rcpp::stop("252: unsupported CSR offset convention");
  for (int j=1; j<off.size(); ++j) {
    if (off[j]<off[j-1]) Rcpp::stop("252: non-monotone offsets");
  }
  const int n_states = off.size()-1;
  const int n_sero = serotiny_species_1based.size();
  const int n_sites = active_positions.size();
  Rcpp::List results(n_sero);
  for (int k=0; k<n_sero; ++k) {
    const int desired = serotiny_species_1based[k]-1+species_id_base;
    const double min_age = maturity[k];
    std::vector<int> out;
    out.reserve(static_cast<size_t>(n_sites/4+1));
    for (int i=0; i<n_sites; ++i) {
      const int p=active_positions[i]-1;
      if(p<0 || p>=cell_state.size()) Rcpp::stop("252: invalid active position");
      const int s=cell_state[p]-1;
      if(s<0 || s>=n_states) Rcpp::stop("252: invalid cell_state index");
      const int from=off[s]-offset_base;
      const int to=off[s+1]-offset_base;
      bool mature=false;
      for (int j=from; j<to; ++j) {
        if (species[j]==desired && ages[j]>=min_age && biom[j]>0) {
          mature=true;
          break;
        }
      }
      if(mature) out.push_back(cell_ids[i]);
    }
    results[k]=Rcpp::wrap(out);
  }
  return results;
}
