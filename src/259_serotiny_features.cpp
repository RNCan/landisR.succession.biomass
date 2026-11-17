// 259 | Serotiny establishment base probabilities directly from LANDIS-R CSR.
// Mirrors bs52_large_state_features52(..., regular_recruitment=FALSE)
// for the ONLY result consumed by post-fire serotiny: $base_probability.
// All dynamic parameter lookups and probability arithmetic are performed
// by the original R helpers, once per ecoregion/shade class.
// No RNG, no cohort mutation, no C++ persistence outside this call.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix basefire_serotiny_baseprob_csr_259(
    Rcpp::List csr,
    Rcpp::IntegerVector state_eco,
    Rcpp::IntegerVector eco_keys,
    Rcpp::List shade_thresholds,
    Rcpp::NumericVector max_biomass_by_eco,
    Rcpp::List class_probability_by_eco) {
  Rcpp::IntegerVector offset=csr["offset"];
  Rcpp::IntegerVector species=csr["species_id"];
  Rcpp::NumericVector age=csr["age"];
  Rcpp::NumericVector biomass=csr["biomass"];
  const int ns=state_eco.size();
  const int ne=eco_keys.size();
  if (ns<1 || ne<1 || offset.size()!=ns+1 ||
      age.size()!=species.size() || biomass.size()!=species.size() ||
      shade_thresholds.size()!=ne || max_biomass_by_eco.size()!=ne ||
      class_probability_by_eco.size()!=ne || offset[0]!=0 ||
      offset[ns]!=age.size()) stop("259: dimensions CSR ou ecoregions incompatibles");
  std::unordered_map<int,int> eco_lookup;
  eco_lookup.reserve(ne*2);
  std::vector<Rcpp::NumericVector> thresholds;
  std::vector<Rcpp::NumericMatrix> prob;
  thresholds.reserve(ne); prob.reserve(ne);
  for (int e=0;e<ne;++e) {
    if (eco_keys[e]==NA_INTEGER || eco_lookup.count(eco_keys[e]))
      stop("259: code ecoregion manquant ou duplique");
    eco_lookup[eco_keys[e]]=e;
    thresholds.push_back(Rcpp::as<Rcpp::NumericVector>(shade_thresholds[e]));
    prob.push_back(Rcpp::as<Rcpp::NumericMatrix>(class_probability_by_eco[e]));
    if (prob.back().nrow()!=thresholds.back().size()+1 ||
        !(max_biomass_by_eco[e]>0) || !R_finite(max_biomass_by_eco[e]))
      stop("259: table des classes d ombre invalide");
  }
  const int nsp=prob[0].ncol();
  for (int e=1;e<ne;++e)
    if (prob[e].ncol()!=nsp) stop("259: nombre d especes non constant");
  Rcpp::NumericMatrix out(ns,nsp);
  for (int s=0;s<ns;++s) {
    const auto it=eco_lookup.find(state_eco[s]);
    if (it==eco_lookup.end()) stop("259: ecoregion d etat introuvable");
    const int e=it->second;
    const int from=offset[s],to=offset[s+1];
    if (from<0 || to<from || to>age.size()) stop("259: CSR offset invalide");
    // Match the R expression sum(biomass[age > 5]), in CSR row order.
    double shade_biomass=0.0;
    for (int c=from;c<to;++c) {
      if (age[c]>5.0) shade_biomass += biomass[c];
    }
    const double maxeco=max_biomass_by_eco[e];
    const double relative=std::min(maxeco,shade_biomass)/maxeco;
    int shade_class=0;
    const Rcpp::NumericVector &th=thresholds[e];
    for (int k=0;k<th.size();++k) {
      if (relative>=th[k]) shade_class=k+1;
    }
    const Rcpp::NumericMatrix &tab=prob[e];
    for (int j=0;j<nsp;++j) out(s,j)=tab(shade_class,j);
  }
  return out;
}
