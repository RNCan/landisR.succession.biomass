// 256: fused per-site post-fire serotiny cohort reconstruction in C++.
// Mirrors LANDIS-II's access-by-site cohort update pattern without copying
// the entire forest into a list of per-cell R data.frames.
// No RNG calls, all probability calculations and RNG remain in the original R.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

struct Cohort256 {
  int sp;
  double age;
  double bio;
};
struct Record256 {
  int cell;
  int sp;
  double arrival;
  double probability;
  double biomass;
};

// [[Rcpp::export(rng = false)]]
Rcpp::List basefire_serotiny_build_256(
    Rcpp::List csr, Rcpp::IntegerVector cell_state,
    Rcpp::IntegerVector active_ids, Rcpp::IntegerVector active_eco,
    Rcpp::IntegerMatrix recruited_pairs, Rcpp::IntegerVector active_sero,
    Rcpp::CharacterVector species_names, Rcpp::NumericMatrix arrival,
    Rcpp::NumericMatrix probability, Rcpp::IntegerVector eco_keys,
    Rcpp::NumericMatrix anpp_by_eco, Rcpp::NumericVector ecomax_by_eco,
    int timestep) {
  const IntegerVector offsets=csr["offset"], sp_codes=csr["species_id"];
  const NumericVector ages=csr["age"], biomass=csr["biomass"];
  int nsite=active_ids.size(), nstates=offsets.size()-1;
  int nsp=species_names.size();
  if(nsite!=cell_state.size() || nsite!=active_eco.size() || nstates<1 ||
     ages.size()!=sp_codes.size() || biomass.size()!=sp_codes.size() ||
     offsets[0]!=0 || offsets[nstates]!=sp_codes.size() ||
     arrival.nrow()!=nsite || probability.nrow()!=nsite ||
     arrival.ncol()!=active_sero.size() ||
     probability.ncol()!=active_sero.size() ||
     anpp_by_eco.nrow()!=eco_keys.size() ||
     anpp_by_eco.ncol()!=nsp ||
     ecomax_by_eco.size()!=eco_keys.size() ||
     recruited_pairs.ncol()!=2) stop("256: dimensions CSR/recrutement incoherentes");
  for(int s=0;s<nstates;++s)
    if(offsets[s]<0 || offsets[s+1]<offsets[s] || offsets[s+1]>sp_codes.size())
      stop("256: offsets CSR non monotones");
  std::map<int,int> eco_index;
  for(int k=0;k<eco_keys.size();++k)
    eco_index[eco_keys[k]]=k;
  // R's split(pairs[,2],pairs[,1]) sorts groups by active position;
  // within each group columns retain their original pair order.
  std::map<int,std::vector<int> > sites;
  for(int k=0;k<recruited_pairs.nrow();++k) {
    int pos=recruited_pairs(k,0)-1;
    int col=recruited_pairs(k,1)-1;
    if(pos<0 || pos>=nsite || col<0 || col>=active_sero.size())
      stop("256: paire active position / espèce invalide");
    sites[pos].push_back(col);
  }
  std::vector<int> out_cell, out_sp;
  std::vector<double> out_age, out_bio;
  std::vector<Record256> records;
  records.reserve(recruited_pairs.nrow());
  for(const auto& site: sites) {
    int pos=site.first, cell=active_ids[pos], state=cell_state[pos]-1;
    if(state<0 || state>=nstates) stop("256: cellule référence un état inexistant");
    auto it=eco_index.find(active_eco[pos]);
    if(it==eco_index.end()) stop("256: écorégion absente des paramètres dynamiques");
    int e=it->second;
    int a=offsets[state], b=offsets[state+1];
    std::vector<Cohort256> v;
    v.reserve(static_cast<size_t>(b-a+site.second.size()));
    // Original pre-recruitment biomass, excluding seedlings and young cohorts.
    double site_biomass=0;
    for(int j=a;j<b;++j) {
      int sp=sp_codes[j]; // landisR beta.3 CSR species coding is 0-based.
      if(sp<0 || sp>=nsp) stop("256: code d'espèce CSR invalide");
      v.push_back({sp,ages[j],biomass[j]});
      if(ages[j]>=timestep && ages[j]>1) site_biomass+=biomass[j];
    }
    const double eco_max=ecomax_by_eco[e];
    for(int col:site.second) {
      int sp=active_sero[col]-1; // R species indexes are 1-based.
      if(sp<0 || sp>=nsp) stop("256: espèce recrutée absente");
      double max_anpp=anpp_by_eco(e,sp);
      if(!std::isfinite(max_anpp) || !std::isfinite(eco_max) || eco_max<=0)
        stop("256: paramètre ANPP / biomasse max non valide");
      // Exactly the same arithmetic grouping and truncation as the R source.
      double bi=std::trunc(std::max(2.0,
                  std::min(max_anpp,max_anpp * std::exp(-1.6*site_biomass/eco_max))));
      bool updated=false;
      for(auto& cohort:v) {
        if(cohort.sp==sp && cohort.age==1.0) {
          cohort.bio+=bi;
          updated=true;
          break;
        }
      }
      if(!updated) v.push_back({sp,1.0,bi});
      records.push_back({cell,sp,arrival(pos,col),probability(pos,col),bi});
    }
    std::stable_sort(v.begin(),v.end(),[](const Cohort256& x,const Cohort256& y){
      if(x.sp!=y.sp) return x.sp<y.sp;
      return x.age<y.age;
    });
    for(const auto& c:v) {
      out_cell.push_back(cell);out_sp.push_back(c.sp);
      out_age.push_back(c.age);out_bio.push_back(c.bio);
    }
  }
  const int n=out_cell.size(), nr=records.size();
  IntegerVector ids(n),rids(nr);
  CharacterVector spp(n),rspp(nr);
  NumericVector age(n),bio(n),rbio(nr),rarri(nr),rprob(nr);
  for(int i=0;i<n;++i) {
    ids[i]=out_cell[i]; age[i]=out_age[i];bio[i]=out_bio[i];
    SET_STRING_ELT(spp,i,STRING_ELT(species_names,out_sp[i]));
  }
  for(int i=0;i<nr;++i) {
    rids[i]=records[i].cell;rarri[i]=records[i].arrival;
    rprob[i]=records[i].probability;rbio[i]=records[i].biomass;
    SET_STRING_ELT(rspp,i,STRING_ELT(species_names,records[i].sp));
  }
  // data.frame column names/order match exactly R's do.call(rbind, ...).
  Rcpp::DataFrame replacement=Rcpp::DataFrame::create(
    _["cell_id"]=ids,_["species"]=spp,_["age"]=age,_["biomass"]=bio);
  Rcpp::DataFrame recruitment=Rcpp::DataFrame::create(
    _["cell_id"]=rids,_["species"]=rspp,
    _["seed_probability"]=rarri,_["establishment_probability"]=rprob,
    _["biomass"]=rbio);
  return Rcpp::List::create(_["replacement"]=replacement,
                            _["recruitment"]=recruitment,
                            _["cells"]=static_cast<int>(sites.size()),
                            _["cohorts"]=nr);
}
