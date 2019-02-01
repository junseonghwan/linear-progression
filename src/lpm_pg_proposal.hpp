//
//  lpm_pg_proposal.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2019-01-27.
//

#ifndef lpm_pg_proposal_hpp
#define lpm_pg_proposal_hpp

#include <vector>

#include <spf/pg_proposal.hpp>

#include "lpm_state.hpp"
#include "lpm_params.hpp"

using namespace std;

class LPMParamProposal : public PGProposal<LinearProgressionState, LinearProgressionParameters>
{
    size_t n_patients;
    size_t n_pathways;
    size_t n_genes;
    size_t n_mh_iter;
    gsl_matrix_view &obs_view;
    vector<size_t> &row_sum;
    vector<size_t> *stages = 0;
    double a;
    double b;
    double mh_proposal_sd = 0.02;

public:
    LPMParamProposal(gsl_matrix_view &obs_view, vector<size_t> &row_sum, size_t n_patients, size_t n_pathways, size_t n_genes, size_t n_mh_iters, double a = 0, double b = 0.3);
    LinearProgressionParameters *sample_from_prior(gsl_rng *random);
    LinearProgressionParameters *propose(gsl_rng *random, LinearProgressionParameters *curr, ParticleGenealogy<LinearProgressionState> *genealogy);
    double log_prior(LinearProgressionParameters *curr);
    inline vector<size_t> &get_stages() { return *stages; }
};

#endif /* lpm_pg_proposal_hpp */
