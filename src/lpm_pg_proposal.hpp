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
    size_t n_mh_iter;
    vector<size_t> *stages = 0;
    double fbp_max;
    double bgp_max;
    double mh_proposal_sd = 0.02;

public:
    LPMParamProposal(size_t n_genes, size_t n_mh_iters, double fbp_max, double bgp_max);
    LinearProgressionParameters *sample_from_prior(gsl_rng *random);
    LinearProgressionParameters *propose(gsl_rng *random, LinearProgressionParameters *curr, ParticleGenealogy<LinearProgressionState> *genealogy);
    double log_prior(LinearProgressionParameters *curr);
    inline vector<size_t> &get_stages() { return *stages; }
};

#endif /* lpm_pg_proposal_hpp */
