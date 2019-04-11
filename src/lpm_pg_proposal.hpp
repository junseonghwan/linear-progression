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
    size_t n_mh_iter;
    double fbp_max;
    double bgp_max;
    double mh_proposal_sd = 0.1;
//    const gsl_matrix &obs_matrix;
//    const vector<size_t> &row_sums;

    void sample_separately(gsl_rng *random, const LinearProgressionState &state, LinearProgressionParameters &new_param);
    void sample_together(gsl_rng *random, const LinearProgressionState &state, LinearProgressionParameters &new_param);

    // store the parameters
    vector<double> bgps;
    vector<double> fbps;

public:
//    LPMParamProposal(const gsl_matrix &obs_matrix, const vector<size_t> &row_sums, size_t n_mh_iters, double fbp_max, double bgp_max);
    LPMParamProposal(size_t n_mh_iters, double fbp_max, double bgp_max);
    shared_ptr<LinearProgressionParameters> sample_from_prior(gsl_rng *random);
    shared_ptr<LinearProgressionParameters> propose(gsl_rng *random, const LinearProgressionParameters &curr, shared_ptr<ParticleGenealogy<LinearProgressionState>> genealogy);
    double log_prior(const LinearProgressionParameters &curr);
    inline vector<double> &get_bgps() { return bgps; }
    inline vector<double> &get_fbps() { return fbps; }
};

#endif /* lpm_pg_proposal_hpp */
