//
//  lpm_params.hpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#ifndef lpm_params_hpp
#define lpm_params_hpp

#include <vector>

using namespace std;

class LinearProgressionParameters
{
    vector<size_t> *stages = 0;
    double fbp;
    double bgp;
public:
    LinearProgressionParameters(double fbp, double bgp);
    LinearProgressionParameters(double fbp, double bgp, vector<size_t> &stages);
    bool has_patient_stages();
    vector<size_t> &get_patient_progression_stages();
    void set_patient_progression_stages(vector<size_t> &stages);
    double get_fbp();
    double get_bgp();
    void set_fbp(double fbp);
    void set_bgp(double bgp);
};

#endif /* lpm_params_hpp */
