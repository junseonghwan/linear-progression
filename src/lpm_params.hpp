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
    double fbp;
    double bgp;
public:
    LinearProgressionParameters(double fbp, double bgp);
    double get_fbp() const;
    double get_bgp() const;

    void set_fbp(double fbp);
    void set_bgp(double bgp);
};

#endif /* lpm_params_hpp */
