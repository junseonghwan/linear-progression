//
//  lpm_params.cpp
//  lpm
//
//  Created by Seong-Hwan Jun on 2018-12-26.
//

#include "lpm_params.hpp"

LinearProgressionParameters::LinearProgressionParameters(double fbp, double bgp) :
fbp(fbp), bgp(bgp)
{
    
}

LinearProgressionParameters::LinearProgressionParameters(double fbp, double bgp, vector<size_t> &stages) :
LinearProgressionParameters(fbp, bgp)
{
    this->stages = &stages;
}

bool LinearProgressionParameters::has_patient_stages()
{
    return (this->stages != 0);
}

vector<size_t> &LinearProgressionParameters::get_patient_progression_stages()
{
    return *stages;
}

void LinearProgressionParameters::set_patient_progression_stages(vector<size_t> &stages)
{
    this->stages = &stages;
}

double LinearProgressionParameters::get_fbp()
{
    return fbp;
}

double LinearProgressionParameters::get_bgp()
{
    return bgp;
}

void LinearProgressionParameters::set_fbp(double fbp)
{
    this->fbp = fbp;
}

void LinearProgressionParameters::set_bgp(double bgp)
{
    this->bgp = bgp;
}
