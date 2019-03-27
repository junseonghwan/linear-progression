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

double LinearProgressionParameters::get_fbp() const
{
    return fbp;
}

double LinearProgressionParameters::get_bgp() const
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
