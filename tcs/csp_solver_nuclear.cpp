/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse 
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES 
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "csp_solver_nuclear.h"
#include "csp_solver_core.h"

#include "Ambient.h"
#include "definitions.h"

C_nuclear::C_nuclear()
{
    m_dummy_area = 42.0;
    m_mode_initial = C_csp_collector_receiver::E_csp_cr_modes::ON;
    m_mode = C_csp_collector_receiver::E_csp_cr_modes::ON;
    m_mode_prev = C_csp_collector_receiver::E_csp_cr_modes::ON;

    m_nuclear_su_delay = 0.0;
    m_nuclear_qf_delay = 0.0;
}

void C_nuclear::init()
{
	m_T_htf_hot_des += 273.15;	//[K] Convert from input in [C]
	m_T_htf_cold_des += 273.15;	//[K] Convert from input in [C]
	m_q_dot_nuc_res *= 1.E6;	    //[W] Convert from input in [MW]    
}

int C_nuclear::get_operating_state()
{
	return m_mode_prev;
}

double C_nuclear::get_startup_time()
{
    return m_nuclear_su_delay * 3600.; // sec
}

double C_nuclear::get_startup_energy()
{
    return m_nuclear_qf_delay * m_q_rec_des * 1.e-6;  // MWh
}

double C_nuclear::get_remaining_startup_energy()
{
	return m_E_su_prev;
}

void C_nuclear::call(const C_csp_weatherreader::S_outputs &weather, 
	const C_csp_solver_htf_1state &htf_state_in,
	const C_nuclear::S_inputs &inputs,
	const C_csp_solver_sim_info &sim_info)
{

}

void C_nuclear::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_solver_sim_info &sim_info)
{
    
}

void C_nuclear::converged()
{

}

void C_nuclear::calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc, double &WdotPump_calc)
{

}

double C_nuclear::get_pumping_parasitic_coef()
{

}

double C_nuclear::area_proj()
{
    return m_A_sf; //[m^2] projected or aperture area of the receiver
}
