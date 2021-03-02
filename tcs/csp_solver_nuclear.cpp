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
    m_ncall = -1;
    m_itermode = -1;
    
    m_dummy_area = 42.0;
    m_mode_initial = C_csp_collector_receiver::E_csp_cr_modes::ON;
    m_mode = C_csp_collector_receiver::E_csp_cr_modes::ON;
    m_mode_prev = C_csp_collector_receiver::E_csp_cr_modes::ON;
    
    m_T_salt_hot_target = std::numeric_limits<double>::quiet_NaN();   
    m_m_dot_htf_max     = std::numeric_limits<double>::quiet_NaN();
    m_od_control        = std::numeric_limits<double>::quiet_NaN();

    m_nuclear_su_delay = 0.0;
    m_nuclear_qf_delay = 0.0;
}

void C_nuclear::init()
{
	m_T_htf_hot_des += 273.15;	//[K] Convert from input in [C]
	m_T_htf_cold_des += 273.15;	//[K] Convert from input in [C]
	m_q_dot_nuc_res *= 1.E6;	    //[W] Convert from input in [MW] 

    m_T_salt_hot_target += 273.15;	//[K] Convert from input in [C]   
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
    return m_nuclear_qf_delay * m_q_dot_nuc_res * 1.e-6;  // MWh
}

double C_nuclear::get_remaining_startup_energy()
{
	return m_E_su_prev;
}

void C_nuclear::call(const C_csp_weatherreader::S_outputs &weather, 
	const C_csp_solver_htf_1state &htf_state_in,
	const C_pt_receiver::S_inputs &inputs,
	const C_csp_solver_sim_info &sim_info)
{
    // Increase call-per-timestep counter
	// Converge() sets it to -1, so on first call this line will adjust it = 0
	m_ncall++;
	
    // When this function is called from TCS solver, input_operation_mode should always be == 2
	C_csp_collector_receiver::E_csp_cr_modes input_operation_mode = inputs.m_input_operation_mode;

	// Get sim info 
	double step = sim_info.ms_ts.m_step;			//[s]
	double time = sim_info.ms_ts.m_time;	//[s]

	// Get applicable htf state info
	double T_salt_cold_in = htf_state_in.m_temp;		//[C]

	// Complete necessary conversions/calculations of input variables
	T_salt_cold_in += 273.15;				//[K] Cold salt inlet temp, convert from C
	double P_amb = weather.m_pres*100.0;	//[Pa] Ambient pressure, convert from mbar
	double hour = time / 3600.0;			//[hr] Hour of the year
	double T_dp = weather.m_tdew + 273.15;	//[K] Dewpoint temperature, convert from C
	double T_amb = weather.m_tdry + 273.15;	//[K] Dry bulb temperature, convert from C
	// **************************************************************************************

	// Read in remaining weather inputs from weather output structure
	double zenith = weather.m_solzen;
	double azimuth = weather.m_solazi;
	double v_wind_10 = weather.m_wspd;
	double I_bn = weather.m_beam;

	double T_sky = CSP::skytemp(T_amb, T_dp, hour);

	// Set current timestep stored values to NaN so we know that code solved for them
	m_mode = C_csp_collector_receiver::OFF;
	m_E_su = std::numeric_limits<double>::quiet_NaN();
	m_t_su = std::numeric_limits<double>::quiet_NaN();

	m_itermode = 1;

	double v_wind = log((m_h_tower + m_h_rec / 2) / 0.003) / log(10.0 / 0.003)*v_wind_10;

	double c_p_coolant, rho_coolant, f, u_coolant, q_conv_sum, q_rad_sum, q_dot_inc_sum, q_dot_piping_loss, q_dot_inc_min_panel;
	c_p_coolant = rho_coolant = f = u_coolant = q_conv_sum = q_rad_sum = q_dot_inc_sum = q_dot_piping_loss = q_dot_inc_min_panel = std::numeric_limits<double>::quiet_NaN();
	double eta_therm, m_dot_salt_tot, T_salt_hot, m_dot_salt_tot_ss, T_salt_hot_rec;
	eta_therm = m_dot_salt_tot = T_salt_hot = m_dot_salt_tot_ss = T_salt_hot_rec = std::numeric_limits<double>::quiet_NaN();
	double clearsky = std::numeric_limits<double>::quiet_NaN();
	
	bool nuc_is_off = false;
	bool nuc_is_defocusing = false;

	double T_coolant_prop = (m_T_salt_hot_target + T_salt_cold_in) / 2.0;		//[K] The temperature at which the coolant properties are evaluated. Validated as constant (mjw)
	c_p_coolant = field_htfProps.Cp(T_coolant_prop)*1000.0;						//[J/kg-K] Specific heat of the coolant

	double m_dot_htf_max = m_m_dot_htf_max;

    // need an analogous "field eff" for nuclear to suggest defocus? maybe not needed
	if (m_od_control < 1.0)
	{	// Suggests controller applied defocus, so reset *controller* defocus
		m_od_control = fmin(m_od_control, 1.0);
	}
    
	// Initialize steady state solutions with current weather, DNI, field efficiency, and inlet conditions
	s_steady_state_soln soln;
	soln.hour = time / 3600.0;
	soln.T_amb = weather.m_tdry + 273.15;
	soln.T_dp = weather.m_tdew + 273.15;
	soln.v_wind_10 = weather.m_wspd;
	soln.p_amb = weather.m_pres * 100.0;

	soln.dni = I_bn;
	soln.T_salt_cold_in = T_salt_cold_in;	
	soln.od_control = m_od_control;         // Initial defocus control (may be adjusted during the solution)
    soln.mode = input_operation_mode;
    soln.itermode = m_itermode;
	soln.nuc_is_off = nuc_is_off;
    

    //--- Solve for mass flow at actual and/or clear-sky DNI extremes
    if (use_previous_solution(soln, m_mflow_soln_prev))  // Same conditions were solved in the previous call to this method
        soln = m_mflow_soln_prev;
    else
        solve_for_mass_flow_and_defocus(soln, m_dot_htf_max, flux_map_input);

    m_mflow_soln_prev = soln;


    //--- Set mass flow and calculate final solution
    soln.q_dot_inc = m_q_dot_nuc_res;  // Absorbed flux profiles at actual DNI and clear-sky defocus
    calculate_steady_state_soln(soln, 0.00025);  // Solve energy balances at clearsky mass flow rate and actual DNI conditions

    



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


bool C_mspt_receiver_222::use_previous_solution(const s_steady_state_soln& soln, const s_steady_state_soln& soln_prev)
{
	// Are these conditions identical to those used in the last solution?
	if (!soln_prev.nuc_is_off && 
		soln.dni == soln_prev.dni &&
		soln.T_salt_cold_in == soln_prev.T_salt_cold_in &&
		soln.od_control == soln_prev.od_control &&
		soln.T_amb == soln_prev.T_amb && 
		soln.T_dp == soln_prev.T_dp &&
		soln.v_wind_10 == soln_prev.v_wind_10 &&
		soln.p_amb == soln_prev.p_amb)
	{
		return true;
	}
	else
		return false;
}
