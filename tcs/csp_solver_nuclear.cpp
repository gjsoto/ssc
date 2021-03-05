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
    // Unit Conversions
	m_od_tube /= 1.E3;			//[m] Convert from input in [mm]
	m_th_tube /= 1.E3;			//[m] Convert from input in [mm]
	m_T_htf_hot_des += 273.15;	//[K] Convert from input in [C]
	m_T_htf_cold_des += 273.15;	//[K] Convert from input in [C]
	m_q_dot_nuc_des *= 1.E6;	    //[W] Convert from input in [MW] 

    m_T_salt_hot_target += 273.15;	//[K] Convert from input in [C]   
    m_id_tube = m_od_tube - 2 * m_th_tube;			//[m] Inner diameter of receiver tube
    
    m_mode = m_mode_initial;					//[-] 0 = requires startup, 1 = starting up, 2 = running
	m_itermode = 1;			//[-] 1: Solve for design temp, 2: solve to match mass flow restriction
	m_od_control = 1.0;			//[-] Additional defocusing for over-design conditions
	m_tol_od = 0.001;		//[-] Tolerance for over-design iteration
    
    m_m_mixed = 3.2;
    m_LoverD = m_h_rec / m_id_tube;
	m_RelRough = (4.5e-5) / m_id_tube;	//[-] Relative roughness of the tubes. http:www.efunda.com/formulae/fluids/roughness.cfm


    set_material_properties();
    
    double c_htf_des = field_htfProps.Cp((m_T_htf_hot_des + m_T_htf_cold_des) / 2.0)*1000.0;		//[J/kg-K] Specific heat at design conditions
	m_m_dot_htf_des = m_q_dot_nuc_des / (c_htf_des*(m_T_htf_hot_des - m_T_htf_cold_des));					//[kg/s]
	m_q_dot_inc_min = m_q_dot_nuc_des * m_f_rec_min ;	//[W] Minimum receiver thermal power
    
    m_m_dot_htf_max = m_m_dot_htf_max_frac * m_m_dot_htf_des;	//[kg/s]

}

void C_nuclear::set_material_properties()
{
    ambient_air.SetFluid(ambient_air.Air);

	// Declare instance of fluid class for FIELD fluid
	if( m_field_fl != HTFProperties::User_defined && m_field_fl < HTFProperties::End_Library_Fluids )
	{
		if( !field_htfProps.SetFluid( m_field_fl ) )
		{
			throw(C_csp_exception("Nuclear HTF code is not recognized", "Nuclear Island"));
		}
	}
	else if( m_field_fl == HTFProperties::User_defined )
	{

			error_msg = util::format("The user defined HTF table not implemented in Nuclear Island model.");
			throw(C_csp_exception(error_msg, "Nuclear Island"));
	}
	else
	{
		throw(C_csp_exception("Nuclear HTF code is not recognized", "Nuclear Island"));
	}

	
	// Declare instance of htf class for receiver tube material
	if( m_mat_tube == HTFProperties::Stainless_AISI316 || m_mat_tube == HTFProperties::T91_Steel ||
        m_mat_tube == HTFProperties::N06230 || m_mat_tube == HTFProperties::N07740)
	{
		if( !tube_material.SetFluid(m_mat_tube) )
		{
			throw(C_csp_exception("Tube material code not recognized", "Nuclear Island"));
		}
	}
	else if( m_mat_tube == HTFProperties::User_defined )
	{
		throw(C_csp_exception("Nuclear material currently does not accept user defined properties", "Nuclear Island"));
	}
	else
	{
		error_msg = util::format("Nuclear material code, %d, is not recognized", m_mat_tube);
		throw(C_csp_exception(error_msg, "Nuclear Island"));
	}

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
    return m_nuclear_qf_delay * m_q_dot_nuc_des * 1.e-6;  // MWh
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
	double I_bn = 1000;

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
    
    double q_thermal_ss = 0.0;


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
        solve_for_mass_flow_and_defocus(soln, m_dot_htf_max);

    m_mflow_soln_prev = soln;

    //--- Set mass flow and calculate final solution
    soln.q_dot_inc = m_q_dot_nuc_des;  // Absorbed flux profiles at actual DNI and clear-sky defocus
    calculate_steady_state_soln(soln, 0.00025);  // Solve energy balances at clearsky mass flow rate and actual DNI conditions

    // Set variables for use in the rest of the solution
	nuc_is_off = soln.nuc_is_off;
	m_mode = soln.mode;
	m_itermode = soln.itermode;
	m_od_control = soln.od_control;
    
	m_dot_salt_tot = soln.m_dot_salt_tot;
	T_salt_hot = soln.T_salt_hot;	
	T_salt_hot_rec = soln.T_salt_hot_rec;
	eta_therm = soln.eta_therm;

	u_coolant = soln.u_salt;
	f = soln.f;
	T_coolant_prop = (T_salt_hot + T_salt_cold_in) / 2.0;
	c_p_coolant = field_htfProps.Cp(T_coolant_prop)*1000.0;
	rho_coolant = field_htfProps.dens(T_coolant_prop, 1.0);

	q_conv_sum = soln.Q_conv_sum;
	q_rad_sum = soln.Q_rad_sum;
	q_dot_piping_loss = soln.Q_dot_piping_loss;
	q_dot_inc_sum = soln.Q_inc_sum;
	q_dot_inc_min_panel = soln.Q_inc_min;

	m_T_s = soln.T_s;
	m_T_panel_in = soln.T_panel_in;
	m_T_panel_out = soln.T_panel_out;
	m_T_panel_ave = soln.T_panel_ave;

	m_q_dot_conv = soln.q_dot_conv;
	m_q_dot_rad = soln.q_dot_rad;
	m_q_dot_loss = soln.q_dot_conv + soln.q_dot_rad;
	m_q_dot_abs = soln.q_dot_abs;
	m_q_dot_inc = soln.q_dot_inc;
    
	// Calculate total absorbed solar energy and minimum absorbed per panel if needed
	q_dot_inc_min_panel = q_dot_inc_sum;

	double q_thermal_steadystate = soln.Q_thermal;

	double DELTAP, Pres_D, W_dot_pump, q_thermal, q_startup;
	DELTAP = Pres_D = W_dot_pump = q_thermal = q_startup = std::numeric_limits<double>::quiet_NaN();

	q_startup = 0.0;

	double time_required_su = step/3600.0;

	if( !nuc_is_off )
	{
		m_dot_salt_tot_ss = m_dot_salt_tot;

		switch( input_operation_mode )
		{
		case C_csp_collector_receiver::STARTUP:
            {
                throw(C_csp_exception("STARTUP mode node allowed", "Nuclear Island"));
			}	

			nuc_is_off = true;
            
			break;

		case C_csp_collector_receiver::ON:
			
			if( m_E_su_prev > 0.0 || m_t_su_prev > 0.0 )
			{
				throw(C_csp_exception("Can't have excess startup energy, STARTUP mode node allowed", "Nuclear Island"));
			}
			else
			{
				m_E_su = m_E_su_prev;
				m_t_su = m_t_su_prev;
				m_mode = C_csp_collector_receiver::ON;
				q_startup = 0.0;

				q_thermal = m_dot_salt_tot*c_p_coolant*(T_salt_hot - T_salt_cold_in);

			}
			break;

		case C_csp_collector_receiver::STEADY_STATE:

            throw(C_csp_exception("STARTUP mode node allowed", "Nuclear Island"));
            
			break;
		
		}	// End switch() on input_operation_mode

		// Pressure drop calculations
        calc_pump_performance(rho_coolant, m_dot_salt_tot, f, Pres_D, W_dot_pump);

		q_thermal = m_dot_salt_tot*c_p_coolant*(T_salt_hot - T_salt_cold_in);
		q_thermal_ss = m_dot_salt_tot_ss*c_p_coolant*(T_salt_hot - T_salt_cold_in);

		// After convergence, determine whether the mass flow rate falls below the lower limit
        if (q_dot_inc_sum < m_q_dot_inc_min && (m_dot_salt_tot < m_f_rec_min * m_m_dot_htf_des))
		{
			// GOTO 900
			// Steady State always reports q_thermal (even when much less than min) because model is letting receiver begin startup with this energy
			// Should be a way to communicate to controller that q_thermal is less than q_min without losing this functionality
			throw(C_csp_exception("Mass flow rate too low", "Nuclear Island"));
		}
	}
	else
	{	// If receiver was off BEFORE startup deductions
		throw(C_csp_exception("Nuclear can't be OFF", "Nuclear Island"));
	}
    
	outputs.m_m_dot_salt_tot = m_dot_salt_tot*3600.0;		//[kg/hr] convert from kg/s
	outputs.m_eta_therm = eta_therm;							//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
	outputs.m_W_dot_pump = W_dot_pump / 1.E6;				//[MW] convert from W
	outputs.m_q_conv_sum = q_conv_sum / 1.E6;				//[MW] convert from W
	outputs.m_q_rad_sum = q_rad_sum / 1.E6;					//[MW] convert from W
	outputs.m_Q_thermal = q_thermal / 1.E6;					//[MW] convert from W
	outputs.m_T_salt_hot = T_salt_hot - 273.15;				//[C] convert from K
    outputs.m_T_salt_hot_rec = T_salt_hot_rec - 273.15;     // [C] convert from K[-]
	outputs.m_component_defocus = m_od_control;				//[-]
	outputs.m_q_dot_rec_inc = q_dot_inc_sum / 1.E6;			//[MW] convert from W
	outputs.m_q_startup = q_startup/1.E6;					//[MW-hr] convert from W-hr
	outputs.m_dP_receiver = DELTAP/ 1.E5;	//[bar] receiver pressure drop, convert from Pa
	outputs.m_dP_total = Pres_D*10.0;						//[bar] total pressure drop, convert from MPa
	outputs.m_vel_htf = u_coolant;							//[m/s]
	outputs.m_T_salt_cold = T_salt_cold_in - 273.15;			//[C] convert from K
	outputs.m_m_dot_ss = m_dot_salt_tot_ss*3600.0;			//[kg/hr] convert from kg/s
	outputs.m_q_dot_ss = q_thermal_ss / 1.E6;				//[MW] convert from W
	outputs.m_time_required_su = time_required_su*3600.0;	//[s], convert from hr in code
	if(q_thermal > 0.0)
		outputs.m_q_dot_piping_loss = q_dot_piping_loss/1.E6;	//[MWt]
	else
		outputs.m_q_dot_piping_loss = 0.0;		//[MWt]
    outputs.m_q_heattrace = 0.0;

	outputs.m_clearsky = clearsky;  // W/m2
	outputs.m_Q_thermal_ss = q_thermal_steadystate / 1.e6; //[MWt]

    ms_outputs = outputs;

}

void C_nuclear::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_solver_sim_info &sim_info)
{
    
}

void C_nuclear::converged()
{
	// Check HTF props?
	//!MJW 9.8.2010 :: Call the property range check subroutine with the inlet and outlet HTF temps to make sure they're in the valid range
	//call check_htf(Coolant,T_salt_hot)
	//call check_htf(Coolant,T_salt_cold)

	if( m_mode == C_csp_collector_receiver::STEADY_STATE )
	{
		throw(C_csp_exception("Receiver should only be run at STEADY STATE mode for estimating output. It must be run at a different mode before exiting a timestep",
			"MSPT receiver converged method"));
	}

	if( m_mode == C_csp_collector_receiver::OFF )
	{
		m_E_su = m_q_rec_des * m_rec_qf_delay;
		m_t_su = m_rec_su_delay;
	}

	m_mode_prev = m_mode;
	m_E_su_prev = m_E_su;
	m_t_su_prev = m_t_su;

	m_itermode = 1;
	m_od_control = 1.0;

	m_ncall = -1;

    ms_outputs = outputs;
}


void C_nuclear::calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc, double &WdotPump_calc)
{
    // Pressure drop calculations
	double u_coolant = mdot / (rho_f * m_id_tube * m_id_tube * 0.25 * CSP::pi);	//[m/s] Average velocity of the coolant through the receiver tubes
    
    double fudge = 0;
	double L_e_45 = 16.0;						// The equivalent length produced by the 45 degree bends in the tubes - Into to Fluid Mechanics, Fox et al.
	double L_e_90 = 30.0;						// The equivalent length produced by the 90 degree bends in the tubes
	double DELTAP_tube = rho_f*(ffact*m_h_rec / m_id_tube*pow(u_coolant, 2) / 2.0);	//[Pa] Pressure drop across the tube, straight length
	double DELTAP_45 = rho_f*(ffact*L_e_45*pow(u_coolant, 2) / 2.0);					//[Pa] Pressure drop across 45 degree bends
	double DELTAP_90 = rho_f*(ffact*L_e_90*pow(u_coolant, 2) / 2.0);					//[Pa] Pressure drop across 90 degree bends
	double DELTAP = DELTAP_tube + 2 * DELTAP_45 + 4 * DELTAP_90;						//[Pa] Total pressure drop across the tube with (4) 90 degree bends, (2) 45 degree bends
	double DELTAP_h_tower = rho_f*m_h_tower*CSP::grav;						//[Pa] The pressure drop from pumping up to the receiver
	double DELTAP_net = fudge * (DELTAP + DELTAP_h_tower);		//[Pa] The new pressure drop across the receiver panels
	PresDrop_calc = DELTAP_net*1.E-6;			//[MPa]
	double est_load = fmax(0.25, mdot / m_m_dot_htf_des) * 100;		//[%] Relative pump load. Limit to 25%
	double eta_pump_adj = m_eta_pump*(-2.8825E-9*pow(est_load, 4) + 6.0231E-7*pow(est_load, 3) - 1.3867E-4*pow(est_load, 2) + 2.0683E-2*est_load);	//[-] Adjusted pump efficiency
	WdotPump_calc = DELTAP_net*mdot / rho_f / eta_pump_adj;
}

double C_nuclear::get_pumping_parasitic_coef()
{
    double Tavg = (m_T_htf_cold_des + m_T_htf_hot_des) / 2.;

    double mu_coolant = field_htfProps.visc(Tavg);				//[kg/m-s] Absolute viscosity of the coolant
    double k_coolant = field_htfProps.cond(Tavg);				//[W/m-K] Conductivity of the coolant
    double rho_coolant = field_htfProps.dens(Tavg, 1.0);        //[kg/m^3] Density of the coolant
    double c_p_coolant = field_htfProps.Cp(Tavg)*1e3;           //[J/kg-K] Specific heat

    double m_dot_salt = m_q_dot_nuc_des / (c_p_coolant * (m_T_htf_hot_des - m_T_htf_cold_des));

    double n_t = (int)(CSP::pi*m_d_rec / (m_od_tube));   // The number of tubes per panel, as a function of the number of panels and the desired diameter of the receiver
    double id_tube = m_od_tube - 2 * m_th_tube;                 //[m] Inner diameter of receiver tube


    double u_coolant = m_dot_salt / (n_t*rho_coolant*pow((id_tube / 2.0), 2)*CSP::pi);	//[m/s] Average velocity of the coolant through the receiver tubes
    double Re_inner = rho_coolant * u_coolant*id_tube / mu_coolant;				        //[-] Reynolds number of internal flow
    double Pr_inner = c_p_coolant * mu_coolant / k_coolant;						        //[-] Prandtl number of internal flow
    double Nusselt_t, f;
    double LoverD = m_h_rec / id_tube;
    double RelRough = (4.5e-5) / id_tube;   //[-] Relative roughness of the tubes. http:www.efunda.com/formulae/fluids/roughness.cfm
    CSP::PipeFlow(Re_inner, Pr_inner, LoverD, RelRough, Nusselt_t, f);

    double deltap, wdot;
    calc_pump_performance(rho_coolant, m_dot_salt, f, deltap, wdot);

    return wdot / m_q_dot_nuc_des;
}

double C_nuclear::area_proj()
{
    return m_A_sf; //[m^2] projected or aperture area of the receiver
}

bool C_nuclear::use_previous_solution(const s_steady_state_soln& soln, const s_steady_state_soln& soln_prev)
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

// Calculate mass flow rate and defocus needed to achieve target outlet temperature given DNI
void C_nuclear::solve_for_mass_flow_and_defocus(s_steady_state_soln &soln, double m_dot_htf_max)
{
	bool nuc_is_defocusing = true;
	double err_od = 999.0;

	while (nuc_is_defocusing)
	{
		if (soln.nuc_is_off)
			break;

		soln.q_dot_inc = m_q_dot_nuc_des;  // Calculate flux profiles
		solve_for_mass_flow(soln);	// Iterative calculation of mass flow to produce target outlet temperature

		if (soln.nuc_is_off)
			break;

        double m_dot_salt_tot = soln.m_dot_salt_tot;

		// Limit the HTF mass flow rate to the maximum, if needed
		nuc_is_defocusing = false;
		if ((m_dot_salt_tot > m_dot_htf_max) || soln.itermode == 2)
		{
			double err_od = (m_dot_salt_tot - m_dot_htf_max) / m_dot_htf_max;
			if (err_od < m_tol_od)
			{
				soln.itermode = 1;
				soln.od_control = 1.0;
				nuc_is_defocusing = false;
			}
			else
			{
				soln.od_control = soln.od_control * pow((m_dot_htf_max / m_dot_salt_tot), 0.8);	//[-] Adjust the over-design defocus control by modifying the current value
				soln.itermode = 2;
				nuc_is_defocusing = true;
			}
		}
	}

	return;
}

// Calculate mass flow rate needed to achieve target outlet temperature (m_T_salt_hot_target) given incident flux profiles
void C_nuclear::solve_for_mass_flow(s_steady_state_soln &soln)
{

	bool soln_exists = (soln.m_dot_salt_tot == soln.m_dot_salt_tot);

	soln.T_salt_props = (m_T_salt_hot_target + soln.T_salt_cold_in) / 2.0;		//[K] The temperature at which the coolant properties are evaluated. Validated as constant (mjw)
	double c_p_coolant = field_htfProps.Cp(soln.T_salt_props)*1000.0;				//[J/kg-K] Specific heat of the coolant

	double m_dot_salt_guess;
	if (soln_exists)  // Use existing solution as intial guess
	{
		m_dot_salt_guess = soln.m_dot_salt_tot;
	}
	else  // Set inital guess for mass flow solution
	{
        
		double q_dot_inc_sum = soln.q_dot_inc;		//[kW] Total power absorbed by receiver

		double c_guess = field_htfProps.Cp((m_T_salt_hot_target + soln.T_salt_cold_in) / 2.0)*1000.;	//[kJ/kg-K] Estimate the specific heat of the fluid in receiver

        m_dot_salt_guess = q_dot_inc_sum / (c_guess*(m_T_salt_hot_target - soln.T_salt_cold_in));	//[kg/s] Mass flow rate for each flow path			

	}


	// Set solution tolerance
	double T_salt_hot_guess = 9999.9;		//[K] Initial guess value for error calculation
	double err = -999.9;					//[-] Relative outlet temperature error
	double tol = std::numeric_limits<double>::quiet_NaN();
	tol = 0.00025;

	int qq_max = 50;
	int qq = 0;
    double m_dot_salt_tot;
	bool converged = false;
	while (!converged)
	{
		qq++;

		// if the problem fails to converge after 50 iterations, then the power is likely negligible and the zero set can be returned
		if (qq > qq_max)
		{
			soln.mode = C_csp_collector_receiver::OFF;  // Set the startup mode
			soln.nuc_is_off = true;
			break;
		}

		soln.m_dot_salt_tot = m_dot_salt_guess;
		double tolT = tol;
		calculate_steady_state_soln(soln, tolT, 50);   // Solve steady state thermal model

        double err_per_path;
        err = -999.9;
        m_dot_salt_tot = m_dot_salt_guess;
        err_per_path = (soln.T_salt_hot - m_T_salt_hot_target) / m_T_salt_hot_target;
        err = fmax(err, fabs(err_per_path));

		if (soln.nuc_is_off)  // SS solution was unsuccessful or resulted in an infeasible exit temperature -> remove outlet T for solution to start next iteration from the default intial guess
			soln.T_salt_hot = std::numeric_limits<double>::quiet_NaN();

		if (fabs(err) > tol)
		{
            double m_dot_salt_min = 1.e10;
            m_dot_salt_guess = (soln.Q_abs_sum - soln.Q_dot_piping_loss) / (c_p_coolant * (m_T_salt_hot_target - soln.T_salt_cold_in));
            m_dot_salt_min = fmin(m_dot_salt_min, m_dot_salt_guess);

			if (m_dot_salt_min < 1.E-5)
			{
				soln.mode = C_csp_collector_receiver::OFF;
				soln.nuc_is_off = true;
				break;
			}
		}
        else
        {
            converged = true;
            if (err_per_path > 0.0)   // Solution has converged but outlet T is above target.  CSP solver seems to perform better with slighly under-design temperature than with slighly over-design temperatures.
            {
                converged = false;
                m_dot_salt_guess *= (soln.T_salt_hot - soln.T_salt_cold_in) / ((1.0 - 0.5 * tol) * m_T_salt_hot_target - soln.T_salt_cold_in);

            }
        }
	}

	soln.m_dot_salt_tot = m_dot_salt_tot;

	return;
}

// Calculate steady state temperature and heat loss profiles for a given mass flow and incident flux
void C_nuclear::calculate_steady_state_soln(s_steady_state_soln &soln, double tol, int max_iter)
{

	double P_amb = soln.p_amb;	
	double hour = soln.hour;			
	double T_dp = soln.T_dp;	
	double T_amb = soln.T_amb;	
	double v_wind_10 = soln.v_wind_10;
	double T_sky = CSP::skytemp(T_amb, T_dp, hour);
	double v_wind = log((m_h_tower + m_h_rec / 2) / 0.003) / log(10.0 / 0.003)*v_wind_10;

	double T_s_guess;
	double T_panel_out_guess;
	double T_panel_in_guess;
	double T_film;

	bool soln_exists = (soln.T_salt_hot == soln.T_salt_hot);

	// Set initial guess
	double T_salt_hot_guess;
	if (soln_exists)    // Use existing solution as the inital guess
	{
		T_salt_hot_guess = soln.T_salt_hot;    // Initial guess for outlet T
		T_s_guess = soln.T_s;
		T_panel_out_guess = soln.T_panel_out;
		T_panel_in_guess = soln.T_panel_in;
	}
	else // Initialize solution from scratch
	{
		T_salt_hot_guess = m_T_salt_hot_target;    // Initial guess for outlet T
        T_s_guess = m_T_salt_hot_target;			//[K] Guess the temperature for the surface nodes
        T_panel_out_guess = soln.T_salt_cold_in;	//[K] Guess values for the fluid temp coming out of the control volume
        T_panel_in_guess = soln.T_salt_cold_in;		//[K] Guess values for the fluid temp coming into the control volume
	}


	// Temperature solution iterations
	for (int q = 0; q < max_iter; q++)
	{
		double T_coolant_prop;
		if (soln.T_salt_props == soln.T_salt_props)   // Temperature for property evaluation exists (calling from within loop over salt mass flow rate)
			T_coolant_prop = soln.T_salt_props;
		else    
			T_coolant_prop = (T_salt_hot_guess + soln.T_salt_cold_in) / 2.0;
		double c_p_coolant = field_htfProps.Cp(T_coolant_prop)*1000.;	

		soln.T_s = T_s_guess;
        soln.T_panel_out= T_panel_out_guess;
        soln.T_panel_in = T_panel_in_guess;
        soln.T_panel_ave = (soln.T_panel_in + soln.T_panel_out) / 2.0;		//[K] The average coolant temperature in each control volume
        T_film = (soln.T_s+ T_amb) / 2.0;									//[K] Film temperature

		// Calculate the average surface temperature
		double T_film_ave = (T_amb + T_salt_hot_guess) / 2.0;

        
		// Convective coefficient for external forced convection using Siebers & Kraabel
		double k_film = ambient_air.cond(T_film_ave);				//[W/m-K] The conductivity of the ambient air
		double mu_film = ambient_air.visc(T_film_ave);				//[kg/m-s] Dynamic viscosity of the ambient air
		double rho_film = ambient_air.dens(T_film_ave, P_amb);		//[kg/m^3] Density of the ambient air
		double c_p_film = ambient_air.Cp(T_film_ave);				//[kJ/kg-K] Specific heat of the ambient air
		double Re_for = rho_film * v_wind*m_d_rec / mu_film;		//[-] Reynolds number
		double ksD = (m_od_tube / 2.0) / m_d_rec;					//[-] The effective roughness of the cylinder [Siebers, Kraabel 1984]
		double Nusselt_for = CSP::Nusselt_FC(ksD, Re_for);			//[-] S&K
		double h_for = Nusselt_for * k_film / m_d_rec * m_hl_ffact;	//[W/m^2-K] Forced convection heat transfer coefficient

		// Convection coefficient for external natural convection using Siebers & Kraabel
		// Note: This relationship applies when the surrounding properties are evaluated at ambient conditions [S&K]
		double beta = 1.0 / T_amb;													//[1/K] Volumetric expansion coefficient
		double nu_amb = ambient_air.visc(T_amb) / ambient_air.dens(T_amb, P_amb);	//[m^2/s] Kinematic viscosity		

        //====== Here there was a loop. 
        
        // Natural convection
        double Gr_nat = fmax(0.0, CSP::grav*beta*(soln.T_s - T_amb)*pow(m_h_rec, 3) / pow(nu_amb, 2));	//[-] Grashof Number at ambient conditions
        double Nusselt_nat = 0.098*pow(Gr_nat, (1.0 / 3.0))*pow(soln.T_s / T_amb, -0.14);				//[-] Nusselt number
        double h_nat = Nusselt_nat * ambient_air.cond(T_amb) / m_h_rec * m_hl_ffact;							//[W/m^-K] Natural convection coefficient

        // Mixed convection
        double h_mixed = pow((pow(h_for, m_m_mixed) + pow(h_nat, m_m_mixed)), 1.0 / m_m_mixed)*4.0;		//(4.0) is a correction factor to match convection losses at Solar II (correspondance with G. Kolb, SNL)
        
        
        soln.q_dot_conv = 0.0;			//[W] Convection losses per node

        // Radiation from the receiver - Calculate the radiation node by node
        soln.q_dot_rad = 0.0;	//[W] Total radiation losses per node
        soln.q_dot_loss = soln.q_dot_rad + soln.q_dot_conv;			//[W] Total overall losses per node
        soln.q_dot_abs = soln.q_dot_inc - soln.q_dot_loss;			//[W] Absorbed flux at each node
        
        
        // Calculate the temperature drop across the receiver tube wall... assume a cylindrical thermal resistance
        double T_wall = (soln.T_s + soln.T_panel_ave) / 2.0;				//[K] The temperature at which the conductivity of the wall is evaluated
        double k_tube = tube_material.cond(T_wall);											//[W/m-K] The conductivity of the wall
        double R_tube_wall = m_th_tube / (k_tube*m_h_rec*m_d_rec*pow(CSP::pi, 2) / 2.0 );	//[K/W] The thermal resistance of the wall

        // Calculations for the inside of the tube						
        double mu_coolant = field_htfProps.visc(T_coolant_prop);							//[kg/m-s] Absolute viscosity of the coolant
        double k_coolant = field_htfProps.cond(T_coolant_prop);								//[W/m-K] Conductivity of the coolant
        double rho_coolant = field_htfProps.dens(T_coolant_prop, 1.0);						//[kg/m^3] Density of the coolant
        double u_coolant = soln.m_dot_salt_tot / (rho_coolant*pow((m_id_tube / 2.0), 2)*CSP::pi);	//[m/s] Average velocity of the coolant through the receiver tubes
        double Re_inner = rho_coolant * u_coolant*m_id_tube / mu_coolant;					//[-] Reynolds number of internal flow
        double Pr_inner = c_p_coolant * mu_coolant / k_coolant;								//[-] Prandtl number of internal flow

        double Nusselt_t, f;
        CSP::PipeFlow(Re_inner, Pr_inner, m_LoverD, m_RelRough, Nusselt_t, f);
        if (Nusselt_t <= 0.0)
        {
            soln.mode = C_csp_collector_receiver::OFF;	
            break;
        }
        double h_inner = Nusselt_t * k_coolant / m_id_tube;								//[W/m^2-K] Convective coefficient between the inner tube wall and the coolant
        double R_conv_inner = 1.0 / (h_inner*CSP::pi*m_id_tube / 2.0*m_h_rec);	//[K/W] Thermal resistance associated with this value

        soln.u_salt = u_coolant;
        soln.f = f;
        
        
        // Update panel inlet/outlet temperature guess
        T_panel_in_guess = soln.T_salt_cold_in;
        
        T_panel_out_guess = T_panel_in_guess + soln.q_dot_abs/ (soln.m_dot_salt_tot*c_p_coolant);		//[K] Energy balance for each node		
        //double Tavg = (T_panel_out_guess + T_panel_in_guess) / 2.0;											//[K] Panel average temperature
        //T_s_guess = Tavg + soln.q_dot_abs*(R_conv_inner + R_tube_wall);										//[K] Surface temperature based on the absorbed heat
        /*
        if (T_s_guess < 1.0)
        {
            soln.mode = C_csp_collector_receiver::OFF;
            break;
        }
        */

		//======== Here ended the giant loop

		// Calculate average receiver outlet temperature
		soln.T_salt_hot = T_panel_out_guess;  // no need for Mass-weighted average exit temperature
        soln.T_salt_hot_rec = soln.T_salt_hot;


		// Calculate outlet temperature after piping losses
        /*
		soln.Q_dot_piping_loss = 0.0;
        soln.delta_T_piping = 0.0;
		if (m_Q_dot_piping_loss > 0.0)
		{
			if (m_piping_loss_coeff != m_piping_loss_coeff)   // Calculate piping loss from constant loss per m (m_Q_dot_piping_loss)
				soln.Q_dot_piping_loss = m_Q_dot_piping_loss;
			else
			{
				double riser_loss = 2.0*CSP::pi * (soln.T_salt_cold_in - T_amb) / m_Rtot_riser; //[W/m]
				double downc_loss = 2.0*CSP::pi * (soln.T_salt_hot - T_amb) / m_Rtot_downc; //[W/m]
				soln.Q_dot_piping_loss = 0.5*(riser_loss + downc_loss) * (m_h_tower*m_pipe_length_mult + m_pipe_length_add); // Total piping thermal loss [W]
			}
			soln.delta_T_piping = soln.Q_dot_piping_loss / (m_dot_salt_tot*c_p_coolant);	//[K]
			soln.T_salt_hot -= soln.delta_T_piping;	//[K]
		}
		*/
		
		// Check convergence
		double err = (soln.T_salt_hot - T_salt_hot_guess) / T_salt_hot_guess;
		T_salt_hot_guess = soln.T_salt_hot;
		if (fabs(err) < tol && q>0)
			break;

	} // End iterations

    // Removing this to allow outlet < inlet for receiver standby (clear-sky control with little DNI)
	//if (soln.T_salt_hot < soln.T_salt_cold_in)  
	//	soln.mode = C_csp_collector_receiver::OFF;


	// Save overall energy loss (can delete later)
    soln.Q_inc_sum = soln.q_dot_inc;
    soln.Q_conv_sum = soln.q_dot_conv;
    soln.Q_rad_sum = soln.q_dot_rad;
    soln.Q_abs_sum = soln.q_dot_abs;
    soln.Q_inc_min = fmin(soln.Q_inc_min, soln.q_dot_inc);

	soln.Q_thermal = soln.Q_abs_sum - soln.Q_dot_piping_loss;

	if (soln.Q_inc_sum > 0.0)
		soln.eta_therm = soln.Q_abs_sum / soln.Q_inc_sum;
	else
		soln.eta_therm = 0.0;

	soln.nuc_is_off = false;
	if (soln.mode == C_csp_collector_receiver::OFF)
		soln.nuc_is_off = true;

	// Save final temperature profile solution
	if (!soln.nuc_is_off)
	{
		soln.T_s = T_s_guess;
		soln.T_panel_out = T_panel_out_guess;
		soln.T_panel_in = T_panel_in_guess;
		soln.T_panel_ave = (soln.T_panel_in + soln.T_panel_out) / 2.0;
	}

	return;

}
