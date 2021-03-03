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

#include "csp_solver_nuclear_plant.h"
#include "sam_csp_util.h"
#include <algorithm>

static C_csp_reported_outputs::S_output_info S_output_info[] =
{
	{C_csp_nuclear_plant::E_FIELD_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_FIELD_ETA_OPT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_FIELD_ADJUST, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    {C_csp_nuclear_plant::E_IS_FIELD_TRACKING_FINAL, C_csp_reported_outputs::TS_LAST},

	{C_csp_nuclear_plant::E_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_ETA_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_Q_DOT_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_M_DOT_HTF, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_Q_DOT_STARTUP, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_T_HTF_IN, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_T_HTF_OUT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    {C_csp_nuclear_plant::E_T_HTF_OUT_REC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_Q_DOT_PIPE_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_nuclear_plant::E_Q_DOT_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    // from transient model:
	{ C_csp_nuclear_plant::E_P_HEATTRACE, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_nuclear_plant::E_T_HTF_OUT_END, C_csp_reported_outputs::TS_LAST },
    { C_csp_nuclear_plant::E_T_HTF_OUT_END_REC, C_csp_reported_outputs::TS_LAST },
	{ C_csp_nuclear_plant::E_T_HTF_OUT_MAX, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_nuclear_plant::E_T_HTF_PANEL_OUT_MAX, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_nuclear_plant::E_T_WALL_INLET, C_csp_reported_outputs::TS_LAST },
	{ C_csp_nuclear_plant::E_T_WALL_OUTLET, C_csp_reported_outputs::TS_LAST },
	{ C_csp_nuclear_plant::E_T_RISER, C_csp_reported_outputs::TS_LAST },
	{ C_csp_nuclear_plant::E_T_DOWNC, C_csp_reported_outputs::TS_LAST },

	{ C_csp_nuclear_plant::E_CLEARSKY, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_nuclear_plant::E_Q_DOT_THERMAL_CSKY_SS, C_csp_reported_outputs::TS_WEIGHTED_AVE },
    { C_csp_nuclear_plant::E_Q_DOT_THERMAL_SS, C_csp_reported_outputs::TS_WEIGHTED_AVE },
    { C_csp_nuclear_plant::E_REC_OP_MODE_FINAL, C_csp_reported_outputs::TS_LAST},
    { C_csp_nuclear_plant::E_REC_STARTUP_TIME_REMAIN_FINAL, C_csp_reported_outputs::TS_LAST },
    { C_csp_nuclear_plant::E_REC_STARTUP_ENERGY_REMAIN_FINAL, C_csp_reported_outputs::TS_LAST },
	csp_info_invalid	
};

C_csp_nuclear_plant::C_csp_nuclear_plant(C_nuclear &nuclear):
    mc_nuclear(nuclear)
{
	mc_reported_outputs.construct(S_output_info);
}

C_csp_nuclear_plant::~C_csp_nuclear_plant()
{}

void C_csp_nuclear_plant::init(const C_csp_collector_receiver::S_csp_cr_init_inputs init_inputs, 
				C_csp_collector_receiver::S_csp_cr_solved_params & solved_params)
{
    mc_nuclear.init();
    solved_params.m_A_aper_total = mc_nuclear.m_dummy_area;	//[m^2]
    solved_params.m_T_htf_cold_des = mc_nuclear.m_T_htf_cold_des;       //[K]
    solved_params.m_T_htf_hot_des = mc_nuclear.m_T_htf_hot_des;         //[K]
    solved_params.m_q_dot_rec_des = mc_nuclear.m_q_rec_des / 1.E6;		//[MW]
	return;
}

int C_csp_nuclear_plant::get_operating_state()
{
	return mc_nuclear.get_operating_state();
}

double C_csp_nuclear_plant::get_startup_time()
{
    return mc_nuclear.get_startup_time();   //[s]
}

double C_csp_nuclear_plant::get_startup_energy()
{
    return mc_nuclear.get_startup_energy(); //[MWh]
}

double C_csp_nuclear_plant::get_pumping_parasitic_coef()  //MWe/MWt
{
    return mc_nuclear.get_pumping_parasitic_coef();
}

double C_csp_nuclear_plant::get_min_power_delivery()    //MWt
{
    return mc_nuclear.m_f_rec_min * mc_nuclear.m_q_dot_nuc_des*1.e-6;
}


double C_csp_nuclear_plant::get_tracking_power()
{
	return 0.0;	//MWe
}

double C_csp_nuclear_plant::get_col_startup_power()
{
	return 0.0;	//MWe-hr
}

double C_csp_nuclear_plant::get_remaining_startup_energy()
{
	//return mc_mspt_receiver_222.get_remaining_startup_energy()/1000.;  //kWjt
    return 0.0;  //kWjt
}


void C_csp_nuclear_plant::call(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_collector_receiver::S_csp_cr_inputs &inputs,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
	// What about catching errors here?
	

	// Get heliostat field outputs and set corresponding receiver inputs
	C_pt_receiver::S_inputs receiver_inputs;

	receiver_inputs.m_input_operation_mode = inputs.m_input_operation_mode;
	mc_nuclear.call(weather, htf_state_in, receiver_inputs, sim_info);
		

	cr_out_solver.m_q_thermal = mc_nuclear.ms_outputs.m_Q_thermal;				//[MW]
	cr_out_solver.m_q_startup = mc_nuclear.ms_outputs.m_q_startup;				//[MWt-hr]
    
	cr_out_solver.m_m_dot_salt_tot = mc_nuclear.ms_outputs.m_m_dot_salt_tot;		//[kg/hr]
	cr_out_solver.m_T_salt_hot = mc_nuclear.ms_outputs.m_T_salt_hot;				//[C]
	
	cr_out_solver.m_component_defocus = mc_nuclear.ms_outputs.m_component_defocus;	//[-]
	
	cr_out_solver.m_W_dot_htf_pump = mc_nuclear.ms_outputs.m_W_dot_pump;			//[MWe]
	cr_out_solver.m_W_dot_col_tracking = 0.0;		//[MWe]

	cr_out_solver.m_time_required_su = 0.0;	//[s]
	cr_out_solver.m_q_rec_heattrace = 0.0;		//[MWt])


	mc_reported_outputs.value(E_FIELD_Q_DOT_INC, 0.0);	//[MWt]
	mc_reported_outputs.value(E_FIELD_ETA_OPT, 0.0);			//[-]
	mc_reported_outputs.value(E_FIELD_ADJUST, 0.0);			//[-]

	mc_reported_outputs.value(E_Q_DOT_INC, mc_nuclear.ms_outputs.m_q_dot_rec_inc);	//[MWt]
	mc_reported_outputs.value(E_ETA_THERMAL, mc_nuclear.ms_outputs.m_eta_therm);		//[-]
	mc_reported_outputs.value(E_Q_DOT_THERMAL, mc_nuclear.ms_outputs.m_Q_thermal);	//[MWt]
	mc_reported_outputs.value(E_M_DOT_HTF, mc_nuclear.ms_outputs.m_m_dot_salt_tot);	//[kg/hr]
		// If startup, then timestep may have changed (why not report this from 222 in MWt?)
	mc_reported_outputs.value(E_Q_DOT_STARTUP, 0.0);		//[MWt])
	mc_reported_outputs.value(E_T_HTF_IN, htf_state_in.m_temp);									//[C]
	mc_reported_outputs.value(E_T_HTF_OUT, mc_nuclear.ms_outputs.m_T_salt_hot);		//[C]
    mc_reported_outputs.value(E_T_HTF_OUT_REC, mc_nuclear.ms_outputs.m_T_salt_hot_rec);		//[C]
	mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, mc_nuclear.ms_outputs.m_q_dot_piping_loss);	//[MWt]
    mc_reported_outputs.value(E_Q_DOT_LOSS, mc_nuclear.ms_outputs.m_q_rad_sum + mc_nuclear.ms_outputs.m_q_conv_sum ); //MWt
    // from transient model:
	mc_reported_outputs.value(E_P_HEATTRACE, 0.0);		//[MWt])
	mc_reported_outputs.value(E_T_HTF_OUT_END, mc_nuclear.ms_outputs.m_inst_T_salt_hot);	//[C]
    mc_reported_outputs.value(E_T_HTF_OUT_END_REC, mc_nuclear.ms_outputs.m_inst_T_salt_hot_rec);		//[C]
	mc_reported_outputs.value(E_T_HTF_OUT_MAX, mc_nuclear.ms_outputs.m_max_T_salt_hot);	//[C]
	mc_reported_outputs.value(E_T_HTF_PANEL_OUT_MAX, mc_nuclear.ms_outputs.m_max_rec_tout);	//[C]
	
	mc_reported_outputs.value(E_T_WALL_INLET, mc_nuclear.ms_outputs.m_Twall_inlet);	//[C]
	mc_reported_outputs.value(E_T_WALL_OUTLET, mc_nuclear.ms_outputs.m_Twall_outlet);	//[C]
	mc_reported_outputs.value(E_T_RISER, mc_nuclear.ms_outputs.m_Triser);	//[C]
	mc_reported_outputs.value(E_T_DOWNC, mc_nuclear.ms_outputs.m_Tdownc);	//[C]

	mc_reported_outputs.value(E_CLEARSKY, mc_nuclear.ms_outputs.m_clearsky);
	mc_reported_outputs.value(E_Q_DOT_THERMAL_CSKY_SS, mc_nuclear.ms_outputs.m_Q_thermal_csky_ss); //[MWt]
	mc_reported_outputs.value(E_Q_DOT_THERMAL_SS, mc_nuclear.ms_outputs.m_Q_thermal_ss); //[MWt]
}

void C_csp_nuclear_plant::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
	// First call heliostat field class
	// In OFF call, looking specifically for weather STOW parasitics apply

	// Set collector/receiver parent class outputs from field model
	//cr_out_report.m_eta_field = mc_pt_heliostatfield.ms_outputs.m_eta_field;				//[-]
    //cr_out_report.m_sf_adjust_out = mc_pt_heliostatfield.ms_outputs.m_sf_adjust_out;
	//cr_out_report.m_q_dot_field_inc = mc_pt_heliostatfield.ms_outputs.m_q_dot_field_inc;	//[MWt]
	cr_out_solver.m_W_dot_col_tracking = 0.0;			//[MWe]

	// Now, call the tower-receiver model
	mc_nuclear.off(weather, htf_state_in, sim_info);

	// Set collector/receiver parent class outputs from field model
	//cr_out_report.m_q_dot_rec_inc = mc_heat_input.ms_outputs.m_q_dot_rec_inc;		 //[MWt]
	//cr_out_report.m_eta_thermal = mc_heat_input.ms_outputs.m_eta_therm;				 //[-]
	cr_out_solver.m_q_thermal = mc_nuclear.ms_outputs.m_Q_thermal;				 //[MW]
	cr_out_solver.m_q_startup = mc_nuclear.ms_outputs.m_q_startup;				 //[MWt-hr]
	//cr_out_report.m_q_dot_piping_loss = mc_heat_input.ms_outputs.m_q_dot_piping_loss; //[MWt]
	cr_out_solver.m_m_dot_salt_tot = mc_nuclear.ms_outputs.m_m_dot_salt_tot;		 //[kg/hr]
	cr_out_solver.m_T_salt_hot = mc_nuclear.ms_outputs.m_T_salt_hot;				 //[C]
	cr_out_solver.m_component_defocus = 1.0;	//[-]
	cr_out_solver.m_W_dot_htf_pump = mc_nuclear.ms_outputs.m_W_dot_pump;			 //[MWe]
		// Not sure that we want 'startup time required' calculated in 'off' call
	cr_out_solver.m_time_required_su = mc_nuclear.ms_outputs.m_time_required_su;	 //[s]
	cr_out_solver.m_q_rec_heattrace = mc_nuclear.ms_outputs.m_q_heattrace / (mc_nuclear.ms_outputs.m_time_required_su / 3600.0);		//[MWt])

	mc_reported_outputs.value(E_FIELD_Q_DOT_INC, 0.0);	//[MWt]
	mc_reported_outputs.value(E_FIELD_ETA_OPT, 0.0);			//[-]
	mc_reported_outputs.value(E_FIELD_ADJUST, 0.0);			//[-]

	mc_reported_outputs.value(E_Q_DOT_INC, mc_nuclear.ms_outputs.m_q_dot_rec_inc);	//[MWt]
	mc_reported_outputs.value(E_ETA_THERMAL, mc_nuclear.ms_outputs.m_eta_therm);		//[-]
	mc_reported_outputs.value(E_Q_DOT_THERMAL, mc_nuclear.ms_outputs.m_Q_thermal);	//[MWt]
	mc_reported_outputs.value(E_M_DOT_HTF, mc_nuclear.ms_outputs.m_m_dot_salt_tot);	//[kg/hr]
		// Should not be startup energy in OFF, but timestep may be subhourly/nonuniform (why not report this from 222 in MWt?)
	mc_reported_outputs.value(E_Q_DOT_STARTUP, mc_nuclear.ms_outputs.m_q_startup / (mc_nuclear.ms_outputs.m_time_required_su / 3600.0));		//[MWt])
	mc_reported_outputs.value(E_T_HTF_IN, htf_state_in.m_temp);									//[C]
	mc_reported_outputs.value(E_T_HTF_OUT, mc_nuclear.ms_outputs.m_T_salt_hot);		//[C]
    mc_reported_outputs.value(E_T_HTF_OUT_REC, mc_nuclear.ms_outputs.m_T_salt_hot_rec);		//[C]
	mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, mc_nuclear.ms_outputs.m_q_dot_piping_loss);	//[MWt]
    mc_reported_outputs.value(E_Q_DOT_LOSS, mc_nuclear.ms_outputs.m_q_rad_sum + mc_nuclear.ms_outputs.m_q_conv_sum ); //MWt
    // from transient model:
	mc_reported_outputs.value(E_P_HEATTRACE, mc_nuclear.ms_outputs.m_q_heattrace / (mc_nuclear.ms_outputs.m_time_required_su / 3600.0));		//[MWt])
	mc_reported_outputs.value(E_T_HTF_OUT_END, mc_nuclear.ms_outputs.m_inst_T_salt_hot);	//[C]
    mc_reported_outputs.value(E_T_HTF_OUT_END_REC, mc_nuclear.ms_outputs.m_inst_T_salt_hot_rec);		//[C]
	mc_reported_outputs.value(E_T_HTF_OUT_MAX, mc_nuclear.ms_outputs.m_max_T_salt_hot);	//[C]
	mc_reported_outputs.value(E_T_HTF_PANEL_OUT_MAX, mc_nuclear.ms_outputs.m_max_rec_tout);	//[C]

	mc_reported_outputs.value(E_T_WALL_INLET, mc_nuclear.ms_outputs.m_Twall_inlet);	//[C]
	mc_reported_outputs.value(E_T_WALL_OUTLET, mc_nuclear.ms_outputs.m_Twall_outlet);	//[C]
	mc_reported_outputs.value(E_T_RISER, mc_nuclear.ms_outputs.m_Triser);	//[C]
	mc_reported_outputs.value(E_T_DOWNC, mc_nuclear.ms_outputs.m_Tdownc);	//[C]

	mc_reported_outputs.value(E_CLEARSKY, mc_nuclear.ms_outputs.m_clearsky);
	mc_reported_outputs.value(E_Q_DOT_THERMAL_CSKY_SS, mc_nuclear.ms_outputs.m_Q_thermal_csky_ss); //[MWt]
	mc_reported_outputs.value(E_Q_DOT_THERMAL_SS, mc_nuclear.ms_outputs.m_Q_thermal_ss); //[MWt]

	return;
}

void C_csp_nuclear_plant::startup(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	const C_csp_solver_sim_info &sim_info)
{
	// For now, define startup(...) shell that calls call() with operation mode defined.
	// Should eventually develop a startup method for the collector receiver

	// Set heliostat field call() parameters and solve
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::STARTUP;
	inputs.m_field_control = 1.0;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_nuclear_plant::on(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	double field_control,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{

	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::ON;
	inputs.m_field_control = field_control;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_nuclear_plant::estimates(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_est_out &est_out,
	const C_csp_solver_sim_info &sim_info)
{
	// For now, define estimates(...) shell that calls call() with operation mode defined.
	// Should eventually develop an estimate(...) method for the MSPT
	
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::ON;
	inputs.m_field_control = 1.0;

	C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver;

	call(weather, htf_state_in, inputs, cr_out_solver, sim_info);

	int mode = get_operating_state();

	if( mode == C_csp_collector_receiver::ON )
	{
		est_out.m_q_dot_avail = cr_out_solver.m_q_thermal;			//[MWt]
		est_out.m_m_dot_avail = cr_out_solver.m_m_dot_salt_tot;		//[kg/hr]
		est_out.m_T_htf_hot = cr_out_solver.m_T_salt_hot;			//[C]
		est_out.m_q_startup_avail = 0.0;
	}
	else
	{
		est_out.m_q_startup_avail = cr_out_solver.m_q_thermal;		//[MWt]
		est_out.m_q_dot_avail = 0.0;
		est_out.m_m_dot_avail = 0.0;
		est_out.m_T_htf_hot = 0.0;
	}
}

double C_csp_nuclear_plant::calculate_optical_efficiency( const C_csp_weatherreader::S_outputs &weather, const C_csp_solver_sim_info &sim )
{
    /*
    Evaluate optical efficiency. This is a required function for the parent class, 
    but doesn't do much other than simply call the optical efficiency model in this case.
    */
    return 1.0;
}

double C_csp_nuclear_plant::get_collector_area()
{
    //C_pt_heliostatfield::S_params *p = &mc_pt_heliostatfield.ms_params;

    //return p->m_dens_mirror * p->m_helio_height * p->m_helio_width * (double)p->m_helio_positions.nrows();

    //return mc_heat_input.m_A_sf;
    return mc_nuclear.m_dummy_area;
}

double C_csp_nuclear_plant::calculate_thermal_efficiency_approx( const C_csp_weatherreader::S_outputs &weather, double q_inc )
{

    return 1.0;

}


void C_csp_nuclear_plant::converged()
{
	mc_nuclear.converged();

    // Set reported heliostat converged value
    mc_reported_outputs.value(E_IS_FIELD_TRACKING_FINAL, 0);			//[-]

    // Set reported receiver converged value
    C_csp_collector_receiver::E_csp_cr_modes rec_op_mode_final;
    double rec_startup_time_remain_final, rec_startup_energy_remain_final;
    rec_startup_time_remain_final = rec_startup_energy_remain_final = std::numeric_limits<double>::quiet_NaN();
    mc_nuclear.get_converged_values(rec_op_mode_final,
        rec_startup_energy_remain_final, rec_startup_time_remain_final);
    mc_reported_outputs.value(E_REC_OP_MODE_FINAL, (int)rec_op_mode_final);
    mc_reported_outputs.value(E_REC_STARTUP_TIME_REMAIN_FINAL, rec_startup_time_remain_final);
    mc_reported_outputs.value(E_REC_STARTUP_ENERGY_REMAIN_FINAL, rec_startup_energy_remain_final);

	// Hardcode to test...
	//mc_reported_outputs.set_timestep_output(E_Q_DOT_THERMAL, mc_heat_input.ms_outputs.m_Q_thermal);	//[MWt]
	mc_reported_outputs.set_timestep_outputs();
}

void C_csp_nuclear_plant::write_output_intervals(double report_time_start,
	const std::vector<double> & v_temp_ts_time_end, double report_time_end)
{
	mc_reported_outputs.send_to_reporting_ts_array(report_time_start,
		v_temp_ts_time_end, report_time_end);
}
