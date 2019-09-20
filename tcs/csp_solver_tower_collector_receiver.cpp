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

#include "csp_solver_tower_collector_receiver.h"
#include "sam_csp_util.h"
#include <algorithm>

static C_csp_reported_outputs::S_output_info S_output_info[] =
{
	{C_csp_tower_collector_receiver::E_FIELD_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_FIELD_ETA_OPT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_FIELD_ADJUST, C_csp_reported_outputs::TS_WEIGHTED_AVE},

	{C_csp_tower_collector_receiver::E_Q_DOT_INC, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_ETA_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_Q_DOT_THERMAL, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_M_DOT_HTF, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_Q_DOT_STARTUP, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_T_HTF_IN, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_T_HTF_OUT, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_Q_DOT_PIPE_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
	{C_csp_tower_collector_receiver::E_Q_DOT_LOSS, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    // from transient model:
	{ C_csp_tower_collector_receiver::E_P_HEATTRACE, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_tower_collector_receiver::E_T_HTF_OUT_END, C_csp_reported_outputs::TS_LAST },
	{ C_csp_tower_collector_receiver::E_T_HTF_OUT_MAX, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_tower_collector_receiver::E_T_HTF_PANEL_OUT_MAX, C_csp_reported_outputs::TS_WEIGHTED_AVE },
	{ C_csp_tower_collector_receiver::E_T_WALL_INLET, C_csp_reported_outputs::TS_LAST },
	{ C_csp_tower_collector_receiver::E_T_WALL_OUTLET, C_csp_reported_outputs::TS_LAST },
	{ C_csp_tower_collector_receiver::E_T_RISER, C_csp_reported_outputs::TS_LAST },
	{ C_csp_tower_collector_receiver::E_T_DOWNC, C_csp_reported_outputs::TS_LAST },

    { C_csp_tower_collector_receiver::E_Q_DOT_INC1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_Q_DOT_INC2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_Q_DOT_INC3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HTF1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HTF2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HTF3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_IN1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_IN2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_IN3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_OUT1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_OUT2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HTF_OUT3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_Q_DOT_HX1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_Q_DOT_HX2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_Q_DOT_HX3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HX1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HX2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_M_DOT_HX3, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HX_OUT1, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HX_OUT2, C_csp_reported_outputs::TS_WEIGHTED_AVE},
    { C_csp_tower_collector_receiver::E_T_HX_OUT3, C_csp_reported_outputs::TS_WEIGHTED_AVE},

	csp_info_invalid	
};

C_csp_tower_collector_receiver::C_csp_tower_collector_receiver(std::vector<C_csp_mspt_collector_receiver> & collector_receivers):
    collector_receivers(collector_receivers)
{
	mc_reported_outputs.construct(S_output_info);
    hxs.assign(collector_receivers.size(), C_heat_exchanger());
}

C_csp_tower_collector_receiver::~C_csp_tower_collector_receiver()
{}

void C_csp_tower_collector_receiver::init(const C_csp_collector_receiver::S_csp_cr_init_inputs init_inputs, 
				C_csp_collector_receiver::S_csp_cr_solved_params & solved_params)
{
    // Declare instance of fluid class for FIELD fluid
    // Set fluid number and copy over fluid matrix if it makes sense
    if (m_field_fl != HTFProperties::User_defined && m_field_fl < HTFProperties::End_Library_Fluids)
    {
        if (!mc_field_htfProps.SetFluid(m_field_fl))
        {
            throw(C_csp_exception("Field HTF code is not recognized", "Tower Initialization"));
        }
    }
    else if (m_field_fl == HTFProperties::User_defined)
    {
        int n_rows = (int)m_field_fl_props.nrows();
        int n_cols = (int)m_field_fl_props.ncols();
        if (n_rows > 2 && n_cols == 7)
        {
            if (!mc_field_htfProps.SetUserDefinedFluid(m_field_fl_props))
            {
                error_msg = util::format(mc_field_htfProps.UserFluidErrMessage(), n_rows, n_cols);
                throw(C_csp_exception(error_msg, "Tower Initialization"));
            }
        }
        else
        {
            error_msg = util::format("The user defined field HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
            throw(C_csp_exception(error_msg, "Tower Initialization"));
        }
    }
    else
    {
        throw(C_csp_exception("Field HTF code is not recognized", "Tower Initialization"));
    }

    // Declare instance of fluid class for STORAGE fluid.
    // Set fluid number and copy over fluid matrix if it makes sense.
    if (m_tes_fl != HTFProperties::User_defined && m_tes_fl < HTFProperties::End_Library_Fluids)
    {
        if (!mc_store_htfProps.SetFluid(m_tes_fl))
        {
            throw(C_csp_exception("Storage HTF code is not recognized", "Tower Initialization"));
        }
    }
    else if (m_tes_fl == HTFProperties::User_defined)
    {
        int n_rows = (int)m_tes_fl_props.nrows();
        int n_cols = (int)m_tes_fl_props.ncols();
        if (n_rows > 2 && n_cols == 7)
        {
            if (!mc_store_htfProps.SetUserDefinedFluid(m_tes_fl_props))
            {
                error_msg = util::format(mc_store_htfProps.UserFluidErrMessage(), n_rows, n_cols);
                throw(C_csp_exception(error_msg, "Tower Initialization"));
            }
        }
        else
        {
            error_msg = util::format("The user defined storage HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
            throw(C_csp_exception(error_msg, "Tower Initialization"));
        }
    }
    else
    {
        throw(C_csp_exception("Storage HTF code is not recognized", "Tower Initialization"));
    }

    T_rec_hot_des += 273.15;    //[K] convert from C
    T_hx_cold_des += 273.15;    //[K] convert from C

    C_csp_collector_receiver::S_csp_cr_solved_params _solved_params;
    double q_dot_rec_des = 0.;
    double A_aper_total = 0.;
    double dP_sf = 0.;

    for (std::vector<int>::size_type i = 0; i != collector_receivers.size(); i++) {
        collector_receivers.at(i).init(init_inputs, _solved_params);
        
        if (i == 0) {
            solved_params.m_T_htf_cold_des = _solved_params.m_T_htf_cold_des;
            solved_params.m_P_cold_des = _solved_params.m_P_cold_des;
            solved_params.m_x_cold_des = _solved_params.m_x_cold_des;
        }
        q_dot_rec_des += _solved_params.m_q_dot_rec_des;
        A_aper_total += _solved_params.m_A_aper_total;
        dP_sf += _solved_params.m_dP_sf;
        
        hxs.at(i).init(mc_field_htfProps, mc_store_htfProps, hx_duty, m_dt_hot, T_rec_hot_des, T_hx_cold_des + m_dt_hot);
    }

    solved_params.m_T_htf_hot_des = _solved_params.m_T_htf_hot_des;     // of last receiver
    solved_params.m_q_dot_rec_des = q_dot_rec_des;
    solved_params.m_A_aper_total = A_aper_total;
    solved_params.m_dP_sf = dP_sf;

	return;
}

int C_csp_tower_collector_receiver::get_operating_state()
{
    //  Truth Table
    //  (for a two receiver tower or when one of the three receivers is off)
    //  +--------------+--------------+--------------+--------------+--------------+
    //  |              | OFF          | STARTUP      | ON           | STEADY_STATE |
    //  +--------------+--------------+--------------+--------------+--------------+
    //  | OFF          | OFF          | STARTUP      | ON           | STEADY_STATE |
    //  | STARTUP      | STARTUP      | STARTUP      | ON           | ON           |
    //  | ON           | ON           | ON           | ON           | ON           |
    //  | STEADY_STATE | STEADY_STATE | ON           | ON           | STEADY_STATE |
    //  +--------------+--------------+--------------+--------------+--------------+

    int operating_state = std::numeric_limits<int>::quiet_NaN();
    std::vector<int> operating_states;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        operating_states.push_back(it->get_operating_state());
    }

    std::vector<int>::iterator it;
    it = find(operating_states.begin(), operating_states.end(), C_csp_collector_receiver::ON);
    if (it != operating_states.end()) {     // one of the receivers is ON
        operating_state = C_csp_collector_receiver::ON;     // so tower is thus ON
    }
    else {
        it = find(operating_states.begin(), operating_states.end(), C_csp_collector_receiver::STEADY_STATE);
        if (it != operating_states.end()) {     // no receivers are ON and one of the receivers is at STEADY_STATE
            it = find(operating_states.begin(), operating_states.end(), C_csp_collector_receiver::STARTUP);
            if (it != operating_states.end()) { // no receivers are ON, one is at STEADY_STATE and one is at STARTUP
                operating_state = C_csp_collector_receiver::ON;
            }
            else {  // no receivers are ON, one is at STEADY_STATE and none are at STARTUP
                operating_state = C_csp_collector_receiver::STEADY_STATE;
            }
        }
        else {      // no receivers are ON and no receivers are at STEADY_STATE
            it = find(operating_states.begin(), operating_states.end(), C_csp_collector_receiver::STARTUP);
            if (it != operating_states.end()) {     // no receivers are ON, none are at STEADY_STATE and one of the receivers is at STARTUP
                operating_state = C_csp_collector_receiver::STARTUP;
            }
            else {  // no receivers are ON, none are at STEADY_STATE and none are at STARTUP
                operating_state = C_csp_collector_receiver::OFF;
            }
        }
    }

    return operating_state;
}

double C_csp_tower_collector_receiver::get_startup_time()
{
    // DO THEY STARTUP IN PARALLEL OR SERIES?

    double startup_time = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        startup_time = std::max(it->get_startup_time(), startup_time);         // assuming parallel
    }

    return startup_time;
}

double C_csp_tower_collector_receiver::get_startup_energy()     //MWh
{
    double startup_energy = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        startup_energy += it->get_startup_energy();
    }

    return startup_energy;
}

double C_csp_tower_collector_receiver::get_pumping_parasitic_coef()  //MWe/MWt
{
    // NOT ACCOUNTING FOR HXs
    
    double q_recs_des = 0.;
    double pumping_powers = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        double q_rec_des = it->get_design_thermal_power();
        pumping_powers += it->get_pumping_parasitic_coef() * q_rec_des;
        q_recs_des += q_rec_des;
    }

    return pumping_powers / q_recs_des;
}

double C_csp_tower_collector_receiver::get_min_power_delivery()    //MWt
{
    double min_power_delivery = std::numeric_limits<double>::infinity();

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        min_power_delivery = std::min(it->get_min_power_delivery(), min_power_delivery);
    }

    return min_power_delivery;
}


double C_csp_tower_collector_receiver::get_tracking_power()     //MWe
{
    double tracking_power = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        tracking_power += it->get_tracking_power();
    }

    return tracking_power;
}

double C_csp_tower_collector_receiver::get_col_startup_power()  //MWe-hr
{
    double startup_power = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        startup_power += it->get_col_startup_power();
    }

    return startup_power;
}

void C_csp_tower_collector_receiver::set_tes(C_csp_two_tank_two_hx_tes * tes)
{
    this->tes = tes;
}


void C_csp_tower_collector_receiver::call(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_collector_receiver::S_csp_cr_inputs &inputs,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
    cr_out_solver.m_q_startup = 0.;
    cr_out_solver.m_time_required_su = 0.;
    cr_out_solver.m_m_dot_salt_tot = 0.;
    cr_out_solver.m_q_thermal = 0.;
    cr_out_solver.m_component_defocus = 1.;
    cr_out_solver.m_is_recirculating = false;
    cr_out_solver.m_m_dot_store_tot = 0.;
    cr_out_solver.m_T_store_hot = 0.;
    cr_out_solver.m_E_fp_total = 0.;
    cr_out_solver.m_W_dot_col_tracking = 0.;
    cr_out_solver.m_W_dot_htf_pump = 0.;
    cr_out_solver.m_dP_sf = 0.;
    cr_out_solver.m_q_rec_heattrace = 0.;

    double q_dot_field_inc = 0.;
    double eta_weighted_sum = 0.;
    double sf_adjust_weighted_sum = 0.;
    double q_dot_rec_inc = 0.;
    double q_losses = 0.;
    double q_thermal = 0.;
    double m_dot_salt_tot = std::numeric_limits<double>::quiet_NaN();
    double q_dot_startup = 0.;
    double T_htf_in = std::numeric_limits<double>::quiet_NaN();
    double T_htf_out = std::numeric_limits<double>::quiet_NaN();
    double q_dot_piping_loss = 0.;
    double q_dot_loss = 0.;
    double q_heattrace = 0.;
    double T_htf_out_end = std::numeric_limits<double>::quiet_NaN();
    double T_htf_out_max = std::numeric_limits<double>::quiet_NaN();
    double T_htf_panel_out_max = std::numeric_limits<double>::quiet_NaN();
    double T_wall_inlet = std::numeric_limits<double>::quiet_NaN();
    double T_wall_outlet = std::numeric_limits<double>::quiet_NaN();
    double T_riser = std::numeric_limits<double>::quiet_NaN();
    double T_downc = std::numeric_limits<double>::quiet_NaN();
    double T_store_hot_weighted_sum = 0.;
   
    C_csp_solver_htf_1state htf_state_in_next = htf_state_in;
    C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver_prev;

    for (std::vector<int>::size_type i = 0; i != collector_receivers.size(); i++) {
        collector_receivers.at(i).call(weather, htf_state_in_next, inputs, cr_out_solver_prev, sim_info);

        cr_out_solver.m_q_startup += cr_out_solver_prev.m_q_startup;
        cr_out_solver.m_time_required_su = std::max(cr_out_solver.m_time_required_su, cr_out_solver_prev.m_time_required_su);
        cr_out_solver.m_m_dot_salt_tot = std::max(cr_out_solver.m_m_dot_salt_tot, cr_out_solver_prev.m_m_dot_salt_tot);
        cr_out_solver.m_q_thermal += cr_out_solver_prev.m_q_thermal;
        cr_out_solver.m_component_defocus = std::min(cr_out_solver.m_component_defocus, cr_out_solver_prev.m_component_defocus);  // NEED TO WEIGHT THIS
        cr_out_solver.m_is_recirculating = cr_out_solver.m_is_recirculating || cr_out_solver_prev.m_is_recirculating;  // if any are
        cr_out_solver.m_E_fp_total += cr_out_solver_prev.m_E_fp_total;
        cr_out_solver.m_W_dot_col_tracking += cr_out_solver_prev.m_W_dot_col_tracking;
        cr_out_solver.m_W_dot_htf_pump += cr_out_solver_prev.m_W_dot_htf_pump;
        cr_out_solver.m_dP_sf += cr_out_solver_prev.m_dP_sf;
        cr_out_solver.m_q_rec_heattrace += cr_out_solver_prev.m_q_rec_heattrace;

        double eff, T_cold_rec_K, T_hot_tes_K, q_trans, m_dot_tes;
        hxs.at(i).hx_charge_mdot_field(cr_out_solver_prev.m_T_salt_hot + 273.15, cr_out_solver_prev.m_m_dot_salt_tot / 3600., tes->get_cold_temp(),
            eff, T_cold_rec_K, T_hot_tes_K, q_trans, m_dot_tes);

        cr_out_solver.m_m_dot_store_tot += m_dot_tes * 3600.;
        T_store_hot_weighted_sum += T_hot_tes_K * m_dot_tes * 3600.;

        // Update tanks -> ASSUMING THE COLD TANK CAN TAKE ALL THE CHARGE 
        tes->ms_params.m_is_hx = false;     // charging from the receivers is direct storage
        double T_cold_tes_K;
        C_csp_tes::S_csp_tes_outputs tes_outputs;
        tes->charge(sim_info.ms_ts.m_step, weather.m_tdry + 273.15, m_dot_tes, T_hot_tes_K, T_cold_tes_K, tes_outputs);
        tes->ms_params.m_is_hx = true;

        cr_out_solver.m_T_salt_hot = T_cold_rec_K - 273.15;  // of last receiver

        htf_state_in_next.m_temp = T_cold_rec_K - 273.15;
        htf_state_in_next.m_m_dot = cr_out_solver_prev.m_m_dot_salt_tot / 3600.;    // THIS WILL BE CONTROLLED BY A MEQ VIA VALVING

        C_pt_sf_perf_interp::S_outputs field_outputs = collector_receivers.at(i).get_field_outputs();
        C_pt_receiver::S_outputs receiver_outputs = collector_receivers.at(i).get_receiver_outputs();
        q_dot_field_inc += field_outputs.m_q_dot_field_inc;
        eta_weighted_sum += field_outputs.m_eta_field * collector_receivers.at(i).get_collector_area();
        sf_adjust_weighted_sum += field_outputs.m_sf_adjust_out * collector_receivers.at(i).get_collector_area();
        q_dot_rec_inc += receiver_outputs.m_q_dot_rec_inc;
        q_losses += (1 - receiver_outputs.m_eta_therm) * receiver_outputs.m_q_dot_rec_inc;
        q_thermal += receiver_outputs.m_Q_thermal;
        m_dot_salt_tot = std::max(receiver_outputs.m_m_dot_salt_tot, m_dot_salt_tot);
        q_dot_startup += receiver_outputs.m_q_startup / (receiver_outputs.m_time_required_su / 3600.);
        T_htf_in = htf_state_in.m_temp;
        T_htf_out = T_cold_rec_K - 273.15;              // of last receiver
        q_dot_piping_loss += receiver_outputs.m_q_dot_piping_loss;
        q_dot_loss += receiver_outputs.m_q_rad_sum + receiver_outputs.m_q_conv_sum;
        q_heattrace += receiver_outputs.m_q_heattrace / (receiver_outputs.m_time_required_su / 3600.);
        T_htf_out_end = T_cold_rec_K - 273.15;     // of last receiver
        T_htf_out_max = T_cold_rec_K - 273.15;      // of last receiver
        T_htf_panel_out_max = std::max(receiver_outputs.m_max_rec_tout, T_htf_panel_out_max);
        if (i == 0) T_wall_inlet = receiver_outputs.m_Twall_inlet;       // of first receiver
        T_wall_outlet = receiver_outputs.m_Twall_outlet;        // of last receiver
        if (i == 0) T_riser = receiver_outputs.m_Triser;                 // of first receiver
        T_downc = T_cold_rec_K - 273.15;                    // of last receiver

        if (i == 0) {
            mc_reported_outputs.value(E_Q_DOT_INC1, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF1, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN1, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT1, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX1, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX1, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT1, T_hot_tes_K - 273.15);                           //[C]
        }
        else if (i == 1) {
            mc_reported_outputs.value(E_Q_DOT_INC2, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF2, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN2, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT2, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX2, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX2, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT2, T_hot_tes_K - 273.15);                           //[C]
        }
        else if (i == 2) {
            mc_reported_outputs.value(E_Q_DOT_INC3, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF3, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN3, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT3, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX3, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX3, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT3, T_hot_tes_K - 273.15);                           //[C]
        }
    }

    cr_out_solver.m_T_store_hot = T_store_hot_weighted_sum / cr_out_solver.m_m_dot_store_tot - 273.15; //[C]

	mc_reported_outputs.value(E_FIELD_Q_DOT_INC, q_dot_field_inc);	                                //[MWt]
    double collector_areas = get_collector_area();
	mc_reported_outputs.value(E_FIELD_ETA_OPT, eta_weighted_sum / collector_areas);			        //[-]
	mc_reported_outputs.value(E_FIELD_ADJUST, sf_adjust_weighted_sum / collector_areas);			//[-]
	mc_reported_outputs.value(E_Q_DOT_INC, q_dot_rec_inc);	                                        //[MWt]
	mc_reported_outputs.value(E_ETA_THERMAL, std::max(1. - q_losses / q_dot_rec_inc, 0.));		    //[-]
	mc_reported_outputs.value(E_Q_DOT_THERMAL, q_thermal);	                                        //[MWt]
	mc_reported_outputs.value(E_M_DOT_HTF, m_dot_salt_tot);	                                        //[kg/hr]
	// If startup, then timestep may have changed (why not report this from 222 in MWt?)
	mc_reported_outputs.value(E_Q_DOT_STARTUP, q_dot_startup);		                                //[MWt])
	mc_reported_outputs.value(E_T_HTF_IN, T_htf_in);									            //[C]
	mc_reported_outputs.value(E_T_HTF_OUT, T_htf_out);		                                        //[C]
	mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, q_dot_piping_loss);	                            //[MWt]
    mc_reported_outputs.value(E_Q_DOT_LOSS, q_dot_loss); //MWt
    // from transient model:
	mc_reported_outputs.value(E_P_HEATTRACE, q_heattrace);		                                    //[MWt])
	mc_reported_outputs.value(E_T_HTF_OUT_END, T_htf_out_end);	                                    //[C]
	mc_reported_outputs.value(E_T_HTF_OUT_MAX, T_htf_out_max);	                                    //[C]
	mc_reported_outputs.value(E_T_HTF_PANEL_OUT_MAX, T_htf_panel_out_max);	                        //[C]
	mc_reported_outputs.value(E_T_WALL_INLET, T_wall_inlet);	                                    //[C]
	mc_reported_outputs.value(E_T_WALL_OUTLET, T_wall_outlet);	                                    //[C]
	mc_reported_outputs.value(E_T_RISER, T_riser);	                                                //[C]
	mc_reported_outputs.value(E_T_DOWNC, T_downc);	                                                //[C]

    return;
}

void C_csp_tower_collector_receiver::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
	const C_csp_solver_sim_info &sim_info)
{
    cr_out_solver.m_q_startup = 0.;
    cr_out_solver.m_time_required_su = 0.;
    cr_out_solver.m_m_dot_salt_tot = 0.;
    cr_out_solver.m_q_thermal = 0.;
    cr_out_solver.m_component_defocus = 1.;
    cr_out_solver.m_is_recirculating = false;
    cr_out_solver.m_m_dot_store_tot = 0.;
    cr_out_solver.m_T_store_hot = 0.;
    cr_out_solver.m_E_fp_total = 0.;
    cr_out_solver.m_W_dot_col_tracking = 0.;
    cr_out_solver.m_W_dot_htf_pump = 0.;
    cr_out_solver.m_dP_sf = 0.;
    cr_out_solver.m_q_rec_heattrace = 0.;

    double q_dot_field_inc = 0.;
    double eta_weighted_sum = 0.;
    double sf_adjust_weighted_sum = 0.;
    double q_dot_rec_inc = 0.;
    double q_losses = 0.;
    double q_thermal = 0.;
    double m_dot_salt_tot = std::numeric_limits<double>::quiet_NaN();
    double q_dot_startup = 0.;
    double T_htf_in = std::numeric_limits<double>::quiet_NaN();
    double T_htf_out = std::numeric_limits<double>::quiet_NaN();
    double q_dot_piping_loss = 0.;
    double q_dot_loss = 0.;
    double q_heattrace = 0.;
    double T_htf_out_end = std::numeric_limits<double>::quiet_NaN();
    double T_htf_out_max = std::numeric_limits<double>::quiet_NaN();
    double T_htf_panel_out_max = std::numeric_limits<double>::quiet_NaN();
    double T_wall_inlet = std::numeric_limits<double>::quiet_NaN();
    double T_wall_outlet = std::numeric_limits<double>::quiet_NaN();
    double T_riser = std::numeric_limits<double>::quiet_NaN();
    double T_downc = std::numeric_limits<double>::quiet_NaN();
    double T_store_hot_weighted_sum = 0.;

    C_csp_solver_htf_1state htf_state_in_next = htf_state_in;
    C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver_prev;

    
    for (std::vector<int>::size_type i = 0; i != collector_receivers.size(); i++) {
        // SHOULD THIS BE IN REVERSE ORDER?
        collector_receivers.at(i).off(weather, htf_state_in_next, cr_out_solver_prev, sim_info);

        cr_out_solver.m_q_startup += cr_out_solver_prev.m_q_startup;
        cr_out_solver.m_time_required_su = std::max(cr_out_solver.m_time_required_su, cr_out_solver_prev.m_time_required_su);
        cr_out_solver.m_m_dot_salt_tot = std::max(cr_out_solver.m_m_dot_salt_tot, cr_out_solver_prev.m_m_dot_salt_tot);
        cr_out_solver.m_q_thermal += cr_out_solver_prev.m_q_thermal;
        cr_out_solver.m_component_defocus = std::min(cr_out_solver.m_component_defocus, cr_out_solver_prev.m_component_defocus);  // NEED TO WEIGHT THIS
        cr_out_solver.m_is_recirculating = cr_out_solver.m_is_recirculating || cr_out_solver_prev.m_is_recirculating;  // if any are
        cr_out_solver.m_E_fp_total += cr_out_solver_prev.m_E_fp_total;
        cr_out_solver.m_W_dot_col_tracking += cr_out_solver_prev.m_W_dot_col_tracking;
        cr_out_solver.m_W_dot_htf_pump += cr_out_solver_prev.m_W_dot_htf_pump;
        cr_out_solver.m_dP_sf += cr_out_solver_prev.m_dP_sf;
        cr_out_solver.m_q_rec_heattrace += cr_out_solver_prev.m_q_rec_heattrace;

        double eff, T_cold_rec_K, T_hot_tes_K, q_trans, m_dot_tes;
        hxs.at(i).hx_charge_mdot_field(cr_out_solver_prev.m_T_salt_hot + 273.15, cr_out_solver_prev.m_m_dot_salt_tot / 3600., tes->get_cold_temp(),
            eff, T_cold_rec_K, T_hot_tes_K, q_trans, m_dot_tes);

        cr_out_solver.m_m_dot_store_tot += m_dot_tes * 3600.;
        T_store_hot_weighted_sum += T_hot_tes_K * m_dot_tes * 3600.;

        // Update tanks
        tes->ms_params.m_is_hx = false;     // charging from the receivers is direct storage
        C_csp_tes::S_csp_tes_outputs tes_outputs;
        tes->idle(sim_info.ms_ts.m_step, weather.m_tdry + 273.15, tes_outputs);
        tes->ms_params.m_is_hx = true;

        cr_out_solver.m_T_salt_hot = T_cold_rec_K - 273.15;  // from last receiver

        htf_state_in_next.m_temp = T_cold_rec_K - 273.15;
        htf_state_in_next.m_m_dot = cr_out_solver_prev.m_m_dot_salt_tot / 3600.;    // THIS WILL BE CONTROLLED BY A MEQ VIA VALVING

        C_pt_sf_perf_interp::S_outputs field_outputs = collector_receivers.at(i).get_field_outputs();
        C_pt_receiver::S_outputs receiver_outputs = collector_receivers.at(i).get_receiver_outputs();
        q_dot_field_inc += field_outputs.m_q_dot_field_inc;
        eta_weighted_sum += field_outputs.m_eta_field * collector_receivers.at(i).get_collector_area();
        sf_adjust_weighted_sum += field_outputs.m_sf_adjust_out * collector_receivers.at(i).get_collector_area();
        q_dot_rec_inc += receiver_outputs.m_q_dot_rec_inc;
        q_losses += (1 - receiver_outputs.m_eta_therm) * receiver_outputs.m_q_dot_rec_inc;
        q_thermal += receiver_outputs.m_Q_thermal;
        m_dot_salt_tot = std::max(receiver_outputs.m_m_dot_salt_tot, m_dot_salt_tot);
        q_dot_startup += receiver_outputs.m_q_startup / (receiver_outputs.m_time_required_su / 3600.);
        T_htf_in = htf_state_in.m_temp;
        T_htf_out = T_cold_rec_K - 273.15;              // of last receiver
        q_dot_piping_loss += receiver_outputs.m_q_dot_piping_loss;
        q_dot_loss += receiver_outputs.m_q_rad_sum + receiver_outputs.m_q_conv_sum;
        q_heattrace += receiver_outputs.m_q_heattrace / (receiver_outputs.m_time_required_su / 3600.);
        T_htf_out_end = T_cold_rec_K - 273.15;     // of last receiver
        T_htf_out_max = T_cold_rec_K - 273.15;      // of last receiver
        T_htf_panel_out_max = std::max(receiver_outputs.m_max_rec_tout, T_htf_panel_out_max);
        if (i == 0) T_wall_inlet = receiver_outputs.m_Twall_inlet;       // of first receiver
        T_wall_outlet = receiver_outputs.m_Twall_outlet;        // of last receiver
        if (i == 0) T_riser = receiver_outputs.m_Triser;                 // of first receiver
        T_downc = T_cold_rec_K - 273.15;                    // of last receiver

        if (i == 0) {
            mc_reported_outputs.value(E_Q_DOT_INC1, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF1, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN1, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT1, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX1, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX1, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT1, T_hot_tes_K - 273.15);                           //[C]
        }
        else if (i == 1) {
            mc_reported_outputs.value(E_Q_DOT_INC2, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF2, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN2, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT2, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX2, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX2, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT2, T_hot_tes_K - 273.15);                           //[C]
        }
        else if (i == 2) {
            mc_reported_outputs.value(E_Q_DOT_INC3, receiver_outputs.m_q_dot_rec_inc);	            //[MWt]
            mc_reported_outputs.value(E_M_DOT_HTF3, receiver_outputs.m_m_dot_salt_tot);	            //[kg/hr]
            mc_reported_outputs.value(E_T_HTF_IN3, receiver_outputs.m_T_salt_cold);					//[C]
            mc_reported_outputs.value(E_T_HTF_OUT3, receiver_outputs.m_T_salt_hot);		            //[C]
            mc_reported_outputs.value(E_Q_DOT_HX3, q_trans);		                                //[MWt]
            mc_reported_outputs.value(E_M_DOT_HX3, m_dot_tes * 3600.);                              //[kg/hr]
            mc_reported_outputs.value(E_T_HX_OUT3, T_hot_tes_K - 273.15);                           //[C]
        }
    }

    cr_out_solver.m_T_store_hot = T_store_hot_weighted_sum / cr_out_solver.m_m_dot_store_tot - 273.15; //[C]

    mc_reported_outputs.value(E_FIELD_Q_DOT_INC, q_dot_field_inc);	                                //[MWt]
    double collector_areas = get_collector_area();
    mc_reported_outputs.value(E_FIELD_ETA_OPT, eta_weighted_sum / collector_areas);			        //[-]
    mc_reported_outputs.value(E_FIELD_ADJUST, sf_adjust_weighted_sum / collector_areas);			//[-]
    mc_reported_outputs.value(E_Q_DOT_INC, q_dot_rec_inc);	                                        //[MWt]
    mc_reported_outputs.value(E_ETA_THERMAL, std::max(1. - q_losses / q_dot_rec_inc, 0.));		    //[-]
    mc_reported_outputs.value(E_Q_DOT_THERMAL, q_thermal);	                                        //[MWt]
    mc_reported_outputs.value(E_M_DOT_HTF, m_dot_salt_tot);	                                        //[kg/hr]
    // If startup, then timestep may have changed (why not report this from 222 in MWt?)
    mc_reported_outputs.value(E_Q_DOT_STARTUP, q_dot_startup);		                                //[MWt])
    mc_reported_outputs.value(E_T_HTF_IN, T_htf_in);									            //[C]
    mc_reported_outputs.value(E_T_HTF_OUT, T_htf_out);		                                        //[C]
    mc_reported_outputs.value(E_Q_DOT_PIPE_LOSS, q_dot_piping_loss);	                            //[MWt]
    mc_reported_outputs.value(E_Q_DOT_LOSS, q_dot_loss); //MWt
    // from transient model:
    mc_reported_outputs.value(E_P_HEATTRACE, q_heattrace);		                                    //[MWt])
    mc_reported_outputs.value(E_T_HTF_OUT_END, T_htf_out_end);	                                    //[C]
    mc_reported_outputs.value(E_T_HTF_OUT_MAX, T_htf_out_max);	                                    //[C]
    mc_reported_outputs.value(E_T_HTF_PANEL_OUT_MAX, T_htf_panel_out_max);	                        //[C]
    mc_reported_outputs.value(E_T_WALL_INLET, T_wall_inlet);	                                    //[C]
    mc_reported_outputs.value(E_T_WALL_OUTLET, T_wall_outlet);	                                    //[C]
    mc_reported_outputs.value(E_T_RISER, T_riser);	                                                //[C]
    mc_reported_outputs.value(E_T_DOWNC, T_downc);	                                                //[C]

    return;
}

void C_csp_tower_collector_receiver::startup(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	const C_csp_solver_sim_info &sim_info)
{
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::STARTUP;
	inputs.m_field_control = 1.0;

    call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_tower_collector_receiver::on(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	double field_control,
	C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
	const C_csp_solver_sim_info &sim_info)
{
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::ON;
	inputs.m_field_control = field_control;

    call(weather, htf_state_in, inputs, cr_out_solver, sim_info);
}

void C_csp_tower_collector_receiver::estimates(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	C_csp_collector_receiver::S_csp_cr_est_out &est_out,
	const C_csp_solver_sim_info &sim_info)
{	
	C_csp_collector_receiver::S_csp_cr_inputs inputs;
	inputs.m_input_operation_mode = C_csp_collector_receiver::STEADY_STATE;
	inputs.m_field_control = 1.0;

    C_csp_collector_receiver::S_csp_cr_out_solver cr_out_solver;

    call(weather, htf_state_in, inputs, cr_out_solver, sim_info);

   	int mode = get_operating_state();

	if( mode == C_csp_collector_receiver::ON ) {
		est_out.m_q_startup_avail = 0.;
		est_out.m_q_dot_avail = cr_out_solver.m_q_thermal;			    //[MWt]
		est_out.m_m_dot_avail = cr_out_solver.m_m_dot_salt_tot;         //[kg/hr]
		est_out.m_T_htf_hot = cr_out_solver.m_T_salt_hot;		    	//[C], last receiver
	}
	else {
		est_out.m_q_startup_avail = cr_out_solver.m_q_thermal;  		//[MWt]
		est_out.m_q_dot_avail = 0.;
		est_out.m_m_dot_avail = 0.;
		est_out.m_T_htf_hot = std::numeric_limits<double>::quiet_NaN();
	}
}

double C_csp_tower_collector_receiver::calculate_optical_efficiency( const C_csp_weatherreader::S_outputs &weather, const C_csp_solver_sim_info &sim )
{
    // simple weighted average of the different collector_receiver optical efficiencies

    double collector_areas = 0.;
    double eff_weighted_sum = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        double collector_area = it->get_collector_area();
        eff_weighted_sum += it->calculate_optical_efficiency(weather, sim) * collector_area;
        collector_areas += collector_area;
    }

    return eff_weighted_sum / collector_areas;
}

double C_csp_tower_collector_receiver::get_collector_area()
{
    double A_sfs = 0.;
    
    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        A_sfs += it->get_collector_area();
    }

    return A_sfs;
}

double C_csp_tower_collector_receiver::calculate_thermal_efficiency_approx( const C_csp_weatherreader::S_outputs &weather, double q_inc )
{
    /* 
    A very approximate thermal efficiency used for quick optimization performance projections
    */

    double q_inc_rec = q_inc / collector_receivers.size();      // assume equal amount of flux on the different receivers
    double q_losses = 0.;

    for (std::vector<C_csp_mspt_collector_receiver>::iterator it = collector_receivers.begin(); it != collector_receivers.end(); ++it) {
        double eta_therm = it->calculate_thermal_efficiency_approx(weather, q_inc_rec);
        q_losses += (1 - eta_therm) * q_inc_rec;
    }

    return std::max(1. - q_losses / q_inc, 0.);
}


void C_csp_tower_collector_receiver::converged()
{
    for (std::vector<int>::size_type i = 0; i != collector_receivers.size(); i++) {
        collector_receivers.at(i).converged();
        hxs.at(i).converged();
    }

	mc_reported_outputs.set_timestep_outputs();
}

void C_csp_tower_collector_receiver::write_output_intervals(double report_time_start,
	const std::vector<double> & v_temp_ts_time_end, double report_time_end)
{
	mc_reported_outputs.send_to_reporting_ts_array(report_time_start,
		v_temp_ts_time_end, report_time_end);
}