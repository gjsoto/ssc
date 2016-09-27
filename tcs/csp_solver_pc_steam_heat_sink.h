#ifndef __csp_solver_pc_steam_heat_sink_
#define __csp_solver_pc_steam_heat_sink_

#include "csp_solver_core.h"
#include "csp_solver_util.h"

#include "water_properties.h"

class C_pc_steam_heat_sink : public C_csp_power_cycle
{
private:

	double m_max_frac;		//[-]

	void check_double_params_are_set();

public:
	// Class to save messages for up stream classes
	C_csp_messages mc_csp_messages;

	water_state mc_water_props;

	struct S_params
	{
		double m_x_hot_des;		//[-]
		double m_T_hot_des;		//[C]
		double m_P_hot_des;		//[kPa]
		double m_T_cold_des;	//[C]
		double m_dP_frac_des;	//[-]
		double m_q_dot_des;		//[MWt]
		double m_m_dot_max_frac;	//[-]
	
		double m_pump_eta_isen;	//[-] Isentropic efficiency of pump

		S_params()
		{
			m_T_hot_des = m_P_hot_des = m_T_cold_des = 
				m_dP_frac_des = m_q_dot_des = m_m_dot_max_frac =
				
				m_pump_eta_isen = std::numeric_limits<double>::quiet_NaN();
		}
	};

	S_params ms_params;

	C_pc_steam_heat_sink()
	{
		m_max_frac = 100.0;

		m_is_sensible_htf = false;	//[-] STEAM
	}

	~C_pc_steam_heat_sink(){};

	virtual void init(C_csp_power_cycle::S_solved_params &solved_params);

	virtual int get_operating_state();

	virtual double get_cold_startup_time();
	virtual double get_warm_startup_time();
	virtual double get_hot_startup_time();
	virtual double get_standby_energy_requirement();    //[MW]
	virtual double get_cold_startup_energy();    //[MWh]
	virtual double get_warm_startup_energy();    //[MWh]
	virtual double get_hot_startup_energy();    //[MWh]
	virtual double get_max_thermal_power();     //MW
	virtual double get_min_thermal_power();     //MW
	virtual double get_efficiency_at_TPH(double T_degC, double P_atm, double relhum_pct);
	virtual double get_efficiency_at_load(double load_frac);

	// This can vary between timesteps for Type224, depending on remaining startup energy and time
	virtual double get_max_q_pc_startup();		//[MWt]

	virtual void call(const C_csp_weatherreader::S_outputs &weather,
		C_csp_solver_htf_1state &htf_state_in,
		const C_csp_power_cycle::S_control_inputs &inputs,
		C_csp_power_cycle::S_csp_pc_out_solver &out_solver,
		//C_csp_power_cycle::S_csp_pc_out_report &out_report,
		const C_csp_solver_sim_info &sim_info);

	virtual void converged();

	virtual void write_output_intervals(double report_time_start,
		const std::vector<double> & v_temp_ts_time_end, double report_time_end);

	virtual void assign(int index, float *p_reporting_ts_array, int n_reporting_ts_array);
};

#endif // __csp_solver_pc_steam_heat_sink_