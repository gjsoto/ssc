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

#ifndef __csp_solver_nuclear_
#define __csp_solver_nuclear_

#include "csp_solver_util.h"

#include "htf_props.h"
#include "ngcc_powerblock.h"
#include "csp_solver_core.h"
#include "csp_solver_pt_receiver.h"

class C_nuclear : public C_pt_receiver
{
// The steady-state receiver (as opposed to the transient, for example)

public:
	// Class to save messages for up stream classes
	C_csp_messages csp_messages;
    
    int m_ncall;
    int m_itermode;

	// Data
    double m_dummy_area;
    double m_q_dot_nuc_res;
    double m_T_salt_hot_target;	    //[C], convert to K in init() call
    double m_m_dot_htf_max;			//[kg/s];		
    double m_od_control;	
    
	int m_n_panels;					//[-]
	double m_d_rec;					//[m]
	double m_h_rec;					//[m]

    double m_A_sf;                      //[m2] "solar field" area, should delete this later
    
    double m_nuclear_su_delay;      //[hr] required startup time
    double m_nuclear_qf_delay;      //[-] required startup energy as fraction of design thermal output
    
    struct s_steady_state_soln
	{
		C_csp_collector_receiver::E_csp_cr_modes mode;
		bool nuc_is_off;
		int itermode;

		double hour;				// Hour of the year 
		double T_amb;				// Dry bulb temperature (K)
		double T_dp;				// Dewpoint temperature (K)
		double v_wind_10;			// Wind speed at 10m (m/s)
		double p_amb;				// Ambient pressure (Pa)

		double dni;					// DNI for this solution
		double od_control;          // Defocus control

		
		double m_dot_salt_tot;      // Total salt mass flow (kg/s)
		double T_salt_cold_in;		// Cold salt inlet temperature (K)
		double T_salt_hot;			// Salt @ nuclear hot outlet T including piping loss (K)
        double T_salt_hot_rec;      // Salt @ nuclear hot outlet T before piping loss (K)
		double T_salt_props;		// Temperature at which salt properties are evaluated

		double u_salt;				// Salt velocity (m/s)
		double f;					// Friction factor

		double Q_inc_sum;			// Total absorbed solar energy (W)
		double Q_conv_sum;			// Total convection loss (W)
		double Q_rad_sum;			// Total radiation loss (W)
		double Q_abs_sum;			// Total energy transferred to HTF, not including piping loss (W)
		double Q_dot_piping_loss;   // Piping loss (W)
		double Q_inc_min;			// Minimum absorbed solar energy on any panel (W)
		double Q_thermal;			// Thermal power delivered to fluid (less piping loss) (W)

		double eta_therm;			// Receiver thermal efficiency (energy to HTF not including piping loss / Absorbed solar energy)
        double delta_T_piping;      // Temperature change from thermal loss in piping (K)

        std::vector<double> m_dot_salt_path;	// Salt mass flow per path (kg/s)
        std::vector<double> T_salt_hot_rec_path;      // Receiver flow path outlet T before piping loss (K)
        std::vector<double> Q_abs_path; // Total energy transferred to HTF per path, not including piping loss (W)

		util::matrix_t<double> T_s;			// Average external tube T (K)
		util::matrix_t<double> T_panel_out; // Panel HTF outlet T (K)
		util::matrix_t<double> T_panel_in;	// Panel HTF inlet T (K)
		util::matrix_t<double> T_panel_ave; // Panel average HTF T (k)

		util::matrix_t<double> q_dot_inc;  // Panel absorbed solar energy (W)
		util::matrix_t<double> q_dot_conv; // Panel convection loss (W)
		util::matrix_t<double> q_dot_rad;  // Panel radiation loss (W)
		util::matrix_t<double> q_dot_loss; // Panel convection + radiation loss (W)
		util::matrix_t<double> q_dot_abs;  // Panel energy to HTF (W)

		s_steady_state_soln()
		{
			clear();
		}

		void clear()
		{
			hour = T_amb = T_dp = v_wind_10 = p_amb = std::numeric_limits<double>::quiet_NaN();
			dni = od_control  = m_dot_salt_tot = T_salt_cold_in = T_salt_hot = T_salt_hot_rec = T_salt_props = std::numeric_limits<double>::quiet_NaN();
			u_salt = f = Q_inc_sum = Q_conv_sum = Q_rad_sum = Q_abs_sum = Q_dot_piping_loss = Q_inc_min = Q_thermal = eta_therm = delta_T_piping = std::numeric_limits<double>::quiet_NaN();

            mode = C_csp_collector_receiver::E_csp_cr_modes::ON;
            itermode = -1;
			nuc_is_off = false;
            m_dot_salt_path.clear();
            T_salt_hot_rec_path.clear();
            Q_abs_path.clear();
		}

	};

	s_steady_state_soln m_mflow_soln_prev;  // Steady state solution using actual DNI from the last call to the model
    bool use_previous_solution(const s_steady_state_soln& soln, const s_steady_state_soln& soln_prev);

	C_csp_collector_receiver::E_csp_cr_modes m_mode_initial;
    double m_E_su_init;             //[W-hr] Initial startup energy
    double m_t_su_init;             //[hr] Startup time requirement
	
	S_outputs outputs;

	// Methods
	C_nuclear();

	~C_nuclear(){};

	virtual void init();
    
    virtual int get_operating_state();

	virtual void call(const C_csp_weatherreader::S_outputs &weather, 
		const C_csp_solver_htf_1state &htf_state_in, 
		const C_pt_receiver::S_inputs &inputs,
		const C_csp_solver_sim_info &sim_info);

	virtual void off(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		const C_csp_solver_sim_info &sim_info);

	virtual void converged();

    virtual double get_remaining_startup_energy();

    void calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc, double &WdotPump_calc);

    virtual double get_pumping_parasitic_coef();

    virtual double area_proj();
    
    virtual double get_startup_time();
    
    virtual double get_startup_energy();

};

#endif // __csp_solver_nuclear_
