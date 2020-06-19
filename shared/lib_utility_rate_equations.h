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

#ifndef _LIB_UTILITY_RATE_EQUATIONS_H_
#define _LIB_UTILITY_RATE_EQUATIONS_H_

#include <vector>

#include "../ssc/core.h"
#include "lib_util.h"

class ur_month
{
public:
	// period numbers
	std::vector<int> ec_periods;
	std::vector<int> dc_periods;
	// track period numbers at 12a, 6a, 12p and 6p for rollover applications. Weekdays only considered
	std::vector<int> ec_rollover_periods;
	// monthly values
	// net energy use per month
	ssc_number_t energy_net;
	// hours per period per month
	int hours_per_month;
	// energy use period and tier
	util::matrix_t<ssc_number_t> ec_energy_use;
	// handle changing period tiers on monthly basis if kWh/kW
	std::vector<std::vector<int> >  ec_periods_tiers; // tier numbers
	// energy surplus - extra generated by system that is either sold or curtailed.
	util::matrix_t<ssc_number_t> ec_energy_surplus;
	// peak demand per period
	std::vector<ssc_number_t> dc_tou_peak;
	std::vector<int> dc_tou_peak_hour;
	ssc_number_t dc_flat_peak;
	int dc_flat_peak_hour;
	// energy tou charges
	util::matrix_t<ssc_number_t>  ec_tou_ub_init;
	util::matrix_t<ssc_number_t>  ec_tou_br_init;
	util::matrix_t<ssc_number_t>  ec_tou_sr_init;
	// may change based on units and year
	util::matrix_t<ssc_number_t>  ec_tou_ub;
	util::matrix_t<ssc_number_t>  ec_tou_br;
	util::matrix_t<ssc_number_t>  ec_tou_sr;
	util::matrix_t<int>  ec_tou_units;
	// calculated charges per period and tier
	util::matrix_t<ssc_number_t>  ec_charge;
	// demand tou charges
	util::matrix_t<ssc_number_t>  dc_tou_ub;
	util::matrix_t<ssc_number_t>  dc_tou_ch;
	// demand flat charges
	std::vector<ssc_number_t>  dc_flat_ub;
	std::vector<ssc_number_t>  dc_flat_ch;
	// calculated charges per period
	std::vector<double>  dc_tou_charge;
	ssc_number_t dc_flat_charge;

	ur_month();
	ur_month(const ur_month& tmp);

	// Runs each step
	void update_net_and_peak(double energy, double power, int step);
};

class rate_data {

public:
	// schedule outputs
	std::vector<int> m_ec_tou_sched;
	std::vector<int> m_dc_tou_sched;
	std::vector<ur_month> m_month;
	std::vector<int> m_ec_periods; // period number

	// time step rates
	std::vector<ssc_number_t> m_ec_ts_sell_rate;
	std::vector<ssc_number_t> m_ec_ts_buy_rate;

	// track initial values - may change based on units
	std::vector<std::vector<int> >  m_ec_periods_tiers_init; // tier numbers
	std::vector<int> m_dc_tou_periods; // period number
	std::vector<std::vector<int> >  m_dc_tou_periods_tiers; // tier numbers
	std::vector<std::vector<int> >  m_dc_flat_tiers; // tier numbers for each month of flat demand charge
	size_t m_num_rec_yearly;

	std::vector<ssc_number_t> rate_scale;
	std::vector<ssc_number_t> dc_hourly_peak;
	std::vector<ssc_number_t> monthly_dc_fixed;
	std::vector<ssc_number_t> monthly_dc_tou;

	bool tou_demand_single_peak;

    bool enable_nm; // 0 or 1
    bool nm_credits_w_rollover; // rate option 0 only
    int net_metering_credit_month;
    double nm_credit_sell_rate;

	rate_data();
	rate_data(const rate_data& tmp);

	void init_energy_rates(bool gen_only);

	void init(int num_rec_yearly);
	void setup_time_series(size_t cnt, ssc_number_t* ts_sr, ssc_number_t* ts_br);
	void setup_energy_rates(ssc_number_t* ec_weekday, ssc_number_t* ec_weekend, size_t ec_tou_rows, ssc_number_t* ec_tou_in, bool sell_eq_buy);
	void setup_demand_charges(ssc_number_t* dc_weekday, ssc_number_t* dc_weekend, size_t dc_tou_rows, ssc_number_t* dc_tou_in, size_t dc_flat_rows, ssc_number_t* dc_flat_in);

	// Runs each step
	void sort_energy_to_periods(int month, double energy, int step); // Net metering only
	void find_dc_tou_peak(int month, double power, int step);
	int get_tou_row(int year_one_index, int month);

	// Runs each month
	void init_dc_peak_vectors(int month);
	ssc_number_t get_demand_charge(int month, int year); 
    // Returns error codes so compute module can print errors. 0: no error, 10x: error in previous month, 20x: error in current month. x is the period where the error occured
    int transfer_surplus(ur_month& curr_month, ur_month& prev_month);
    void compute_surplus(ur_month& curr_month);
};

#endif // _LIB_UTILITY_RATE_EQUATIONS_H_
