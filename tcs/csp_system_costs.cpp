#include "csp_system_costs.h"
#include "csp_solver_util.h"

#include <math.h>

void C_mspt_system_costs::check_parameters_are_set()
{
	if( ms_par.A_sf_refl != ms_par.A_sf_refl ||
		ms_par.site_improv_spec_cost != ms_par.site_improv_spec_cost ||
		ms_par.heliostat_spec_cost != ms_par.heliostat_spec_cost ||
		ms_par.heliostat_fixed_cost != ms_par.heliostat_fixed_cost ||
		
		ms_par.h_tower != ms_par.h_tower ||
		ms_par.h_rec != ms_par.h_rec ||
		ms_par.h_helio != ms_par.h_helio ||
		ms_par.tower_fixed_cost != ms_par.tower_fixed_cost ||
		ms_par.tower_cost_scaling_exp != ms_par.tower_cost_scaling_exp ||
		
		ms_par.A_rec != ms_par.A_rec ||
		ms_par.rec_ref_cost != ms_par.rec_ref_cost ||
		ms_par.A_rec_ref != ms_par.A_rec_ref ||
		ms_par.rec_cost_scaling_exp != ms_par.rec_cost_scaling_exp ||
		
		ms_par.Q_storage != ms_par.Q_storage ||
		ms_par.tes_spec_cost != ms_par.tes_spec_cost ||
		
		ms_par.W_dot_design != ms_par.W_dot_design ||
		ms_par.power_cycle_spec_cost != ms_par.power_cycle_spec_cost ||
		
		ms_par.bop_spec_cost != ms_par.bop_spec_cost ||
		
		ms_par.fossil_backup_spec_cost != ms_par.fossil_backup_spec_cost ||
		
		ms_par.contingency_rate != ms_par.contingency_rate ||
		
		ms_par.total_land_area != ms_par.total_land_area ||
		ms_par.plant_net_capacity != ms_par.plant_net_capacity ||
		ms_par.EPC_land_spec_cost != ms_par.EPC_land_spec_cost ||
		ms_par.EPC_land_perc_direct_cost != ms_par.EPC_land_perc_direct_cost ||
		ms_par.EPC_land_per_power_cost != ms_par.EPC_land_per_power_cost ||
		ms_par.EPC_land_fixed_cost != ms_par.EPC_land_fixed_cost ||
		ms_par.total_land_spec_cost != ms_par.total_land_spec_cost ||
		ms_par.total_land_perc_direct_cost != ms_par.total_land_perc_direct_cost ||
		ms_par.total_land_per_power_cost != ms_par.total_land_per_power_cost ||
		ms_par.total_land_fixed_cost != ms_par.total_land_fixed_cost ||
		ms_par.sales_tax_basis != ms_par.sales_tax_basis ||
		ms_par.sales_tax_rate != ms_par.sales_tax_rate )
	{
		std::string msg = "C_mspt_system_costs initialization failed because not all required parameters were defined"
		"before calculate_costs() was called";
		C_csp_exception(msg,0);
	}

	return;
}


void C_mspt_system_costs::calculate_costs()
{
	check_parameters_are_set();

	ms_out.site_improvement_cost = 
		N_mspt::site_improvement_cost(ms_par.A_sf_refl, ms_par.site_improv_spec_cost);

	ms_out.heliostat_cost = 
		N_mspt::heliostat_cost(ms_par.A_sf_refl, ms_par.heliostat_spec_cost, ms_par.heliostat_fixed_cost);

	ms_out.tower_cost = 
		N_mspt::tower_cost(ms_par.h_tower, ms_par.h_rec, ms_par.h_helio, ms_par.tower_fixed_cost, ms_par.tower_cost_scaling_exp);

	ms_out.receiver_cost = 
		N_mspt::receiver_cost(ms_par.A_rec, ms_par.rec_ref_cost, ms_par.A_rec_ref, ms_par.rec_cost_scaling_exp);

	ms_out.tes_cost = 
		N_mspt::tes_cost(ms_par.Q_storage, ms_par.tes_spec_cost);

	ms_out.power_cycle_cost = 
		N_mspt::power_cycle_cost(ms_par.W_dot_design, ms_par.power_cycle_spec_cost);

	ms_out.bop_cost = 
		N_mspt::bop_cost(ms_par.W_dot_design, ms_par.bop_spec_cost);

	ms_out.fossil_backup_cost = 
		N_mspt::fossil_backup_cost(ms_par.W_dot_design, ms_par.fossil_backup_spec_cost);

	ms_out.direct_capital_precontingency_cost = 
		N_mspt::direct_capital_precontingency_cost(
			ms_out.site_improvement_cost,
			ms_out.heliostat_cost,
			ms_out.tower_cost,
			ms_out.receiver_cost,
			ms_out.tes_cost,
			ms_out.power_cycle_cost,
			ms_out.bop_cost,
			ms_out.fossil_backup_cost);

	ms_out.contingency_cost = 
		N_mspt::contingency_cost(ms_par.contingency_rate, ms_out.direct_capital_precontingency_cost);

	ms_out.total_direct_cost = 
		N_mspt::total_direct_cost(ms_out.direct_capital_precontingency_cost, ms_out.contingency_cost);

	ms_out.total_land_cost = 
		N_mspt::total_land_cost(ms_par.total_land_area, ms_out.total_direct_cost, ms_par.plant_net_capacity,
			ms_par.total_land_spec_cost, ms_par.total_land_perc_direct_cost, ms_par.total_land_per_power_cost, ms_par.total_land_fixed_cost);

	ms_out.epc_and_owner_cost = 
		N_mspt::epc_and_owner_cost(ms_par.total_land_area, ms_out.total_direct_cost, ms_par.plant_net_capacity,
			ms_par.EPC_land_spec_cost, ms_par.EPC_land_perc_direct_cost, ms_par.EPC_land_per_power_cost, ms_par.EPC_land_fixed_cost);

	ms_out.sales_tax_cost = 
		N_mspt::sales_tax_cost(ms_out.total_direct_cost, ms_par.sales_tax_basis, ms_par.sales_tax_rate);

	ms_out.total_indirect_cost = 
		N_mspt::total_indirect_cost(ms_out.total_land_cost, ms_out.epc_and_owner_cost, ms_out.sales_tax_cost);

	ms_out.total_installed_cost = 
		N_mspt::total_installed_cost(ms_out.total_direct_cost, ms_out.total_indirect_cost);

	ms_out.estimated_installed_cost_per_cap = 
		N_mspt::estimated_installed_cost_per_cap(ms_out.total_installed_cost, ms_par.plant_net_capacity);

	return;
}

double N_mspt::site_improvement_cost(double A_refl /*m^2*/, double site_improv_spec_cost /*$/m^2_reflect*/)
{
	return A_refl*site_improv_spec_cost;		//[$]
}

double N_mspt::heliostat_cost(double A_refl /*m^2*/, double heliostat_spec_cost /*$/m^2*/, double heliostate_fixed_cost /*$*/)
{
	return A_refl*heliostat_spec_cost + heliostate_fixed_cost;	//[$]
}

double N_mspt::tower_cost(double h_tower /*m*/, double h_rec /*m*/, double h_helio /*m*/, double tower_fixed_cost /*$*/, double tower_cost_scaling_exp /*-*/)
{
	return tower_fixed_cost * exp(tower_cost_scaling_exp * (h_tower - h_rec / 2.0 + h_helio / 2.0));	//[$]
}

double N_mspt::receiver_cost(double A_rec /*m^2*/, double rec_ref_cost /*$*/, double rec_ref_area /*m^2*/, double rec_cost_scaling_exp /*-*/)
{
	return rec_ref_cost*pow(A_rec/rec_ref_area, rec_cost_scaling_exp);	//[$]
}

double N_mspt::tes_cost(double Q_storage /*MWt-hr*/, double tes_spec_cost /*$/kWt-hr*/)
{
	return Q_storage*1.E3*tes_spec_cost;		//[$]
}

double N_mspt::power_cycle_cost(double W_dot_design /*MWe*/, double power_cycle_spec_cost /*$/kWe*/)
{
	return W_dot_design*1.E3*power_cycle_spec_cost;		//[$]
}

double N_mspt::bop_cost(double W_dot_design /*MWe*/, double bop_spec_cost /*$/kWe*/)
{
	return W_dot_design*1.E3*bop_spec_cost;				//[$]
}

double N_mspt::fossil_backup_cost(double W_dot_design /*MWe*/, double fossil_backup_spec_cost /*$/kWe*/)
{
	return W_dot_design*1.E3*fossil_backup_spec_cost;	//[$]
}

double N_mspt::direct_capital_precontingency_cost(double site_improvement_cost /*$*/,
	double heliostat_cost /*$*/,
	double tower_cost /*$*/,
	double receiver_cost /*$*/,
	double tes_cost /*$*/,
	double power_cycle_cost /*$*/,
	double bop_cost /*$*/,
	double fossil_backup_cost /*$*/)
{
	return site_improvement_cost +
		heliostat_cost +
		tower_cost +
		receiver_cost +
		tes_cost +
		power_cycle_cost +
		bop_cost +
		fossil_backup_cost;	//[$]
}

double N_mspt::contingency_cost(double contingency_rate /*%*/, double direct_capital_precontingency_cost /*$*/)
{
	return contingency_rate/100.0*direct_capital_precontingency_cost;	//[$]
}

double N_mspt::total_direct_cost(double direct_capital_precontingency_cost /*$*/, double contingency_cost /*$*/)
{
	return direct_capital_precontingency_cost + contingency_cost;		//[$]
}

double N_mspt::total_land_cost(double total_land_area /*acres*/, double total_direct_cost /*$*/, double plant_net_capacity /*MWe*/,
	double land_spec_cost /*$/acre*/, double land_perc_direct_cost /*%*/, double land_spec_per_power_cost /*$/We*/, double land_fixed_cost /*$*/)
{
	return total_land_area*land_spec_cost +
		total_direct_cost*land_perc_direct_cost/100.0 +
		plant_net_capacity*1.E6*land_spec_per_power_cost +
		land_fixed_cost;	//[$]
}

double N_mspt::epc_and_owner_cost(double total_land_area /*acres*/, double total_direct_cost /*$*/, double plant_net_capacity /*MWe*/,
	double land_spec_cost /*$/acre*/, double land_perc_direct_cost /*%*/, double land_spec_per_power_cost /*$/We*/, double land_fixed_cost /*$*/)
{
	return total_land_area*land_spec_cost +
		total_direct_cost*land_perc_direct_cost / 100.0 +
		plant_net_capacity*1.E6*land_spec_per_power_cost +
		land_fixed_cost;	//[$]
}

double N_mspt::sales_tax_cost(double total_direct_cost /*$*/, double sales_tax_basis /*% of tot. direct cost*/, double sales_tax_rate /*%*/)
{
	return total_direct_cost*(sales_tax_basis/100.0)*(sales_tax_rate/100.0);
}

double N_mspt::total_indirect_cost(double total_land_cost /*$*/, double epc_and_owner_cost /*$*/, double sales_tax_cost /*$*/)
{
	return total_land_cost + epc_and_owner_cost + sales_tax_cost;	//[$]
}

double N_mspt::total_installed_cost(double total_direct_cost /*$*/, double total_indirect_cost /*$*/)
{
	return total_direct_cost + total_indirect_cost;		//[$]
}

double N_mspt::estimated_installed_cost_per_cap(double total_installed_cost /*$*/, double plant_net_capacity /*$*/)
{
	return total_installed_cost/(plant_net_capacity*1.E3);
}