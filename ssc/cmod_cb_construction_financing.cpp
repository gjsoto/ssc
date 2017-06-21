#include "core.h"

#include "csp_system_costs.h"
#include "csp_solver_util.h"


static var_info _cm_vtab_cb_construction_financing[] = {

	{ SSC_INPUT,  SSC_NUMBER,   "total_installed_cost",                "Total installed cost",                        "$",    "",     "system costs",           "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_interest_rate1",            "Interest rate, loan 1",                       "%",    "",     "financial parameters",   "*",   "",  "" },   
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_interest_rate2",            "Interest rate, loan 2",                       "%",    "",     "financial parameters",   "*",   "",  "" },   
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_interest_rate3",            "Interest rate, loan 3",                       "%",    "",     "financial parameters",   "*",   "",  "" },   
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_interest_rate4",            "Interest rate, loan 4",                       "%",    "",     "financial parameters",   "*",   "",  "" },   
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_interest_rate5",            "Interest rate, loan 5",                       "%",    "",     "financial parameters",   "*",   "",  "" },   
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_months1",                   "Months prior to operation, loan 1",           "",     "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_months2",                   "Months prior to operation, loan 2",           "",     "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_months3",                   "Months prior to operation, loan 3",           "",     "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_months4",                   "Months prior to operation, loan 4",           "",     "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_months5",                   "Months prior to operation, loan 5",           "",     "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_percent1",                  "Percent of tot. installed cost, loan 1",      "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_percent2",                  "Percent of tot. installed cost, loan 2",      "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_percent3",                  "Percent of tot. installed cost, loan 3",      "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_percent4",                  "Percent of tot. installed cost, loan 4",      "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_percent5",                  "Percent of tot. installed cost, loan 5",      "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_upfront_rate1",             "Upfront fee on principal, loan 1",            "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_upfront_rate2",             "Upfront fee on principal, loan 2",            "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_upfront_rate3",             "Upfront fee on principal, loan 3",            "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_upfront_rate4",             "Upfront fee on principal, loan 4",            "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_INPUT,  SSC_NUMBER,   "const_per_upfront_rate5",             "Upfront fee on principal, loan 5",            "%",    "",     "financial parameters",   "*",   "",  "" },
																												      
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal1",                "Principal, loan 1",                           "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal2",                "Principal, loan 2",                           "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal3",                "Principal, loan 3",                           "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal4",                "Principal, loan 4",                           "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal5",                "Principal, loan 5",                           "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest1",                 "Interest cost, loan 1",                       "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest2",                 "Interest cost, loan 2",                       "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest3",                 "Interest cost, loan 3",                       "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest4",                 "Interest cost, loan 4",                       "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest5",                 "Interest cost, loan 5",                       "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_total1",                    "Total financing cost, loan 1",                "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_total2",                    "Total financing cost, loan 2",                "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_total3",                    "Total financing cost, loan 3",                "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_total4",                    "Total financing cost, loan 4",                "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_total5",                    "Total financing cost, loan 5",                "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_percent_total",	           "Total percent of installed costs, all loans", "%",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_principal_total",           "Total principal, all loans",				  "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "const_per_interest_total",	           "Total interest costs, all loans",			  "$",    "",     "financial parameters",   "*",   "",  "" },
	{ SSC_OUTPUT, SSC_NUMBER,   "construction_financing_cost",         "Total construction financing cost",           "$",    "",     "financial parameters",   "*",   "",  "" },

	var_info_invalid};


class cm_cb_construction_financing : public compute_module
{
public:

	cm_cb_construction_financing()
	{
		add_var_info(_cm_vtab_cb_construction_financing);
	}

	void exec() throw(general_error)
	{
		
		double total_installed_cost     = as_double("total_installed_cost");    
		double const_per_interest_rate1 = as_double("const_per_interest_rate1");
		double const_per_interest_rate2 = as_double("const_per_interest_rate2");
		double const_per_interest_rate3 = as_double("const_per_interest_rate3");
		double const_per_interest_rate4 = as_double("const_per_interest_rate4");
		double const_per_interest_rate5 = as_double("const_per_interest_rate5");
		double const_per_months1        = as_double("const_per_months1");
		double const_per_months2        = as_double("const_per_months2");       
		double const_per_months3        = as_double("const_per_months3");       
		double const_per_months4        = as_double("const_per_months4");       
		double const_per_months5        = as_double("const_per_months5");       
		double const_per_percent1       = as_double("const_per_percent1");      
		double const_per_percent2       = as_double("const_per_percent2");      
		double const_per_percent3       = as_double("const_per_percent3");      
		double const_per_percent4       = as_double("const_per_percent4");      
		double const_per_percent5       = as_double("const_per_percent5");      
		double const_per_upfront_rate1  = as_double("const_per_upfront_rate1"); 
		double const_per_upfront_rate2  = as_double("const_per_upfront_rate2"); 
		double const_per_upfront_rate3  = as_double("const_per_upfront_rate3"); 
		double const_per_upfront_rate4  = as_double("const_per_upfront_rate4"); 
		double const_per_upfront_rate5  = as_double("const_per_upfront_rate5"); 

		double const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5;
		double const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5;
		double const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5;
		double const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost;

		const_per_principal1 = const_per_principal2 = const_per_principal3 = const_per_principal4 = const_per_principal5 =
			const_per_interest1 = const_per_interest2 = const_per_interest3 = const_per_interest4 = const_per_interest5 = 
			const_per_total1 = const_per_total2 = const_per_total3 = const_per_total4 = const_per_total5 =
			const_per_percent_total = const_per_principal_total = const_per_interest_total = construction_financing_cost =
			std::numeric_limits<double>::quiet_NaN();

		N_financial_parameters::construction_financing_total_cost(total_installed_cost,
			const_per_interest_rate1, const_per_interest_rate2, const_per_interest_rate3, const_per_interest_rate4, const_per_interest_rate5,
			const_per_months1, const_per_months2, const_per_months3, const_per_months4, const_per_months5,
			const_per_percent1, const_per_percent2, const_per_percent3, const_per_percent4, const_per_percent5,
			const_per_upfront_rate1, const_per_upfront_rate2, const_per_upfront_rate3, const_per_upfront_rate4, const_per_upfront_rate5,
			const_per_principal1, const_per_principal2, const_per_principal3, const_per_principal4, const_per_principal5,
			const_per_interest1, const_per_interest2, const_per_interest3, const_per_interest4, const_per_interest5,
			const_per_total1, const_per_total2, const_per_total3, const_per_total4, const_per_total5,
			const_per_percent_total, const_per_principal_total, const_per_interest_total, construction_financing_cost);

		assign("const_per_principal1",        const_per_principal1);         
		assign("const_per_principal2",        const_per_principal2);      
		assign("const_per_principal3",        const_per_principal3);      
		assign("const_per_principal4",        const_per_principal4);      
		assign("const_per_principal5",        const_per_principal5);      
		assign("const_per_interest1",         const_per_interest1);       
		assign("const_per_interest2",         const_per_interest2);       
		assign("const_per_interest3",         const_per_interest3);       
		assign("const_per_interest4",         const_per_interest4);       
		assign("const_per_interest5",         const_per_interest5);       
		assign("const_per_total1",            const_per_total1);          
		assign("const_per_total2",            const_per_total2);          
		assign("const_per_total3",            const_per_total3);          
		assign("const_per_total4",            const_per_total4);          
		assign("const_per_total5",            const_per_total5);          
		assign("const_per_percent_total",	  const_per_percent_total);	 
		assign("const_per_principal_total",   const_per_principal_total); 
		assign("const_per_interest_total",	  const_per_interest_total);	 
		assign("construction_financing_cost", construction_financing_cost);


		/* Useful for lk script:
		outln("const_per_principal1 = ",      var("const_per_principal1")           );  
		outln("const_per_principal2 = ",      var("const_per_principal2")      		);
		outln("const_per_principal3 = ",      var("const_per_principal3")			);
		outln("const_per_principal4 = ",      var("const_per_principal4")			);
		outln("const_per_principal5 = ",      var("const_per_principal5")			);
		outln("const_per_interest1 = ",       var("const_per_interest1")			);
		outln("const_per_interest2 = ",       var("const_per_interest2")			);
		outln("const_per_interest3 = ",       var("const_per_interest3")			);
		outln("const_per_interest4 = ",       var("const_per_interest4")			);
		outln("const_per_interest5 = ",       var("const_per_interest5")			);
		outln("const_per_total1 = ",          var("const_per_total1")				);
		outln("const_per_total2 = ",          var("const_per_total2")				);
		outln("const_per_total3 = ",          var("const_per_total3")				);
		outln("const_per_total4 = ",          var("const_per_total4")				);
		outln("const_per_total5 = ",          var("const_per_total5")				);
		outln("const_per_percent_total = ",	  var("const_per_percent_total")		);
		outln("const_per_principal_total = ",   var("const_per_principal_total")	);
		outln("const_per_interest_total = ",	var("const_per_interest_total")		);
		outln("construction_financing_cost = ", var("construction_financing_cost")	);
		*/
				
	}

};

DEFINE_MODULE_ENTRY(cb_construction_financing, "Construction financing cost calculations", 0)
