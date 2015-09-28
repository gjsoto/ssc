#include "common.h"

var_info vtab_standard_financial[] = {

/*   VARTYPE           DATATYPE         NAME                                         LABEL                              UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT,        SSC_NUMBER,      "analysis_period",                           "Analyis period",                                  "years",  "",                      "Financials",      "?=30",                   "INTEGER,MIN=0,MAX=50",          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "federal_tax_rate",                         "Federal tax rate",                                "%",      "",                      "Financials",      "*",                      "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "state_tax_rate",                           "State tax rate",                                  "%",      "",                      "Financials",      "*",                      "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "property_tax_rate",                        "Property tax rate",                               "%",      "",                      "Financials",      "?=0.0",                  "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,     "prop_tax_cost_assessed_percent",            "Percent of pre-financing costs assessed","%","",			  "Financials",			 "?=95",                     "MIN=0,MAX=100",      			"" },
	{ SSC_INPUT,        SSC_NUMBER,     "prop_tax_assessed_decline",                 "Assessed value annual decline",	"%",	 "",					  "Financials",             "?=5",                     "MIN=0,MAX=100",      			"" },
	{ SSC_INPUT,        SSC_NUMBER,      "sales_tax_rate",                           "Sales tax rate",                                  "%",      "",                      "Financials",      "?=0.0",                  "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "real_discount_rate",                       "Real discount rate",                              "%",      "",                      "Financials",      "*",                      "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inflation_rate",                           "Inflation rate",                                  "%",      "",                      "Financials",      "*",                      "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "insurance_rate",                           "Insurance rate",                                  "%",      "",                      "Financials",      "?=0.0",                  "MIN=0,MAX=100",                 "" },

	{ SSC_INPUT,        SSC_NUMBER,      "system_capacity",                          "System nameplate capacity",                       "kW",     "",                      "System",          "*",                      "POSITIVE",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "system_heat_rate",                         "System heat rate",                                "MMBTus/MWh", "",                  "System",          "?=0.0",                  "MIN=0",                                         "" },

var_info_invalid };

var_info vtab_battery_replacement_cost[] = {

	/*   VARTYPE           DATATYPE         NAME                            LABEL                              UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
		{ SSC_INPUT, SSC_ARRAY, "batt_bank_replacement", "Battery bank replacements per year", "number/year", "", "Battery", "", "", "" },
		{ SSC_INPUT, SSC_ARRAY, "batt_replacement_schedule", "Battery bank replacements per year (user specified)", "number/year", "", "Battery", "", "", "" },
		{ SSC_INPUT, SSC_NUMBER, "en_batt", "Enable battery storage model", "0/1", "", "Battery", "?=0", "", "" },
		{ SSC_INPUT, SSC_NUMBER, "batt_replacement_option", "Enable battery replacement?", "0=none,1=capacity based,2=user schedule", "", "Battery", "?=0", "INTEGER,MIN=0,MAX=2", "" },
		{ SSC_INPUT, SSC_NUMBER, "battery_per_kWh", "Battery cost", "$/kWh", "", "Battery", "?=0.0", "", "" },
		{ SSC_INPUT, SSC_NUMBER, "batt_computed_bank_capacity", "Battery bank capacity", "kWh", "", "Battery", "?=0.0", "", "" },

		{ SSC_INPUT, SSC_ARRAY, "batt_replacement_cost", "Battery bank replacement cost", "$/kWh", "", "Battery", "?=0.0", "", "" },
		{ SSC_INPUT, SSC_NUMBER, "batt_replacement_cost_escal", "Battery bank replacement cost escalation", "%/year", "", "Battery", "?=0.0", "", "" },
		{ SSC_OUTPUT, SSC_ARRAY, "cf_battery_replacement_cost", "Battery replacement cost", "$", "", "Cash Flow", "*", "", "" },
		{ SSC_OUTPUT, SSC_ARRAY, "cf_battery_replacement_cost_schedule", "Battery replacement cost schedule", "$/kWh", "", "Cash Flow", "*", "", "" },

		var_info_invalid };



var_info vtab_standard_loan[] = {

/*   VARTYPE           DATATYPE         NAME                            LABEL                              UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT,        SSC_NUMBER,      "loan_term",					"Loan term",					  "years",  "",                      "Loan",            "?=0",                    "INTEGER,MIN=0,MAX=50",          "" },
	{ SSC_INPUT,        SSC_NUMBER,      "loan_rate",					"Loan rate",					  "%",      "",                      "Loan",            "?=0",                    "MIN=0,MAX=100",                 "" },
	{ SSC_INPUT,        SSC_NUMBER,      "debt_fraction",                   "Debt percentage",                "%",      "",                      "Loan",			"?=0",                    "MIN=0,MAX=100",                 "" },

var_info_invalid };

var_info vtab_oandm[] = {
/*   VARTYPE           DATATYPE         NAME                             LABEL                                UNITS      META                 GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	
	{ SSC_INPUT,        SSC_ARRAY,       "om_fixed",                     "Fixed O&M annual amount",           "$/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_fixed_escal",               "Fixed O&M escalation",              "%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_ARRAY,       "om_production",                "Production-based O&M amount",       "$/MWh",   "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_production_escal",          "Production-based O&M escalation",   "%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_ARRAY,       "om_capacity",                  "Capacity-based O&M amount",         "$/kWcap", "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_capacity_escal",            "Capacity-based O&M escalation",     "%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_ARRAY,		 "om_fuel_cost",                 "Fuel cost",                         "$/MMBtu", "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_fuel_cost_escal",           "Fuel cost escalation",              "%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "annual_fuel_usage",        "Fuel usage",                         "kWht",         "",                      "O&M",      "?=0",                     "MIN=0",                                         "" },

	// optional fuel o and m for Biopower - usage can be in any unit and cost is in $ per usage unit
	{ SSC_INPUT,        SSC_NUMBER,      "om_opt_fuel_1_usage",           "Biomass feedstock usage",              "unit",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_ARRAY,		 "om_opt_fuel_1_cost",                 "Biomass feedstock cost",          "$/unit", "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_opt_fuel_1_cost_escal",           "Biomass feedstock cost escalation","%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_opt_fuel_2_usage",           "Coal feedstock usage",              "unit",  "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_ARRAY,		 "om_opt_fuel_2_cost",                 "Coal feedstock cost",          "$/unit", "",                  "O&M",            "?=0.0",                 "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "om_opt_fuel_2_cost_escal",           "Coal feedstock cost escalation","%/year",  "",                  "O&M",            "?=0.0",                 "",                                         "" },


var_info_invalid };

var_info vtab_depreciation[] = {
/*   VARTYPE           DATATYPE         NAME                              LABEL                                 UNITS     META                                      GROUP             REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/

	{ SSC_INPUT,        SSC_NUMBER,      "depr_fed_type",                "Federal depreciation type",           "",       "0=none,1=macrs_half_year,2=sl,3=custom",  "Depreciation",      "?=0",                     "INTEGER,MIN=0,MAX=3",        "" },
	{ SSC_INPUT,        SSC_NUMBER,      "depr_fed_sl_years",            "Federal depreciation straight-line Years",       "years",  "",                                        "Depreciation",      "depr_fed_type=2",                     "INTEGER,POSITIVE",           "" },
	{ SSC_INPUT,        SSC_ARRAY,       "depr_fed_custom",              "Federal custom depreciation",         "%/year", "",                                        "Depreciation",      "depr_fed_type=3",         "",                           "" },
	{ SSC_INPUT,        SSC_NUMBER,      "depr_sta_type",                "State depreciation type",             "",       "0=none,1=macrs_half_year,2=sl,3=custom",  "Depreciation",      "?=0",                     "INTEGER,MIN=0,MAX=3",        "" },
	{ SSC_INPUT,        SSC_NUMBER,      "depr_sta_sl_years",            "State depreciation straight-line years",         "years",  "",                                        "Depreciation",      "depr_sta_type=2",                     "INTEGER,POSITIVE",           "" },
	{ SSC_INPUT,        SSC_ARRAY,       "depr_sta_custom",              "State custom depreciation",           "%/year", "",                                        "Depreciation",      "depr_sta_type=3",         "",                           "" },

var_info_invalid };
	
var_info vtab_tax_credits[] = {
/*   VARTYPE           DATATYPE         NAME                               LABEL                                                UNITS     META                      GROUP                 REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/

	{ SSC_INPUT,        SSC_NUMBER,       "itc_fed_amount",                 "Federal amount-based ITC amount",                         "$",      "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_fed_amount_deprbas_fed",     "Federal amount-based ITC reduces federal depreciation basis",       "0/1",    "",          "Tax Credit Incentives",      "?=1",                       "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_fed_amount_deprbas_sta",     "Federal amount-based ITC reduces state depreciation basis",       "0/1",    "",          "Tax Credit Incentives",      "?=1",                       "BOOLEAN",                     "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "itc_sta_amount",                 "State amount-based ITC amount",                           "$",      "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_sta_amount_deprbas_fed",     "State amount-based ITC reduces federal depreciation basis",         "0/1",    "",          "Tax Credit Incentives",      "?=0",                       "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_sta_amount_deprbas_sta",     "State amount-based ITC reduces state depreciation basis",         "0/1",    "",          "Tax Credit Incentives",      "?=0",                       "BOOLEAN",                     "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "itc_fed_percent",                "Federal percentage-based ITC percent",                    "%",      "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_fed_percent_maxvalue",       "Federal percentage-based ITC maximum value",                 "$",      "",          "Tax Credit Incentives",      "?=1e99",                    "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_fed_percent_deprbas_fed",    "Federal percentage-based ITC reduces federal depreciation basis",   "0/1",    "",          "Tax Credit Incentives",      "?=1",                       "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_fed_percent_deprbas_sta",    "Federal percentage-based ITC reduces state depreciation basis",   "0/1",    "",          "Tax Credit Incentives",      "?=1",                       "BOOLEAN",                     "" },

	{ SSC_INPUT,        SSC_NUMBER,       "itc_sta_percent",               "State percentage-based ITC percent",                      "%",      "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_sta_percent_maxvalue",       "State percentage-based ITC maximum Value",                   "$",      "",          "Tax Credit Incentives",      "?=1e99",                    "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_sta_percent_deprbas_fed",    "State percentage-based ITC reduces federal depreciation basis",     "0/1",    "",          "Tax Credit Incentives",      "?=0",                       "BOOLEAN",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "itc_sta_percent_deprbas_sta",    "State percentage-based ITC reduces state depreciation basis",     "0/1",    "",          "Tax Credit Incentives",      "?=0",                       "BOOLEAN",                     "" },

	{ SSC_INPUT,        SSC_ARRAY,       "ptc_fed_amount",                 "Federal PTC amount",                                      "$/kWh",  "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ptc_fed_term",                   "Federal PTC term",                                        "years",  "",          "Tax Credit Incentives",      "?=10",                      "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ptc_fed_escal",                  "Federal PTC escalation",                                  "%/year", "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	
	{ SSC_INPUT,        SSC_ARRAY,       "ptc_sta_amount",                 "State PTC amount",                                        "$/kWh",  "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ptc_sta_term",                   "State PTC term",                                          "years",  "",          "Tax Credit Incentives",      "?=10",                      "",                            "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ptc_sta_escal",                  "State PTC escalation",                                    "%/year", "",          "Tax Credit Incentives",      "?=0",                       "",                            "" },
	
var_info_invalid };


var_info vtab_payment_incentives[] = {
	/*   VARTYPE           DATATYPE         NAME                          LABEL                                                  UNITS     META                      GROUP                   REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/

	{ SSC_INPUT,        SSC_NUMBER,       "ibi_fed_amount",                "Federal amount-based IBI amount",                    "$",      "",                      "Payment Incentives",      "?=0",                 "",                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_amount_tax_fed",        "Federal amount-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_amount_tax_sta",        "Federal amount-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_amount_deprbas_fed",    "Federal amount-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_amount_deprbas_sta",    "Federal amount-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                 "BOOLEAN",                       "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "ibi_sta_amount",                "State amount-based IBI amount",                    "$",      "",                      "Payment Incentives",      "?=0",                   "",                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_amount_tax_fed",        "State amount-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_amount_tax_sta",        "State amount-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_amount_deprbas_fed",    "State amount-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_amount_deprbas_sta",    "State amount-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                   "BOOLEAN",                       "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "ibi_uti_amount",                "Utility amount-based IBI amount",                    "$",      "",                      "Payment Incentives",      "?=0",                 "",                      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_amount_tax_fed",        "Utility amount-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_amount_tax_sta",        "Utility amount-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_amount_deprbas_fed",    "Utility amount-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                 "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_amount_deprbas_sta",    "Utility amount-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                 "BOOLEAN",                       "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "ibi_oth_amount",                "Other amount-based IBI amount",                    "$",      "",                      "Payment Incentives",      "?=0",                   "",                      "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_amount_tax_fed",        "Other amount-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_amount_tax_sta",        "Other amount-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_amount_deprbas_fed",    "Other amount-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                   "BOOLEAN",                       "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_amount_deprbas_sta",    "Other amount-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",                   "BOOLEAN",                       "" },

	{ SSC_INPUT,        SSC_NUMBER,       "ibi_fed_percent",               "Federal percentage-based IBI percent",                   "%",      "",                      "Payment Incentives",      "?=0.0",           "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_percent_maxvalue",      "Federal percentage-based IBI maximum value",             "$",      "",                      "Payment Incentives",      "?=1e99",          "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_percent_tax_fed",       "Federal percentage-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_percent_tax_sta",       "Federal percentage-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_percent_deprbas_fed",   "Federal percentage-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_fed_percent_deprbas_sta",   "Federal percentage-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },

	{ SSC_INPUT,        SSC_NUMBER,       "ibi_sta_percent",               "State percentage-based IBI percent",                   "%",      "",                      "Payment Incentives",      "?=0.0",           "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_percent_maxvalue",      "State percentage-based IBI maximum value",             "$",      "",                      "Payment Incentives",      "?=1e99",          "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_percent_tax_fed",       "State percentage-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_percent_tax_sta",       "State percentage-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_percent_deprbas_fed",   "State percentage-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_sta_percent_deprbas_sta",   "State percentage-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },

	{ SSC_INPUT,        SSC_NUMBER,       "ibi_uti_percent",               "Utility percentage-based IBI percent",                   "%",      "",                      "Payment Incentives",      "?=0.0",           "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_percent_maxvalue",      "Utility percentage-based IBI maximum value",             "$",      "",                      "Payment Incentives",      "?=1e99",          "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_percent_tax_fed",       "Utility percentage-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_percent_tax_sta",       "Utility percentage-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_percent_deprbas_fed",   "Utility percentage-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_uti_percent_deprbas_sta",   "Utility percentage-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },

	{ SSC_INPUT,        SSC_NUMBER,       "ibi_oth_percent",               "Other percentage-based IBI percent",                   "%",      "",                      "Payment Incentives",      "?=0.0",           "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_percent_maxvalue",      "Other percentage-based IBI maximum value",             "$",      "",                      "Payment Incentives",      "?=1e99",          "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_percent_tax_fed",       "Other percentage-based IBI federal taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_percent_tax_sta",       "Other percentage-based IBI state taxable",              "0/1",    "",                      "Payment Incentives",      "?=1",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_percent_deprbas_fed",   "Other percentage-based IBI reduces federal depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "ibi_oth_percent_deprbas_sta",   "Other percentage-based IBI reduces state depreciation basis",  "0/1",    "",                      "Payment Incentives",      "?=0",             "BOOLEAN",                                         "" },


	{ SSC_INPUT,        SSC_NUMBER,       "cbi_fed_amount",         "Federal CBI amount",                   "$/Watt", "",                      "Payment Incentives",      "?=0.0",                     "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_fed_maxvalue",       "Federal CBI maximum",                  "$",      "",                      "Payment Incentives",      "?=1e99",                    "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_fed_tax_fed",        "Federal CBI federal taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_fed_tax_sta",        "Federal CBI state taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_fed_deprbas_fed",    "Federal CBI reduces federal depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_fed_deprbas_sta",    "Federal CBI reduces state depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "cbi_sta_amount",         "State CBI amount",                   "$/Watt", "",                      "Payment Incentives",      "?=0.0",                     "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_sta_maxvalue",       "State CBI maximum",                  "$",      "",                      "Payment Incentives",      "?=1e99",                    "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_sta_tax_fed",        "State CBI federal taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_sta_tax_sta",        "State CBI state taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_sta_deprbas_fed",    "State CBI reduces federal depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_sta_deprbas_sta",    "State CBI reduces state depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "cbi_uti_amount",         "Utility CBI amount",                   "$/Watt", "",                      "Payment Incentives",      "?=0.0",                     "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_uti_maxvalue",       "Utility CBI maximum",                  "$",      "",                      "Payment Incentives",      "?=1e99",                    "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_uti_tax_fed",        "Utility CBI federal taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_uti_tax_sta",        "Utility CBI state taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_uti_deprbas_fed",    "Utility CBI reduces federal depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_uti_deprbas_sta",    "Utility CBI reduces state depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_NUMBER,       "cbi_oth_amount",         "Other CBI amount",                   "$/Watt", "",                      "Payment Incentives",      "?=0.0",                     "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_oth_maxvalue",       "Other CBI maximum",                  "$",      "",                      "Payment Incentives",      "?=1e99",                    "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_oth_tax_fed",        "Other CBI federal taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_oth_tax_sta",        "Other CBI state taxable",             "0/1",    "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_oth_deprbas_fed",    "Other CBI reduces federal depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "cbi_oth_deprbas_sta",    "Other CBI reduces state depreciation basis", "0/1",    "",                      "Payment Incentives",      "?=0",                       "BOOLEAN",                                         "" },


	{ SSC_INPUT,        SSC_ARRAY,       "pbi_fed_amount",                "Federal PBI amount",           "$/kWh",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_fed_term",                  "Federal PBI term",             "years",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_fed_escal",                 "Federal PBI escalation",       "%",        "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_fed_tax_fed",               "Federal PBI federal taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_fed_tax_sta",               "Federal PBI state taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_ARRAY,       "pbi_sta_amount",                "State PBI amount",           "$/kWh",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_sta_term",                  "State PBI term",             "years",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_sta_escal",                 "State PBI escalation",       "%",        "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_sta_tax_fed",               "State PBI federal taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_sta_tax_sta",               "State PBI state taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_ARRAY,       "pbi_uti_amount",                "Utility PBI amount",           "$/kWh",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_uti_term",                  "Utility PBI term",             "years",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_uti_escal",                 "Utility PBI escalation",       "%",        "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_uti_tax_fed",               "Utility PBI federal taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_uti_tax_sta",               "Utility PBI state taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	
	{ SSC_INPUT,        SSC_ARRAY,       "pbi_oth_amount",                "Other PBI amount",           "$/kWh",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_oth_term",                  "Other PBI term",             "years",    "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_oth_escal",                 "Other PBI escalation",       "%",        "",                      "Payment Incentives",      "?=0",                       "",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_oth_tax_fed",               "Other PBI federal taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },
	{ SSC_INPUT,        SSC_NUMBER,      "pbi_oth_tax_sta",               "Other PBI state taxable",     "0/1",      "",                      "Payment Incentives",      "?=1",                       "BOOLEAN",                                         "" },

var_info_invalid };


var_info vtab_adjustment_factors[] = {
/*   VARTYPE           DATATYPE         NAME                               LABEL                                       UNITS     META                                     GROUP                 REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/

	{ SSC_INPUT,        SSC_NUMBER,      "adjust:constant",               "Constant loss adjustment",                "%",    "",                                     "Loss Adjustments",      "*",                     "MAX=100",                     "" },
	{ SSC_INPUT,        SSC_ARRAY,       "adjust:hourly",                 "Hourly loss adjustments",                 "%",    "",                                     "Loss Adjustments",      "?",                     "LENGTH=8760",                "" },
	{ SSC_INPUT,        SSC_MATRIX,      "adjust:periods",                "Period-based loss adjustments",           "%",    "n x 3 matrix [ start, end, loss ]",    "Loss Adjustments",      "?",                     "COLS=3",                     "" },
	
var_info_invalid };


var_info vtab_technology_outputs[] = {
	/*   VARTYPE           DATATYPE         NAME                               LABEL                                       UNITS     META                                     GROUP                 REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
// instantaneous power at each timestep - consistent with sun position
		{ SSC_OUTPUT, SSC_ARRAY, "gen", "Power generated by system", "kW", "", "Time Series", "*", "", "" },
		var_info_invalid };


adjustment_factors::adjustment_factors( compute_module *cm )
	: m_cm(cm)
{	
}

//adjustment factors changed from derates to percentages jmf 1/9/15
bool adjustment_factors::setup()
{
	float f = (float)m_cm->as_number( "adjust:constant" );
	f = 1 - f / 100; //convert from percentage to factor
	m_factors.resize( 8760, f );

	if ( m_cm->is_assigned("adjust:hourly") )
	{
		size_t n;
		ssc_number_t *p = m_cm->as_array( "adjust:hourly", &n );
		if ( p != 0 && n == 8760 )
		{
			for( size_t i=0;i<8760;i++ )
				m_factors[i] *= (1 - p[i]/100); //convert from percentages to factors
		}
	}

	if ( m_cm->is_assigned("adjust:periods") )
	{
		size_t nr, nc;
		ssc_number_t *mat = m_cm->as_matrix( "adjust:periods", &nr, &nc );
		if ( mat != 0 && nc == 3 )
		{
			for( size_t r=0;r<nr;r++ )
			{
				int start = (int) mat[ nc*r ];
				int end = (int) mat[ nc*r + 1 ];
				float factor = (float) mat[ nc*r + 2 ];
				
				if ( start < 0 || start >= 8760 || end < start )
				{
					m_error = util::format( "period %d is invalid ( start: %d, end %d )", (int)r, start, end );
					continue;
				}

				if ( end >= 8760 ) end = 8759;

				for( int i=start;i<=end;i++ )
					m_factors[i] *= (1 - factor/100); //convert from percentages to factors
			}
		}
	}

	return m_error.length() == 0;
}

float adjustment_factors::operator()( size_t time )
{
	if ( time < m_factors.size() ) return m_factors[time];
	else return 0.0;
}


shading_factor_calculator::shading_factor_calculator()
{
	m_enAzAlt = false;
	m_diffFactor = 1.0;
}

bool shading_factor_calculator::setup( compute_module *cm, const std::string &prefix )
{
	bool ok = true;
	m_diffFactor = 1.0;
	m_en_shading_db = false;
	m_steps_per_hour = 1;
	m_db8 = NULL;

	if (cm->is_assigned(prefix + "shading:shading_db_lookup"))
		m_en_shading_db = cm->as_boolean(prefix + "shading:shading_db_lookup");
	if (m_en_shading_db)
	{
		m_db8 = new DB8_mpp();
		m_db8->init();
	}

	// initialize to 8760x1 for mxh and change based on shading:timestep
	m_beamFactors.resize_fill(8760, 1, 1.0);

	/*
	m_beamFactors.resize( 8760, 1.0 );

	
	for ( size_t j=0;j<8760;j++)
		m_beamFactors[j] = 1.0;

	if ( cm->is_assigned( prefix+"shading:hourly" ) )
	{
		size_t len = 0;
		ssc_number_t *vals = cm->as_array( prefix+"shading:hourly", &len );
		if ( len == 8760 )
		{
			for ( size_t j=0;j<8760;j++)
				m_beamFactors[j] = 1-vals[j]/100;
		}
		else
		{
			ok = false;
			m_errors.push_back("hourly shading beam losses must have 8760 values");
		}
	}
	*/


	if (cm->is_assigned(prefix + "shading:timestep"))
	{
		size_t nrows, ncols;
		ssc_number_t *mat = cm->as_matrix(prefix + "shading:timestep", &nrows, &ncols);
		if (nrows % 8760 == 0)
		{
			m_beamFactors.resize_fill(nrows, ncols, 1.0);
			if (m_en_shading_db) // use percent shaded to lookup in database
			{
				for (size_t r = 0; r < nrows; r++)
					for (size_t c = 0; c < ncols; c++)
						m_beamFactors.at(r, c) = mat[r*ncols + c]; //entered in % shaded 
			}
			else // use unshaded factors to apply to beam
			{
				for (size_t r = 0; r < nrows; r++)
					for (size_t c = 0; c < ncols; c++)
						m_beamFactors.at(r, c) = 1 - mat[r*ncols + c] / 100; //all other entries must be converted from % to factor unshaded for beam
			}
			m_steps_per_hour = nrows / 8760;
		}
		else
		{
			ok = false;
			m_errors.push_back("hourly shading beam losses must be multiple of 8760 values");
		}
	}

	if (!m_en_shading_db)
	{
		if (cm->is_assigned(prefix + "shading:mxh"))
		{
			size_t nrows, ncols;
			ssc_number_t *mat = cm->as_matrix(prefix + "shading:mxh", &nrows, &ncols);
			if (nrows != 12 || ncols != 24)
			{
				ok = false;
				m_errors.push_back("month x hour shading losses must have 12 rows and 24 columns");
			}
			else
			{
				int c = 0;
				for (int m = 0; m < 12; m++)
					for (int d = 0; d < util::nday[m]; d++)
						for (int h = 0; h < 24; h++)
							for (int jj = 0; jj < m_steps_per_hour; jj++)
								m_beamFactors.at(c++,0) *= 1 - mat[m*ncols + h] / 100;
			}
		}

		m_enAzAlt = false;
		if (cm->is_assigned(prefix + "shading:azal"))
		{
			size_t nrows, ncols;
			ssc_number_t *mat = cm->as_matrix(prefix + "shading:azal", &nrows, &ncols);
			if (nrows < 3 || ncols < 3)
			{
				ok = false;
				m_errors.push_back("azimuth x altitude shading losses must have at least 3 rows and 3 columns");
			}

			m_azaltvals.resize_fill(nrows, ncols, 1.0);
			for (size_t r = 0; r < nrows; r++)
				for (size_t c = 0; c < ncols; c++)
				{
					if (r == 0 || c == 0)
						m_azaltvals.at(r, c) = mat[r*ncols + c]; //first row and column contain azimuth by altitude values
					else
						m_azaltvals.at(r, c) = 1 - mat[r*ncols + c] / 100; //all other entries must be converted from % to factor
				}

			m_enAzAlt = true;
		}
	}


	if (cm->is_assigned(prefix + "shading:diff"))
		m_diffFactor = 1 - cm->as_double(prefix + "shading:diff") / 100;

	return ok;
}

std::string shading_factor_calculator::get_error(size_t i)
{
	if( i < m_errors.size() ) return m_errors[i];
	else return std::string("");
}


size_t shading_factor_calculator::get_row_index_for_input(size_t hour, size_t hour_step, size_t steps_per_hour)
{
	// handle different simulation timesteps and shading input timesteps
	size_t ndx = hour * m_steps_per_hour; // m_beam row index for input hour
	int hr_step = 0;
	if (steps_per_hour > 0)
		hr_step = (int)((int)m_steps_per_hour*(int)hour_step/ (int)steps_per_hour);
	if (hr_step >= m_steps_per_hour) hr_step = m_steps_per_hour - 1;
	if (hr_step < 0) hr_step = 0;
	ndx += hr_step;
	return ndx;
}

double shading_factor_calculator::fbeam(size_t hour, double solalt, double solazi, size_t hour_step, size_t steps_per_hour, double gpoa, double dpoa)
{
	double factor = 1.0;
	size_t irow = get_row_index_for_input(hour,hour_step,steps_per_hour);
	if ((irow >= 0) && (irow < m_beamFactors.nrows()))
	{
		if (m_en_shading_db) // shading database lookup
		{
			std::vector<double> shad_fracs;
			for (size_t icol = 0; icol < m_beamFactors.ncols(); icol++)
				shad_fracs.push_back(m_beamFactors.at(irow, icol));
			factor = 1.0-m_db8->get_shade_loss(gpoa, dpoa, shad_fracs);
		}
		else
			//	if (hour >= 0 && hour < m_beamFactors.size())
		{
			factor = m_beamFactors.at(irow, 0);
			//		double factor = m_beamFactors[hour];
			if (m_enAzAlt)
				factor *= util::bilinear(solalt, solazi, m_azaltvals);
		}
	}
	return factor;
}

double shading_factor_calculator::fdiff()
{
	return m_diffFactor;
}

weatherdata::weatherdata( var_data *data_table )
{
	m_startSec = m_stepSec = m_nRecords = 0;
	m_index = 0;

	if ( data_table->type != SSC_TABLE ) 
	{
		m_error = "solar data must be an SSC table variable with fields: "
			"(numbers): lat, lon, tz, elev, "
			"(arrays): year, month, day, hour minute, gh, dn, df, wspd, wdir, tdry, twet, tdew, rhum, pres, snow, alb, aod";
		return;
	}

	m_hdr.lat = get_number( data_table, "lat" );
	m_hdr.lon = get_number( data_table, "lon" );
	m_hdr.tz = get_number( data_table, "tz" );
	m_hdr.elev = get_number( data_table, "elev" );

	int nrec = 0;
	vec year = get_vector( data_table, "year", &nrec );
	vec month = get_vector( data_table, "month", &nrec );
	vec day = get_vector( data_table, "day", &nrec );
	vec hour = get_vector( data_table, "hour", &nrec );
	vec minute = get_vector( data_table, "minute", &nrec );
	vec gh = get_vector( data_table, "gh", &nrec );
	vec dn = get_vector( data_table, "dn", &nrec );
	vec df = get_vector( data_table, "df", &nrec );
	vec wspd = get_vector( data_table, "wspd", &nrec );
	vec wdir = get_vector( data_table, "wdir", &nrec );
	vec tdry = get_vector( data_table, "tdry", &nrec ); 
	vec twet = get_vector( data_table, "twet", &nrec ); 
	vec tdew = get_vector( data_table, "tdew", &nrec ); 
	vec rhum = get_vector( data_table, "rhum", &nrec ); 
	vec pres = get_vector( data_table, "pres", &nrec ); 
	vec snow = get_vector( data_table, "snow", &nrec ); 
	vec alb = get_vector( data_table, "alb", &nrec ); 
	vec aod = get_vector( data_table, "aod", &nrec ); 
	
	m_nRecords = (size_t)nrec;

	int nmult = nrec / 8760;

	// estimate time step
	if ( nmult * 8760 == nrec )
	{
		m_stepSec = 3600 / nmult;
		m_startSec = m_stepSec / 2;
	}
	else if ( m_nRecords%8784==0 )
	{ 
		// Check if the weather file contains a leap day
		// if so, correct the number of nrecords 
		m_nRecords = m_nRecords/8784*8760;
		nmult = m_nRecords/8760;
		m_stepSec = 3600 / nmult;
		m_startSec = m_stepSec / 2;
	}
	else
	{
		return;
	}

	if ( nrec > 0 && nmult >= 1 )
	{
		m_data.resize( nrec );
		for( size_t i=0;i<nrec;i++ )
		{
			weather_record *r = new weather_record;

			if ( i < year.len ) r->year = year.p[i]; 
			else r->year = 2000;

			if ( i < month.len ) r->month = month.p[i];
			else if ( m_stepSec == 3600 && m_nRecords == 8760 ) {
				r->month = util::month_of(i);
			}

			if ( i < day.len ) r->day = day.p[i];
			else if ( m_stepSec == 3600 && m_nRecords == 8760 ) {
				int month = util::month_of( i );
				r->day = util::day_of_month( month, i );
			}

			if ( i < hour.len ) r->hour = hour.p[i];
			else if ( m_stepSec == 3600 && m_nRecords == 8760 ) {
				int day = i / 24;
				int start_of_day = day * 24;
				r->hour = (float)(i - start_of_day);
			}

			if ( i < minute.len ) r->minute = minute.p[i];
			else r->minute = (double)((m_stepSec / 2) / 60);


			if ( i < gh.len ) r->gh = gh.p[i];
			if ( i < dn.len ) r->dn = dn.p[i];
			if ( i < df.len ) r->df = df.p[i];

			if ( i < wspd.len ) r->wspd = wspd.p[i];
			if ( i < wdir.len ) r->wdir = wdir.p[i];

			if ( i < tdry.len ) r->tdry = tdry.p[i];
			if ( i < twet.len ) r->twet = twet.p[i];
			if ( i < tdew.len ) r->tdew = tdew.p[i];

			if ( i < rhum.len ) r->rhum = rhum.p[i];
			if ( i < pres.len ) r->pres = pres.p[i];

			if ( i < snow.len ) r->snow = snow.p[i];
			if ( i < alb.len ) r->alb = alb.p[i];
			if ( i < aod.len ) r->aod = aod.p[i];

			m_data[i] = r;
		}
	}
}

weatherdata::~weatherdata()
{
	for( size_t i=0;i<m_data.size();i++ )
		delete m_data[i];
}

const char *weatherdata::error( size_t idx )
{
	return ( idx == 0 && m_error.size() > 0 ) ? m_error.c_str() : 0;
}

weatherdata::vec weatherdata::get_vector( var_data *v, const char *name, int *maxlen )
{
	vec x;
	x.p = 0;
	x.len = 0;
	if ( var_data *value = v->table.lookup( name ) )
	{
		if ( value->type == SSC_ARRAY )
		{
			x.len = (int) value->num.length();
			x.p = value->num.data();
			if ( maxlen && x.len > *maxlen )
				*maxlen = x.len;
		}
	}

	return x;
}

ssc_number_t weatherdata::get_number( var_data *v, const char *name )
{
	if ( var_data *value = v->table.lookup( name ) )
	{
		if ( value->type == SSC_NUMBER )
			return value->num;
	}

	return std::numeric_limits<ssc_number_t>::quiet_NaN();
}

bool weatherdata::header( weather_header *h )
{
	*h = m_hdr;
	return true;
}

bool weatherdata::read( weather_record *r )
{
	if ( m_index >= 0 && m_index < m_data.size() )
	{
		*r = *m_data[m_index++];
		return true;
	}
	else
		return false;
}

void weatherdata::rewind()
{
	m_index = 0;
}
