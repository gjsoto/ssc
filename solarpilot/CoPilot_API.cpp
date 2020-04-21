#include <fstream>
#include "rapidxml.hpp"

#include "CoPilot_API.h"
#include "IOUtil.h"
#include "shared/lib_weatherfile.h"


static void EditorOutput(const char *msg)
{
    //SPFrame::Instance().ScriptMessageOutput(msg);
}

//static void APICallback(simulation_info* s, void* data)
//{
//
//    python_callback(void*, sp_data_t*, ...)
//}

struct api_helper
{
    SolarField solarfield;
    var_map variables;
    sim_results results;
    LayoutSimThread* simthread;
    SimControl sim_control;

    api_helper()
    {
        variables.reset();
        solarfield.Create(variables);
        results.clear();
        simthread = 0;
        //solarfield.getSimInfoObject()->setCallbackFunction()
    };
};


SPEXPORT sp_data_t sp_data_create()
{
    return static_cast<sp_data_t> (new api_helper);
}

SPEXPORT void sp_data_free(sp_data_t p_data)
{
    api_helper *mc = static_cast<api_helper*>(p_data);

    if (p_data)
        delete p_data;
}

SPEXPORT void sp_set_value(sp_data_t p_data, const char* name, sp_number_t v)
{
    api_helper *mc = static_cast<api_helper*>(p_data);
    
    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //create a string copy of the variable name
    std::string sname = (std::string)name;

    //if it's a combo, make sure the specified combo choice exists
    if (mc->variables._varptrs.at(sname)->ctype == "combo")
    {
        std::string svalue = my_to_string(v);

        std::vector< std::string > cbchoices = mc->variables._varptrs.at(sname)->combo_get_choices();
        if (std::find(cbchoices.begin(), cbchoices.end(), sname) != cbchoices.end())
        {
            //valid variable and selection
            mc->variables._varptrs.at(sname)->set_from_string(svalue.c_str());
        }
        else
        {
            throw std::invalid_argument("Invalid variable choice for \"" + sname + "\": \"" + my_to_string(v) + "\" is not a valid option.");
        }
    }
    else
    {
        //no problems, just set the variable
        std::string svalue = my_to_string(v);
        mc->variables._varptrs.at(sname)->set_from_string(svalue.c_str());
        return;
    }
}

SPEXPORT void sp_set_string(sp_data_t p_data, const char *name, const char *value)
{
    api_helper *mc = static_cast<api_helper*>(p_data);
    std::string sname = (std::string)name;
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
    {
        throw std::invalid_argument("No such variable: " + sname);
    }
    else
    {
        mc->variables._varptrs.at(sname)->set_from_string(value);
    }
}


/** Assigns value of type SP_VEC_DOUBLE */
SPEXPORT void sp_set_array(sp_data_t p_data, const char *name, sp_number_t *pvalues, int length)
{
    api_helper *mc = static_cast<api_helper*>(p_data);

    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    if (mc->variables._varptrs.at(name)->dattype != SP_DATTYPE::SP_VEC_DOUBLE)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_set_array. " + __FILE__ + ":" + my_to_string(__LINE__));

    //create a string copy of the variable name
    std::string sname = (std::string)name;

    //collect the array data into a string with correct format

    std::stringstream array_string;
    
    for (size_t i = 0; i < length; i++)
        array_string << pvalues[i] << ",";
    
    //assign the array
    mc->variables._varptrs.at(name)->set_from_string(array_string.str().c_str());

}

/** Assigns value of type @a SSC_MATRIX . Matrices are specified as a continuous array, in row-major order.  Example: the matrix [[5,2,3],[9,1,4]] is stored as [5,2,3,9,1,4]. */
SPEXPORT void sp_set_matrix(sp_data_t p_data, const char *name, sp_number_t *pvalues, int nrows, int ncols)
{
    /*
    rows separated by ';'
    cols separated by ','
    */

    api_helper *mc = static_cast<api_helper*>(p_data);

    //make sure the specified variable exists
    if (mc->variables._varptrs.find(name) == mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    if (mc->variables._varptrs.at(name)->dattype != SP_DATTYPE::SP_MATRIX_T)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_set_matrix. " + __FILE__ + ":" + my_to_string(__LINE__));

    //create a string copy of the variable name
    std::string sname = (std::string)name;

    //collect the array data into a string with correct format

    std::stringstream matrix_string;

    for (size_t i = 0; i < nrows; i++)
    {
        for (size_t j = 0; j < ncols; j++)
        {
            matrix_string << pvalues[i*ncols + j] << ",";
        }
        matrix_string << ';';
    }
        
    //assign the matrix
    mc->variables._varptrs.at(name)->set_from_string(matrix_string.str().c_str());

}


SPEXPORT sp_number_t sp_get_number(sp_data_t p_data, const char* name)
{
    /*
    Return value of type double, int, or bool. Bools are returned as 0 or 1.
    */

    api_helper *mc = static_cast<api_helper*>(p_data);

    if (mc->variables._varptrs.find(name) != mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_DOUBLE || dattype != SP_DATTYPE::SP_INT)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_set_matrix. " + __FILE__ + ":" + my_to_string(__LINE__));


    spbase *var = mc->variables._varptrs[name];

    switch (var->dattype)
    {
    case SP_INT:
    {
        spvar<int> *v = static_cast<spvar<int>*>(var);
        return (sp_number_t)v->val;
    }
    case SP_DOUBLE:
    {
        spvar<double> *v = static_cast<spvar<double>*>(var);
        return (sp_number_t)v->val;
    }
    case SP_BOOL:
    {
        spvar<bool> *v = static_cast<spvar<bool>*>(var);
        return (sp_number_t)(v->val ? 1. : 0.);
    }


    
    default:
        break;
    }

}


/** Returns the value of a @a SP_STRING variable with the given name. */
SPEXPORT const char *ssc_data_get_string(sp_data_t p_data, const char *name)
{
    api_helper *mc = static_cast<api_helper*>(p_data);

    if (mc->variables._varptrs.find(name) != mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_DOUBLE || dattype != SP_DATTYPE::SP_INT)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_set_matrix. " + __FILE__ + ":" + my_to_string(__LINE__));


    return mc->variables._varptrs[name]->as_string().c_str();
}


/** Returns the value of a @a SSC_ARRAY variable with the given name. */
SPEXPORT void sp_get_array(sp_data_t p_data, const char *name, sp_number_t* values, int *length)
{
    /*
    Populates 'value_array' with 'length' entries
    */

    api_helper *mc = static_cast<api_helper*>(p_data);

    if (mc->variables._varptrs.find(name) != mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_VEC_DOUBLE || dattype != SP_DATTYPE::SP_VEC_INTEGER)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_get_array. " + __FILE__ + ":" + my_to_string(__LINE__));

    spbase *var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector to a vector<double>
    std::vector<double> Vd;
    spbase::_setv(varstr, Vd);

    //allocate space at the value_array pointer
    values = new sp_number_t[(int)Vd.size()];
    //set length for return
    *length = (int)Vd.size();

    //convert to to return format
    for (size_t i = 0; i < *length; i++)
        values[i] = Vd.at(i);

}


SPEXPORT void sp_get_matrix(sp_data_t p_data, const char *name, sp_number_t* values, int *ncols, int *nrows)
{
    /*
    Populates 'value_array' with 'length' entries
    */

    api_helper *mc = static_cast<api_helper*>(p_data);

    if (mc->variables._varptrs.find(name) != mc->variables._varptrs.end())
        throw std::invalid_argument("No such variable: " + std::string(name));

    //make sure the data type of the variable provided matches the internal data type
    int dattype = mc->variables._varptrs.at(name)->dattype;
    if (dattype != SP_DATTYPE::SP_MATRIX_T)
        throw std::runtime_error("Data type of " + std::string(name) + " is not compatible with sp_get_matrix. " + __FILE__ + ":" + my_to_string(__LINE__));

    spbase *var = mc->variables._varptrs[name];
    std::string varstr = var->as_string();

    //convert the string formatted vector 
    matrix_t<double> Md;
    spbase::_setv(varstr, Md);

    //allocate space at the value_array pointer
    values = new sp_number_t[(int)(Md.nrows()*Md.ncols())];
    //set lengths for return
    *ncols = Md.ncols();
    *nrows = Md.nrows();

    //convert to to return format
    for (size_t i = 0; i < *nrows; i++)
        for (size_t j = 0; j < *ncols; i++)
        values[j + (*ncols)*i] = Md.at(i,j);

}

SPEXPORT void sp_reset_geometry(sp_data_t p_data)
{
    /*
	Reset the system geometry to the SolarPILOT default values, clearing any changes or variable settings.
	Returns: (void):null
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    mc->variables.reset();
}

SPEXPORT int sp_add_receiver(sp_data_t p_data, const char* receiver_name)
{
    /*
	Add a new receiver, returning the unique index
	Returns: (string:name[, boolean:make selection]):integer
    */

    api_helper* mc = static_cast<api_helper*>(p_data);

    std::string tname = std::string(receiver_name);
    var_map* V = &mc->variables;

    // check to make sure this isn't a duplicate. Each item needs a unique name
    bool dupe = false;
    for (unsigned int i = 0; i < V->recs.size(); i++)
    {
        if (tname == V->recs.at(i).rec_name.val)
            dupe = true;
    }
    if (dupe)
    {
        throw std::runtime_error("Please enter a unique name for this geometry.");
        return;
    }

    //Add a receiver
    int ind = V->recs.size();

    V->add_receiver(ind);
    V->recs[ind].rec_name.val = tname;
    
    //Re-create the solar field object
    mc->solarfield.Create(*V);
    
    //F.UpdateReceiverUITemplates();   // TODO: Need to confirm no variables are changed.

    //update the input display
        //F.UpdateCalculatedGUIValues();    // Updates clculated parameters
    mc->solarfield.updateAllCalculatedParameters(*V);

    mc->solarfield.updateCalculatedReceiverPower(*V);  // (line 1783) unsure if this is necessary

    //cxt.result().assign((double)SPFrame::Instance().GetVariablesObject()->recs.back().id.val);
    return V->recs.back().id.val; 
}

SPEXPORT int sp_drop_receiver(sp_data_t p_data, const char* receiver_name)
{
    /*
	Drop a receiver from the current solar field
	Returns: (integer:index)
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    std::string tname = lower_case(std::string(receiver_name));

    for (size_t i = 0; i < V->hels.size(); i++)
    {
        if (tname == lower_case(V->recs.at(i).rec_name.val))
        {
            //delete the item
            V->drop_receiver(V->recs.at(i).id.val);
            mc->solarfield.Create(*V);
            //cxt.result().assign(1.);
            return 1.;
        }
    }

    //cxt.result().assign(0.);
    return 0.;
}

SPEXPORT int sp_add_heliostat_template(sp_data_t p_data, const char* heliostat_name)
{
    /*
	Add a new heliostat template that can be used in the layout.
	Returns: (string:template name):integer
    */

    //Add a heliostat
    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    //string tname = cxt.arg(0).as_string();
    std::string tname = std::string(heliostat_name);
    bool dupe = false;
    for (unsigned int i = 0; i < V->hels.size(); i++)
    {
        if (tname == V->hels.at(i).helio_name.val)
            dupe = true;
    }
    if (dupe)
    {
        throw std::runtime_error("Please enter a unique name for this heliostat template.");
        return;
    }

    int ind = V->hels.size();
    V->add_heliostat(ind);
    V->hels.back().helio_name.val = tname;
    //Re-create the solar field object
    mc->solarfield.Create(*V);
    // F.GetSolarFieldObject()->Create(*V);

    //cxt.result().assign((double)SPFrame::Instance().GetVariablesObject()->hels.back().id.val);
    return V->hels.back().id.val;
}

SPEXPORT int sp_drop_heliostat_template(sp_data_t p_data, const char* heliostat_name)
{
    /*
	Delete (drop) the specified heliostat template from the current setup. Returns true if successful.
	Returns: (string:template name):bool
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    std::string tname = lower_case(std::string(heliostat_name));

    for (size_t i = 0; i < V->hels.size(); i++)
    {
        if (tname == lower_case(V->hels.at(i).helio_name.val))
        {
            //delete the item
            V->drop_heliostat(V->hels.at(i).id.val);
            //SF->Create(*V);
            mc->solarfield.Create(*V);
            //cxt.result().assign(1.);
            return 1.;
        }
    }

    //cxt.result().assign(0.);
    return 0.;
}

SPEXPORT int sp_update_geometry(sp_data_t p_data)
{
    /*	Refresh the solar field, receiver, or ambient condition settings based on the current parameter settings.	Returns: (void):boolean    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;
    SolarField* SF = &mc->solarfield;
    
    if (SF->getHeliostats()->size() == 0)
    {
        //no layout exists, so we should be calling the 'run_layout' method instead
        std::runtime_error("No layout exists, so the 'update_geometry' function cannot be executed. Please first create or import a layout using 'run_layout'.");
        //cxt.result().assign(0.);
        return 0.;
    }

    std::string weatherfile_str = std::string(V->amb.weather_file.val);

    Ambient::readWeatherFile(*V);
    
    //Saving local verison of weather data
    weatherfile wf;
    if (!wf.open(weatherfile_str)) return; //error

    //Update the weather data
    std::string linef = "%d,%d,%d,%.2f,%.1f,%.1f,%.1f";
    char cline[300];

    int nrec = (int)wf.nrecords();

    ArrayString local_wfdat;
    local_wfdat.resize(nrec);

    weather_record wrec;
    for (int i = 0; i < nrec; i++)
    {
        //int year, month, day, hour;
        wf.read(&wrec);
        sprintf(cline, linef.c_str(), wrec.day, wrec.hour, wrec.month, wrec.dn, wrec.tdry, wrec.pres / 1000., wrec.wspd);
        std::string line(cline);
        local_wfdat.at(i) = line;
    }

    //Update the design method box.. this actually updates both the map values and GUI. Probably should fix this sometime..
    //F.UpdateDesignSelect(V->sf.des_sim_detail.mapval(), *V);
        // Function seems to only update var_map with simulation data through GenearateSimulationWeatherData()
    interop::GenerateSimulationWeatherData(*V, V->sf.des_sim_detail.mapval(), local_wfdat);

    //Set up the solar field
    SF->Clean();
    SF->Create(*V);

    try
    {
        SolarField::PrepareFieldLayout(mc->solarfield, 0, true);

        if (mc->solarfield.ErrCheck())
        {
            std::runtime_error("An error occurred when preparing the updated field geometry in the call 'update_geometry'.");
            //cxt.result().assign(0.);
            return 0.;
        }

        SF->calcHeliostatArea();
        SF->updateAllCalculatedParameters(*V);

        double azzen[2];
        mc->solarfield.CalcDesignPtSunPosition(V->sf.sun_loc_des.mapval(), azzen[0], azzen[1]);
        Vect sun = Ambient::calcSunVectorFromAzZen(azzen[0] * D2R, azzen[1] * D2R);

        SF->updateAllTrackVectors(sun);

        if (SF->ErrCheck())
        {
            std::runtime_error("An error occurred when preparing the updated field geometry in the call 'update_geometry'.");
            //cxt.result().assign(0.);
            return 0.;
        }
    }
    catch (std::exception &e)
    {
        std::runtime_error("An error occurred when preparing the updated field geometry in the call 'update_geometry'. Error:\n");
        std::runtime_error(e.what());
        //cxt.result().assign(0.);
        return 0.;
    }
    catch (...)
    {
        std::runtime_error("Unknown error when executing 'update_geometry'.");
        //cxt.result().assign(0.);
        return 0.;
    }

    //cxt.result().assign(1.);
    return 1.;
}

SPEXPORT bool sp_assign_layout(sp_data_t p_data, sp_number_t *pvalues, int nrows, int ncols, int nthreads = 0) //, bool save_detail = true)
{
    /*
    Run layout with specified positions. User specifies layout positions in the following format, where first 4 columns are required:
        "<template (int)> <location X> <location Y> <location Z> <x focal length> <y focal length> <cant i> <cant j> <cant k> <aim X> <aim Y> <aim Z>" (array:positions)
    Returns: boolean
    */
    api_helper* mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    //assign layout positions
    V->sf.layout_data.val.clear();

    std::stringstream heliodata;

    size_t npos = nrows;
    for (size_t i = 0; i < npos; i++)
    {
        for (size_t j = 0; j < 12; j++)
        {
            if (j > ncols - 1)
                heliodata << "NULL";
            else
                heliodata << pvalues[j + ncols * i];
            heliodata << (j < 11 ? "," : ";");
        }

    }
    V->sf.layout_data.val = heliodata.str();

    //user specified layout
    V->sf.layout_method.combo_select_by_mapval(var_solarfield::LAYOUT_METHOD::USERDEFINED);

    //TODO: Will this work?
    bool simok = sp_generate_layout(p_data, nthreads);

    return simok;
}


SPEXPORT bool sp_generate_layout(sp_data_t p_data, int nthreads = 0) //, bool save_detail = true)
{
    /*
    Create a solar field layout. Options include 'nthreads':integer (default All),'save_detail':boolean (default True)",
    run layout without specified positions SolarPILOT generates layout positions ([table:options])
    Returns: boolean
        
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;
    SolarField* SF = &mc->solarfield;

    SimControl* SC = &mc->sim_control;
    LayoutSimThread* SThread = mc->simthread;

    if (nthreads!=0)
        SC->SetThreadCount(nthreads);    

    Ambient* A = SF->getAmbientObject();
    A->readWeatherFile(*V);

    //TODO:Is this needed?  I dont think so. It updates Design page in UI
    //F.UpdateDesignSelect(V->sf.des_sim_detail.mapval(), *V);

    SF->Clean();
    SF->Create(*V);
    
    bool simok = interop::DoManagedLayout(*SC, *SF, *V, SThread);        //Returns TRUE if successful

    return simok;
}

SPEXPORT bool sp_get_layout_info(sp_data_t p_data, sp_number_t *layoutinfo, sp_number_t* nhelio, sp_number_t* ncol)
{
    /*
    Get information regarding the heliostat field layout. Returns matrix with each row corresponding to a heliostat.
        "Information includes: [index, position-x, position-y, position-z, template_id, ranking metric value]
    Returns: (void):table
    */
    api_helper* mc = static_cast<api_helper*>(p_data);
    //var_map* V = &mc->variables;

    SolarField* SF = &mc->solarfield;
    Hvector* hels = SF->getHeliostats();

    htemp_map htemps = *SF->getHeliostatTemplates();

    *nhelio = hels->size();
    unsigned int ncol_i = 6;
    *ncol = (sp_number_t)ncol_i;

    layoutinfo = new sp_number_t[(*nhelio) * (*ncol)];
    //  TODO: Free memory somewhere -> look to see what SSC API does for this

    for (size_t i = 0; i < (int)hels->size(); i++)
    {
        sp_point* loc = hels->at(i)->getLocation();

        layoutinfo[i * ncol_i + 0] = hels->at(i)->getId();
        layoutinfo[i * ncol_i + 1] = loc->x;
        layoutinfo[i * ncol_i + 2] = loc->y;
        layoutinfo[i * ncol_i + 3] = loc->z;
        layoutinfo[i * ncol_i + 4] = hels->at(i)->getMasterTemplate()->getVarMap()->id.val;
        layoutinfo[i * ncol_i + 5] = hels->at(i)->getRankingMetricValue();

    }

    return true;
}

SPEXPORT void var_free_memory(sp_number_t* varptr)
{
    delete[] varptr;
};

SPEXPORT bool sp_simulate(sp_data_t p_data, int nthreads = 1, bool save_detail = true, bool update_aimpoints = true)
//SPEXPORT void sp_simulate(sp_data_t p_data)
{
    /*
    Calculate heliostat field performance. Options include 'nthreads':integer (default All),
    'save_detail':boolean (default True), 'update_aimpoints':boolean (default True)
    Returns: [table:options]):boolean
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    SolarField* SF = &mc->solarfield;
    var_map* V = &mc->variables;
    sim_results* res = &mc->results;
    SimControl* SC = &mc->sim_control;
    
    if (nthreads != 1)
        SC->SetThreadCount(nthreads);
    
    if (update_aimpoints != true)
        V->flux.aim_method.combo_select("Keep existing");

    //Which type of simulation is this?
    int simtype = V->flux.flux_model.mapval();    //0=Delsol, 1=Soltrace

    //Set up field, update aimpoints, and simulate at the performance sun position
    Hvector* helios = SF->getHeliostats();

    if (!interop::PerformanceSimulationPrep(*SF, *helios, simtype)) return false;

    Vect sun = Ambient::calcSunVectorFromAzZen(V->flux.flux_solar_az.Val() * D2R, (90. - V->flux.flux_solar_el.Val())*D2R);

    SF->calcHeliostatShadows(sun);
    if (SF->ErrCheck()) return false;

    res->clear();
    res->resize(1);

    //TODO: Start a timer -> this uses wxTimerRunner
    //F.StartSimTimer();
    std::clock_t start;
    double duration;

    start = std::clock();

    //Which type of simulation?
    bool ok;
    switch (simtype)
    {
    case var_fluxsim::FLUX_MODEL::HERMITE_ANALYTICAL:
        ok = interop::HermiteFluxSimulationHandler(*res, *SF, *helios);
        break;
    case var_fluxsim::FLUX_MODEL::SOLTRACE:
        //ok = F.SolTraceFluxSimulation(*SF, *V, *helios);
        std::runtime_error("SOLTRACE is currently not supported in API.");
        //TODO: Move SolTraceFluxSimulation to interop -> Requires multi-threading
        break;
    default:
        ok = false;
        break;
    }

    //F.StopSimTimer();
    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    // TODO: where does this need to be stored?
        // Based on an inital look through, the time appears only on GUI. 

    SF->getSimInfoObject()->Reset();

    //F.GetFluxPlotObject()->SetPlotData(*SF, *helios, 0);
    //F.GetFieldPlotObject()->SetPlotData(*SF, FIELD_PLOT::EFF_TOT);

    return true;
}

SPEXPORT const char *sp_summary_results(sp_data_t p_data)
{
    /*
	Return an array of tables with summary results from each simulation. The array length is greater than 1 for multiple-receiver simulations.
	Returns: (void):array
    */
    
    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    grid_emulator_base table;
    sim_results* results = &mc->results;

    std::string ret;    //return string

    if (results->size() < 1)
    {
        EditorOutput("No simulation summary results exist. Please simulate first.");
        ret = "Failure";
        return ret.c_str();
    }
    
    // intialize size of string vector
    //std::vector<std::string> r;
    //std::vector<std::vector<std::string>> rvec(V->recs.size(), r);

    // for multiple receivers
    for (size_t i = 0; i < V->recs.size(); i++)
    {
        interop::CreateResultsTable(results->at(i), table);
        //F.CreateResultsTable(results->at(i), table);
        
        unordered_map<std::string, double> res_map;

        for (int j = 0; j < table.GetNumberRows(); j++)
        {
            ret.append(table.GetRowLabelValue(j) + ", " + table.GetCellValue(j, 1) + "\n");
            res_map.insert({ table.GetRowLabelValue(j), std::stod(table.GetCellValue(j, 1)) });
        }

        //add a few more summary results
        bool is_soltrace = res_map.find("Shadowing and Cosine efficiency") != res_map.end();

        double Qwf;
        double Qin = Qwf = res_map.at("Power incident on field");

        if (is_soltrace)
        {
            /*
            soltrace
            for this option, the "Shadowing and Cosine efficiency" is already calculated by the
            program. Just make sure the Shading and Cosine efficiencies aren't double counted.
            */
            ret.append("Shading efficiency, 100.\n");
            ret.append("Cosine efficiency, 100.\n");
            ret.append("Shading loss, 0.\n");
            ret.append("Cosine loss, 0.\n");

            double eta_sc = res_map.at("Shadowing and Cosine efficiency") / 100.;
            Qwf *= eta_sc;
        }
        else
        {
            //hermite
            double eta_sc = res_map.at("Shading efficiency") * res_map.at("Cosine efficiency") / 100.;
            ret.append("Shadowing and Cosine efficiency, " + std::to_string(eta_sc) + "\n");

            double eta_s = res_map.at("Shading efficiency") / 100.;
            Qwf *= eta_s;
            ret.append("Shading loss, " + std::to_string(Qin*(1. - eta_s)) + "\n");
            double eta_c = res_map.at("Cosine efficiency") / 100.;
            ret.append("Cosine loss, " + std::to_string(Qwf*(1 - eta_c)) + "\n");
            Qwf *= eta_c;
        }
        ret.append("Shadowing and Cosine loss, " + std::to_string(Qin - Qwf) + "\n");

        double eta_r = res_map.at("Reflection efficiency") / 100.;
        ret.append("Reflection loss, " + std::to_string(Qwf * (1. - eta_r)) + "\n");
        Qwf *= eta_r;
        double eta_b = res_map.at("Blocking efficiency") / 100.;
        ret.append("Blocking loss, " + std::to_string(Qwf*(1. - eta_b)) + "\n");
        Qwf *= eta_b;
        double eta_i = res_map.at("Image intercept efficiency") / 100.;
        ret.append("Image intercept loss, " + std::to_string(Qwf*(1. - eta_i)) + "\n");
        Qwf *= eta_i;
        double eta_a = res_map.at("Absorption efficiency") / 100.;
        ret.append("Absorption loss, " + std::to_string(Qwf*(1. - eta_a)) + "\n");

        ret.append("Receiver name, " + i == 0 ? "All receivers" : results->at(i).receiver_names.front() + "\n\n\n");
    }

    return ret.c_str();
}

SPEXPORT sp_number_t *sp_detail_results(sp_data_t p_data, int* nrows, int* ncols, const char* header, sp_number_t* selhel = NULL, int nselhel = 0)
{
    /*

    returns a vector with hash entries for each heliostat

    id(integer), location (array), aimpoint (array), tracking_vector (array), layout_metric (double), power_to_receiver (double),
        efficiency (double), cosine (double), intercept (double), reflectance (double), attenuation (double),
        blocking (double), shading (double), clouds (double)


        "[Only valid for Hermite (analytical) simulation engine.]\n"
        "Return an array with detailed heliostat-by-heliostat results from a simulation. "
        "Each entry in the array is a table with entries as follows:\n"
        "{ id(integer), location (array), aimpoint (array), tracking_vector (array), "
        "layout_metric (double), "
        "power_to_receiver (double), "
        "power_reflected (double), "
        "energy (double),"
        "annual efficiency (double),"
        "total efficiency (double), "
        "cosine (double), "
        "intercept (double), "
        "reflectance (double), "
        "attenuation (double), "
        "blocking (double), "
        "shading (double), "
        "clouds (double) }",
        "([array:selected heliostat indices]):array");
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    SolarField* SF = &mc->solarfield;

    std::vector<double>* ret;   //return vector

    if (SF->getHeliostats()->size() > 0)
    {
        Hvector helio_select;

        if (selhel != NULL) //use selected heliostats
        {
            //expecting an int array
            if (typeid(selhel[0]) != typeid(int))
                return NULL;

            std::vector<int> ids;
            ids.reserve(nselhel);

            for (size_t i = 0; i < nselhel; i++)
                ids.push_back(selhel[i]);

            helio_select.reserve(nselhel);

            unordered_map<int, Heliostat*> *heliosid = SF->getHeliostatsByID();

            for (size_t i = 0; i < ids.size(); i++)
                helio_select.push_back(heliosid->at(ids.at(i)));

        }
        else    //use all heliostats
        {
            helio_select.assign(SF->getHeliostats()->begin(), SF->getHeliostats()->end());
        }

        *nrows = helio_select.size();
        *ncols = 23;  //number of results in table

        //loop through selected heliostats, gathering information
        ret->reserve((*nrows) * (*ncols));

        /*
        std::vector<double>* ret;
        ret = new std::vector<double>[(*nrows) * (*ncols)];
        ret->reserve((*nrows) * (*ncols));
        */

        for (size_t i = 0; i < helio_select.size(); i++)
        {
            //select specific heliostat object
            Heliostat* H = helio_select.at(i);

            ret->push_back(H->getId());
            ret->push_back(H->getLocation()->x);
            ret->push_back(H->getLocation()->y);
            ret->push_back(H->getLocation()->z);
            
            ret->push_back(H->getAimPoint()->x);
            ret->push_back(H->getAimPoint()->y);
            ret->push_back(H->getAimPoint()->z);
            
            ret->push_back(H->getTrackVector()->i);
            ret->push_back(H->getTrackVector()->j);
            ret->push_back(H->getTrackVector()->k);
            
            ret->push_back(H->getRankingMetricValue());
            ret->push_back(H->getPowerToReceiver() / 1000.);  //kW
            ret->push_back(H->getArea()
                *H->getEfficiencyCosine()
                *H->getTotalReflectivity()
                *H->getEfficiencyBlock()
                *H->getEfficiencyShading()
                *H->getEfficiencyCloudiness()
                *SF->getVarMap()->flux.flux_dni.val / 1000. //kW
            );
            ret->push_back(H->getEnergyValue()); //kWh -- energy delivered over the simulation time period
            ret->push_back(H->getAnnualEfficiency());
            ret->push_back(H->getEfficiencyTotal());
            ret->push_back(H->getEfficiencyCosine());
            ret->push_back(H->getEfficiencyIntercept());
            ret->push_back(H->getTotalReflectivity());
            ret->push_back(H->getEfficiencyAtten());
            ret->push_back(H->getEfficiencyBlock());
            ret->push_back(H->getEfficiencyShading());
            ret->push_back(H->getEfficiencyCloudiness());

            // TODO: can we do this?
            /*
            if (i == 0)
            {
                *nrows = helio_select.size();
                *ncols = ret->size();
                ret->reserve((*nrows) * (*ncols));
            }
            */
        }

        std::string tab_header;

        // UPDATE: If table changes
        tab_header.append("id,");
        tab_header.append("x_location,");
        tab_header.append("y_location,");
        tab_header.append("z_location,");
        tab_header.append("x_aimpoint,");
        tab_header.append("y_aimpoint,");
        tab_header.append("z_aimpoint,");
        tab_header.append("i_tracking_vector,");
        tab_header.append("j_tracking_vector,");
        tab_header.append("k_tracking_vector,");
        tab_header.append("layout_metric,");
        tab_header.append("power_to_receiver,");
        tab_header.append("power_reflected,");
        tab_header.append("energy,");
        tab_header.append("efficiency_annual,");
        tab_header.append("efficiency,");
        tab_header.append("cosine,");
        tab_header.append("intercept,");
        tab_header.append("reflectance,");
        tab_header.append("attenuation,");
        tab_header.append("blocking,");
        tab_header.append("shading,");
        tab_header.append("clouds");

        // TODO: Will this work?
        header = tab_header.c_str();

        return (sp_number_t*)ret;
    }

    return NULL;
}

SPEXPORT const char* sp_detail_results_header(sp_data_t p_data)
{
    std::string ret;



    return ret.c_str();
}

SPEXPORT sp_number_t *sp_get_fluxmap(sp_data_t p_data, int* nrows, int* ncols, int rec_id = 0)
{

    /*
	Retrieve the receiver fluxmap, optionally specifying the receiver ID to retrieve.
	Returns: ([integer:receiver id]):array
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    SolarField* SF = &mc->solarfield;

    Receiver *rec;

    if (rec_id != 0)
    {
        if (rec_id > SF->getReceivers()->size() - 1)
            return NULL;

        rec = SF->getReceivers()->at(rec_id);
    }
    else
    {
        rec = SF->getReceivers()->front();
    }

    FluxGrid *fg = rec->getFluxSurfaces()->front().getFluxMap();

    *nrows = fg->size();
    *ncols = fg->front.size();
    std::vector<double>* flux_ret;
    flux_ret->reserve((*nrows) * (*ncols));

    std::vector<std::vector<double>> flux_mat(fg->size(), std::vector<double>(fg->front().size(), 0.));

    //TODO: Ask Mike if this is suppose to be reversed, line 894 scripting.cpp
    //rows
    for (size_t i = 0; i < fg->size(); i++)
    {
        //cols
        for (size_t j = 0; j < fg->front().size(); j++)
        {
            flux_ret->push_back(fg->at(i).at(j).flux);
        }
    }
    return (sp_number_t*)flux_ret;
}

//TODO: Skipped this function initially 
SPEXPORT void sp_optimize(sp_data_t p_data, sp_number_t* pvalues, int nvar)
{
    /*
    Execute an optimization run, returning the optimized result and iteration information. 
    Variables to be optimized are passed in a vector, with each row containing a table specifying 
    {variable, step, upbound, lowbound, inital}. The table must include the variable key, others are optional. 
    
    The return table includes the following: 'result':table of variable names and associated optimized values, 
    'objective':number, 'flux':number, 'iterations':array of evaluation point, objective, flux. 
    Optional arguments include maxiterations/tolerance/defaultstep/powerpenalty/nthreads.,
    (vector:variable tables[, table:options])
    Returns: table
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    //get the variable table
    if (cxt.arg_count() < 1 || cxt.arg(0).type() != lk::vardata_t::VECTOR)
        return;

    lk::vardata_t &vartab = cxt.arg(0);

    std::stringstream heliodata;

    size_t nvars = nvar;
    for (size_t i = 0; i < nvars; i++)
    {
        for (size_t j = 0; j < 12; j++)
        {
            if (j > ncols - 1)
                heliodata << "NULL";
            else
                heliodata << pvalues[j + ncols * i];
            heliodata << (j < 11 ? "," : ";");
        }

    }
    V->sf.layout_data.val = heliodata.str();


    //set up options, if provided
    if (cxt.arg_count() == 2)
    {
        //maxiterations/tolerance/defaultstep/powerpenalty/nthreads

        lk::varhash_t *opthash = cxt.arg(1).hash();

        if (opthash->find("nthreads") != opthash->end())
            F.SetThreadCount(opthash->at("nthreads")->as_integer());

        if (opthash->find("maxiterations") != opthash->end())
            V->opt.max_iter.val = opthash->at("maxiterations")->as_integer();

        if (opthash->find("tolerance") != opthash->end())
            V->opt.converge_tol.val = opthash->at("tolerance")->as_number();

        if (opthash->find("defaultstep") != opthash->end())
            V->opt.max_step.val = opthash->at("defaultstep")->as_number();

        if (opthash->find("powerpenalty") != opthash->end())
            V->opt.power_penalty.val = opthash->at("powerpenalty")->as_number();
    }

    //set up variables
    int nv = vartab.vec()->size();
    vector<double*> optvars(nv);
    vector<double> upper(nv);
    vector<double> lower(nv);
    vector<double> stepsize(nv);
    vector<string> names(nv);

    for (size_t i = 0; i < nv; i++)
    {
        //check that the specified variable names exist
        std::string varname = vartab.vec()->at(i).hash()->at("variable")->as_string();

        if (V->_varptrs.find(varname) == V->_varptrs.end())
            throw lk::error_t("Specified variable does not exist: " + varname);

        //handle the variable
        if (V->_varptrs.at(varname)->dattype != SP_DATTYPE::SP_DOUBLE)
            throw lk::error_t("Optimized variable must be of type 'double'; discrete or boolean variables are not supported. Variable: " + varname);

        spvar<double> *varptr = static_cast<spvar<double>*>(V->_varptrs.at(varname));
        optvars.at(i) = &varptr->val;
        vector<string> namedat = split(varname, ".");
        names.at(i) = namedat.back();

        lk::varhash_t *varhash = vartab.vec()->at(i).hash();

        //bounds
        if (varhash->find("lowbound") == varhash->end())
            lower.at(i) = -HUGE_VAL;
        else
            lower.at(i) = varhash->at("lowbound")->as_number();

        if (varhash->find("upbound") == varhash->end())
            upper.at(i) = HUGE_VAL;
        else
            upper.at(i) = varhash->at("upbound")->as_number();

        if (varhash->find("initial") != varhash->end())
            varptr->val = varhash->at("initial")->as_number();

        if (varhash->find("step") == varhash->end())
            stepsize.at(i) = V->opt.max_step.val * varptr->val;
        else
            stepsize.at(i) = varhash->at("step")->as_number();

    }

    int n_threads = F.GetThreadCount();
    ArrayString *local_wfdat = F.GetLocalWeatherDataObject();
    lk::vardata_t iter_vec;
    std::vector< double > obj_vals;
    std::vector< std::vector<double> > flux_vals;
    std::vector< std::vector< double > > eval_points;

    if (n_threads > 1)
    {
        AutoPilot_MT *SFopt_MT = new AutoPilot_MT();

        SFopt_MT->SetSummaryCallback(LKInfoCallback, SF->getSimInfoObject()->getCallbackData());

        //set up the weather data for simulation
        vector<string> wdata;
        for (int i = 0; i < local_wfdat->size(); i++)
            wdata.push_back(local_wfdat->at(i));
        SFopt_MT->GenerateDesignPointSimulations(*V, wdata);

        //Do the expert setup
        SFopt_MT->Setup(*V, true);

        //run the optimization
        SFopt_MT->Optimize(optvars, upper, lower, stepsize, &names);

        //get resulting info
        SFopt_MT->GetOptimizationObject()->getOptimizationSimulationHistory(eval_points, obj_vals, flux_vals);

        try
        {
            delete SFopt_MT;
        }
        catch (...)
        {
        }
    }
    else
    {

        AutoPilot_S *SFopt_S = new AutoPilot_S();
        SFopt_S->SetSummaryCallback(LKInfoCallback, SF->getSimInfoObject()->getCallbackData());

        //set up the weather data for simulation
        vector<string> wdata;
        for (int i = 0; i < local_wfdat->size(); i++)
            wdata.push_back(local_wfdat->at(i));
        SFopt_S->GenerateDesignPointSimulations(*V, wdata);

        //Do the expert setup
        SFopt_S->Setup(*V, true);

        //run the optimization
        SFopt_S->Optimize(optvars, upper, lower, stepsize, &names);

        //get resulting info
        SFopt_S->GetOptimizationObject()->getOptimizationSimulationHistory(eval_points, obj_vals, flux_vals);


        try
        {
            delete SFopt_S;
        }
        catch (...)
        {
        }
    }


    //set up return structure
    //result/objective/flux/iterations
    cxt.result().empty_hash();

    lk::vardata_t res_hash;
    res_hash.empty_hash();

    for (size_t i = 0; i < optvars.size(); i++)
    {
        std::string varname = vartab.vec()->at(i).hash()->at("variable")->as_string();
        spvar<double> *varptr = static_cast<spvar<double>*>(V->_varptrs.at(varname));
        res_hash.hash_item(varname, varptr->val);
    }

    cxt.result().hash_item("result", res_hash);

    iter_vec.empty_vector();

    for (size_t i = 0; i < flux_vals.size(); i++)
    {
        iter_vec.vec()->push_back(lk::vardata_t());
        iter_vec.vec()->at(i).empty_vector();
        for (size_t j = 0; j < eval_points.front().size(); j++)
        {
            iter_vec.vec()->at(i).vec_append(eval_points.at(i).at(j));
        }
        iter_vec.vec()->at(i).vec_append(obj_vals.at(i));
        for (size_t j = 0; j < flux_vals.front().size(); j++)
            iter_vec.vec()->at(i).vec_append(flux_vals.at(i).at(j));
    }

    lk::vardata_t fluxresult;
    fluxresult.empty_vector();
    for (size_t j = 0; j < flux_vals.back().size(); j++)
        fluxresult.vec_append(flux_vals.back().at(j));

    cxt.result().hash_item("objective", obj_vals.back());
    cxt.result().hash_item("flux", fluxresult);
    cxt.result().hash_item("iterations", iter_vec);

}

SPEXPORT void sp_clear_land(sp_data_t p_data, const char* type = NULL)
{
    /*
	Reset the land boundary polygons, clearing any data. Optionally specify 'type' as 'inclusion' or 'exclusion'.
	Returns: ([string:type]):void
    */
    api_helper* mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    bool clear_inclusions = true;
    bool clear_exclusions = true;

    if (type != NULL)
    {
        std::string arg = lower_case(type);

        if (arg.find("inclusion") != std::string::npos)
            clear_exclusions = false;
        else if (arg.find("exclusion") != std::string::npos)
            clear_inclusions = false;
    }

    if (clear_inclusions)
        V->land.inclusions.val.clear();
    if (clear_exclusions)
        V->land.exclusions.val.clear();

    if (clear_inclusions && clear_exclusions)
        V->land.is_bounds_array.val = false;
}

SPEXPORT bool sp_add_land(sp_data_t p_data, const char* type, sp_number_t* polygon_points, int* npts , int* ndim, bool is_append = true)
{
    /*
	Add land inclusion or a land exclusion region within a specified polygon. Specify the type as 'inclusion' or 'exclusion', and optionally append (true) or overwrite (false) the existing regions.
	Returns: (array:polygon, string:type[, boolean:append=true]):boolean
    */
    api_helper* mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;
    var_land& L = V->land;

    std::string type_str = lower_case(type);

    if (type_str.find("incl") != std::string::npos)
        type_str = "inclusion";
    else if (type_str.find("excl") != std::string::npos)
        type_str = "exclusion";
    else
    {
        //invalid argument
        return false;
    }

    //convert the polygon into the required string format
    std::vector< std::string > pt, poly;
    for (size_t i = 0; i < *npts; i++)
    {
        pt.clear();
        for (size_t j = 0; j < *ndim; j++)
            pt.push_back(std::to_string(polygon_points[i*(*ndim) + j]));

        poly.push_back(join(pt, ","));
    }

    std::string spoly = "[POLY]" + join(poly, "[P]");

    if (type_str == "inclusion")
    {
        if (!is_append)
            L.inclusions.val.clear();

        L.inclusions.set_from_string(spoly.c_str());
    }
    else
    {
        if (!is_append)
            L.exclusions.val.clear();

        L.exclusions.set_from_string(spoly.c_str());
    }

    L.is_bounds_array.val = true;

    return true;

}

SPEXPORT sp_number_t *sp_heliostats_by_region(sp_data_t p_data, int* lenret, const char* coor_sys, bool is_returnloc = false,
                                                    sp_number_t* arguments = NULL, int* len_arg = NULL, 
                                                     const char* svgfname_data = NULL, sp_number_t* svg_opt_tab = NULL)
{
    /*
    Returns heliostats that fall within a region. Options are:
    >> all (no additional arguments),
    >> cylindrical (provide [rmin,rmax,azmin,azmax radians]),
    >> cartesian (provide [xmin, xmax, ymin, ymax[, zmin, zmax]]),
    >> polygon (provide [[x1,y1],[x2,y2],...]),
    >> svg (provide string with 'scale-x scale-y;offset-x offset-y;<svg path 1>;<svg path 2>;...',
    >> svgfile (provide string filename, optional table {'offset'=array, 'scale'=array}).
    (string:system, variant:region info[, string:return info - id/location])

    Returns an array of included heliostat ID's or locations. : array
    */

    api_helper* mc = static_cast<api_helper*>(p_data);
    SolarField* SF = &mc->solarfield;
    Hvector *helios = SF->getHeliostats();

    //which coordinate system?
    std::string system = (std::string) coor_sys;

    //return vector -> length is unknown a priori
    std::vector<double>* ret;

    if (helios->size() < 1)
        return NULL;

    double delta = 0.0001;

    if (lower_case(system) == "all")
    {
        for (size_t i = 0; i < helios->size(); i++)
        {
            if (is_returnloc)
            {
                ret->push_back(helios->at(i)->getLocation()->x);
                ret->push_back(helios->at(i)->getLocation()->y);
                ret->push_back(helios->at(i)->getLocation()->z);
            }
            else
                ret->push_back((double)helios->at(i)->getId());
        }
    }
    else if (lower_case(system) == "cylindrical")
    {
        double rmin, rmax, azmin, azmax;
        if (*len_arg != 4)
            return NULL;
     
        rmin = arguments[0] - delta;
        rmax = arguments[1] + delta;
        azmin = arguments[2] - delta;
        azmax = arguments[3] + delta;

        for (size_t i = 0; i < helios->size(); i++)
        {
            double rpos = helios->at(i)->getRadialPos();
            double apos = helios->at(i)->getAzimuthalPos();

            if (rpos > rmin)
                if (rpos < rmax)
                    if (apos > azmin)
                        if (apos < azmax)
                            if (is_returnloc)
                            {
                                ret->push_back(helios->at(i)->getLocation()->x);
                                ret->push_back(helios->at(i)->getLocation()->y);
                                ret->push_back(helios->at(i)->getLocation()->z);
                            }
                            else
                                ret->push_back((double)helios->at(i)->getId());
        }

    }
    else if (lower_case(system) == "cartesian")
    {
        if ((*len_arg != 4) || (*len_arg != 6))
            return NULL;

        double xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = arguments[0] - delta;
        xmax = arguments[1] + delta;
        ymin = arguments[2] - delta;
        ymax = arguments[3] + delta;
        if (*len_arg == 6)
        {
            zmin = arguments[4] - delta;
            zmax = arguments[5] + delta;
        }
        else
        {
            zmin = -9e9;
            zmax = 9e9;
        }


        for (size_t i = 0; i < helios->size(); i++)
        {
            sp_point *loc = helios->at(i)->getLocation();

            if (loc->x > xmin)
                if (loc->x < xmax)
                    if (loc->y > ymin)
                        if (loc->y < ymax)
                            if (loc->z > zmin)
                                if (loc->z < zmax)
                                    if (is_returnloc)
                                    {
                                        ret->push_back(helios->at(i)->getLocation()->x);
                                        ret->push_back(helios->at(i)->getLocation()->y);
                                        ret->push_back(helios->at(i)->getLocation()->z);
                                    }
                                    else
                                        ret->push_back((double)helios->at(i)->getId());
        }
    }
    else if (lower_case(system) == "polygon")
    {
        //construct a polygon from the listed points
        std::vector< sp_point > polygon;
        for (size_t i = 0; i < *len_arg; i++)
        {
            polygon.push_back(sp_point(arguments[i], arguments[i++], 0.));
        }

        for (size_t i = 0; i < helios->size(); i++)
        {
            if (Toolbox::pointInPolygon(polygon, *helios->at(i)->getLocation()))
                if (is_returnloc)
                {
                    ret->push_back(helios->at(i)->getLocation()->x);
                    ret->push_back(helios->at(i)->getLocation()->y);
                    ret->push_back(helios->at(i)->getLocation()->z);
                }
                else
                    ret->push_back((double)helios->at(i)->getId());
        }

    }
    else if ((lower_case(system) == "svg") || (lower_case(system) == "svgfile"))
    {
        /*
        undocumented feature

        Provide string. Coordinates are space separated, points are comma separated, shapes are semicolon separated.
        First entry is x-scale y-scale, Second entry is x-offset y-offset.

        */

        std::vector< std::string > entries;
        std::vector< std::string > scale_s;
        std::vector< std::string > offset_s;


        if (lower_case(system) == "svgfile")
        {

            if (!ioutil::file_exists(svgfname_data))
                throw std::runtime_error("Invalid SVG file - not found.");

            if (svg_opt_tab != NULL)
            {
                scale_s.push_back(std::to_string(svg_opt_tab[0]));
                scale_s.push_back(std::to_string(svg_opt_tab[1]));

                offset_s.push_back(std::to_string(svg_opt_tab[2]));
                offset_s.push_back(std::to_string(svg_opt_tab[3]));
            }
            else
            {
                scale_s.push_back("1.");
                scale_s.push_back("1.");

                offset_s.push_back("0.");
                offset_s.push_back("0.");
            }

            //load the svg file and parse it as an xml document
            using namespace rapidxml;
            //Read in the file to a string
            std::string file;        //contents of the file
            std::string eol;
            ioutil::read_file((std::string) svgfname_data, file, eol);

            char *fstr = new char[file.size() + 1];
            strncpy(fstr, (const char*)file.c_str(), file.size());
            fstr[file.size()] = 0;    //Null terminator

            xml_document<> doc;
            doc.parse<0>(fstr);
            xml_node<> *top_node = doc.first_node();    //<data>

            xml_node<> *node = top_node->first_node("g");
            xml_node<> *tnode = node->first_node("g");
            if (tnode)
                node = tnode;
            node = node->first_node("path");

            //assume that this is consistent with SVG files created by inkscape.. I don't know whether other SVG creators have a consistent XML structure. 
            entries.clear();
            while (node)
            {
                entries.push_back(node->first_attribute("d")->value());
                node = node->next_sibling("path");
            }
        }
        else
        {
            //get the string data and break it up into units
            if (svgfname_data == NULL)
                std::runtime_error("svg data must be provided for the svg option.");
            std::string data = (std::string) svgfname_data;
            entries = split(data, ";");
            scale_s = split(entries.front(), " ");
            offset_s = split(entries.at(1), " ");

            entries.erase(entries.begin(), entries.begin() + 2);
        }



        //get the scale and offset vectors
        double scale_x, scale_y, offset_x, offset_y;

        to_double(scale_s.at(0), &scale_x);
        to_double(scale_s.at(1), &scale_y);
        to_double(offset_s.at(0), &offset_x);
        to_double(offset_s.at(1), &offset_y);

        //allocate the main polygons structure
        std::vector< std::vector< sp_point > > polygons;

        for (size_t i = 0; i < entries.size(); i++)
        {
            polygons.push_back(std::vector< sp_point >());
            std::vector< sp_point > *P = &polygons.back();

            Toolbox::poly_from_svg(entries.at(i), *P, true);

            for (size_t j = 0; j < P->size(); j++)
            {
                P->at(j).x = P->at(j).x * scale_x + offset_x;
                P->at(j).y = P->at(j).y * scale_y + offset_y;
            }
        }

        //check each heliostat to see if it's in any polygon
        for (size_t i = 0; i < helios->size(); i++)
        {
            for (size_t j = 0; j < polygons.size(); j++)
            {
                std::vector< sp_point > *polygon = &polygons.at(j);
                sp_point *loc = helios->at(i)->getLocation();

                if (Toolbox::pointInPolygon(*polygon, *loc))
                {
                    if (is_returnloc)
                    {
                        ret->push_back(loc->x);
                        ret->push_back(loc->y);
                        ret->push_back(loc->z);
                    }
                    else
                        ret->push_back((double)helios->at(i)->getId());

                    //if included, don't need to check other polygons
                    break;
                }
            }
        }
    }
    else
    {
        throw std::runtime_error("invalid region type specified. Expecting one of [cylindrical, cartesian, polygon]");
    }
    // pass back size of return vector
    *lenret = ret->size();
    return (sp_number_t*)ret;
}

SPEXPORT bool sp_modify_heliostats(sp_data_t p_data, sp_number_t* helio_data, int* nhel, int* ncols, const char* table_hdr)
{

    /*
    Modify attributes of a subset of heliostats in the current layout. Modifiable attributes include 
    location/aimpoint/soiling/reflectivity/enabled, 
    and are specified in the input table by variable name and array pairs. 
    The length of each variable array must match the length of the ID array. For example, 
    modify_heliostats([1,2,3], { 'location'=[[1.,10.],[2.,11.],[3.,12.]] } );, (array:heliostat IDs, table:variables)
    
    Returns: none
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    SolarField* SF = &mc->solarfield;
    
    if (SF->getHeliostats()->size() < 1)
        return false;

    //collect all relevant heliostat ID's
    std::vector< int > hids;
    std::vector< std::vector < double >> data_table;

    //pull the ID's and data from the provided array
    for (size_t i = 0; i < (*nhel); i++)
    {
        for (size_t j = 0; j < (*ncols); j++)
        {
            if (j == 0)
                hids.push_back((int)helio_data[i * (*ncols) + j]);
            else
                data_table[i].push_back((double)helio_data[i * (*ncols) + j]);
        }
    }

    //consolidate all heliostat's into a vector by ID
    unordered_map< int, Heliostat* > *hmap = SF->getHeliostatsByID();
    Hvector helios;

    for (size_t i = 0; i < hids.size(); i++)
    {
        try
        {
            helios.push_back(hmap->at(hids.at(i)));
        }
        catch (...)
        {
        }
    }

    //get the variable table header
    std::string hdr_str(table_hdr);
    std::vector<std::string> vars;
    vars = split(hdr_str, ",");

    //these are the supported options
    std::vector< std::string > attrs = {
    "location",
    "aimpoint",
    "soiling",
    "reflectivity",
    "enabled"
    };

    //for each provided option
    int var_i = 0;
    for (std::size_t j = 0; j<*ncols - 1; j++)
    {
        std::string varname = vars[var_i];
        var_i++;

        //first make sure this is a valid attribute
        if (std::find(attrs.begin(), attrs.end(), varname) == attrs.end())
            throw std::runtime_error("Invalid attribute specified: " + varname);

        if (varname == "location")
        {
            /*
            //gather data (assuming the next two columns are y and z)
            std::vector<std::vector<double>>* locvec;
            for (std::size_t i = 0; i < helios.size(); i++)
            {
                locvec->at(i).push_back(data_table[i][j]);
                locvec->at(i).push_back(data_table[i][j+1]);
                locvec->at(i).push_back(data_table[i][j+2]);
            }
            j += 2;  // advance column count
            */

            //locations need to be modified through the layout shell object
            layout_shell *layout = SF->getLayoutShellObject();

            /*
            //make sure the heliostat ID's array is the same length as the location array
            if (locvec->size() != helios.size())
                throw std::runtime_error("The number of locations provided does not match the number of heliostat ID's provided.");
            */

            //assign location(s)
            layout->clear();

            for (size_t i = 0; i < helios.size(); i++)
            {
                //update the layout object
                layout->push_back(layout_obj());
                layout_obj& lobj = layout->back();

                lobj.aim = *helios.at(i)->getAimPoint();
                lobj.cant = *helios.at(i)->getCantVector();
                lobj.focal_x = helios.at(i)->getFocalX();
                lobj.focal_y = helios.at(i)->getFocalY();
                lobj.helio_type = helios.at(i)->getMasterTemplate()->getId();

                //update location
                lobj.location.x = data_table[i][j];
                lobj.location.y = data_table[i][j+1];
                lobj.location.z = data_table[i][j+2];

                //update enabled/in layout statuses
                lobj.is_enabled = helios.at(i)->IsEnabled();
                lobj.is_in_layout = helios.at(i)->IsInLayout();

                /* TODO:
                if (locvec->at(i).size() > 2)
                    lobj.location.z = locvec->at(i).at(2);
                */
            }
            j += 2;

            SF->PrepareFieldLayout(*SF, 0, true);
        }
        else if (varname == "aimpoint")
        {
            /* TODO:
            //make sure the heliostat ID's array is the same length as the aimpoint array
            std::vector< lk::vardata_t > *aimvec = cxt.arg(1).hash()->at("aimpoint")->vec();
            if (aimvec->size() != helios.size())
                throw std::runtime_error("The number of aimpoints provided does not match the number of heliostat ID's provided.");
            */
            
            //assign aimpoint(s)
            //gather data (assuming the next two columns are j and k)
            for (size_t i = 0; i < helios.size(); i++)
            {
                double ii = data_table[i][j];
                double jj = data_table[i][j + 1];
                double kk = data_table[i][j + 2];

                helios.at(i)->setAimPoint(ii, jj, kk);
            }
            j += 2;  // advance column count

        }
        else if (varname == "soiling")
        {
            /* TODO:
            //make sure the heliostat ID's array is the same length as the soiling array
            std::vector< lk::vardata_t > *svec = cxt.arg(1).hash()->at("soiling")->vec();
            if (svec->size() != helios.size())
                throw std::runtime_error("The number of soiling values provided does not match the number of heliostat ID's provided.");
            */
            
            //assign soiling(s)
            for (size_t i = 0; i < helios.size(); i++)
                helios.at(i)->getEfficiencyObject()->soiling = data_table[i][j];

        }
        else if (varname == "reflectivity")
        {
            /*
            //make sure the heliostat ID's array is the same length as the soiling array
            std::vector< lk::vardata_t > *svec = cxt.arg(1).hash()->at("reflectivity")->vec();
            if (svec->size() != helios.size())
                throw std::runtime_error("The number of reflectivity values provided does not match the number of heliostat ID's provided.");
            */
            
            //assign reflectivities
            for (size_t i = 0; i < helios.size(); i++)
                helios.at(i)->getEfficiencyObject()->reflectivity = data_table[i][j];

        }
        else if (varname == "enabled")
        {
            /* TODO:
            //make sure the heliostat ID's array is the same length as the soiling array
            std::vector< lk::vardata_t > *svec = cxt.arg(1).hash()->at("enabled")->vec();
            if (svec->size() != helios.size())
                throw std::runtime_error("The number of 'enabled' values provided does not match the number of heliostat ID's provided.");
            */
            //Enable
            for (size_t i = 0; i < helios.size(); i++)
                if ((int) data_table[i][j] == 1)
                    helios.at(i)->IsEnabled(true);
                else
                    helios.at(i)->IsEnabled(false);
        }

    }
    return true;
}

SPEXPORT bool sp_save_from_script(sp_data_t p_data, const char* sp_fname)
{
    /*
	Save the current case as a SolarPILOT .spt file. Returns true if successful.
	Returns: (string:path):boolean
    */

    api_helper *mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }

    try
    {
        ioutil::saveXMLInputFile(fname, *V, *F.GetParametricDataObject(), *F.GetOptimizationDataObject(), F.GetVersionInfo() + " (lk script)");
        return true;
    }
    catch (...)
    {

    }
    return false;
}

SPEXPORT bool sp_open_from_script(sp_data_t p_data, const char* sp_fname, const char* error)
{
    /*
	Open a SolarPILOT .spt case file. Returns true if successful. Updates the interface.
	Returns: (string:path):boolean
    */

    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }

    api_helper *mc = static_cast<api_helper*>(p_data);

    try
    {
        F.Open(fname, true);
        return true;
    }
    catch (std::exception &e)
    {
        error = e.what();
    }
    return false;
}

SPEXPORT bool sp_dump_varmap(sp_data_t p_data, const char* sp_fname, const char* error)
{
    /*
	Dump the variable structure to a text csv file. Returns true if successful.
	Returns: (string:path):boolean
    */
    api_helper* mc = static_cast<api_helper*>(p_data);
    var_map* V = &mc->variables;

    //valid path?
    std::string fname(sp_fname);
    if (!ioutil::dir_exists(ioutil::path_only(fname).c_str()))
    {
        return false;
    }

    //var_map *V = SPFrame::Instance().GetSolarFieldObject()->getVarMap();

    try
    {
        //create a list of all variable keys
        std::vector< std::string > names;
        names.reserve(V->_varptrs.size());

        for (unordered_map< std::string, spbase* >::iterator it = V->_varptrs.begin(); it != V->_varptrs.end(); it++)
            names.push_back(it->first);

        sort(names.begin(), names.end());   //output the variables alphabetically

        std::ofstream of( fname );
        std::string sep = " ,";

        if (of.is_open())      //check whether the file is accessible
        {
            for (size_t i = 0; i < names.size(); i++)
            {
                spbase *var = V->_varptrs.at(names.at(i));

                std::string val = var->as_string();

                if (val.size() > 30)      //tuncate very long values
                {
                    val.erase(val.begin() + 20, val.end());
                    val += "... (truncated)";
                }

                //replace all commas
                std::string::size_type n = 0;
                while ((n = val.find(",", n)) != std::string::npos)
                {
                    val.replace(n, 1, "'");
                    n += 1;
                }

                std::string units = var->units;
                //replace all commas
                n = 0;
                while ((n = units.find(",", n)) != std::string::npos)
                {
                    units.replace(n, 1, "'");
                    n += 1;
                }

                of << names.at(i) << sep << val << sep << units << sep << var->short_desc << "\n";
            }

            of.close();
            return true;  //success
        }

    }
    catch (std::exception &e)
    {
        error = e.what();
    }
    return false;
}
