#include "CoPilot_API.h"


struct api_helper
{
    SolarField solarfield;
    var_map variables;

    api_helper()
    {
        variables.reset();
        solarfield.Create(variables);
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
    if (mc->variables._varptrs.find(name) != mc->variables._varptrs.end())
    {
        mc->variables._varptrs[name]->set_from_string(to_string(v).c_str());
    }
}

SPEXPORT sp_number_t sp_get_value(sp_data_t p_data, const char* name)
{
    api_helper *mc = static_cast<api_helper*>(p_data);
    return (sp_number_t)mc->get(name);
}


SPEXPORT int sp_test(int v)
{
    return v + 10;
}



