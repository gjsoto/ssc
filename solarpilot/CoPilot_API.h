#ifndef _COPILOT_API_
#define _COPILOT_API_

#include <unordered_map>
#include <string>
#include "SolarField.h"

#if defined(_WINDLL)
#define SPEXPORT __declspec(dllexport)
#else
#define SPEXPORT
#endif
#ifdef __cplusplus
extern "C" {
#endif

    /** An opaque reference to a structure that holds a collection of variables.  This structure can contain any number of variables referenced by name, and can hold strings, numbers, arrays, and matrices.  Matrices are stored in row-major order, where the array size is nrows*ncols, and the array index is calculated by r*ncols+c. An ssc_data_t object holds all input and output variables for a simulation. It does not distinguish between input, output, and input variables - that is handled at the model context level. */
    typedef void* sp_data_t;

    /** The numeric type used in the SolarPILOT API. */ 
    typedef float sp_number_t;

    /** Creates a new data object in memory.  A data object stores a table of named values, where each value can be of any SolarPILOT datatype. */
    SPEXPORT sp_data_t sp_data_create();

    /** Frees the memory associated with a data object, where p_data is the data container to free. */
    SPEXPORT void sp_data_free(sp_data_t p_data);

    SPEXPORT void sp_set_value(sp_data_t p_data, const char* name, sp_number_t v);

    SPEXPORT sp_number_t sp_get_value(sp_data_t p_data, const char* name);




#ifdef __cplusplus
}
#endif


#endif  // _COPILOT_API_
