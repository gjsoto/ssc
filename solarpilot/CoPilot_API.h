#ifndef _COPILOT_API_
#define _COPILOT_API_



#include <unordered_map>
#include <string>
#include <ctime>
#include "SolarField.h"
//#include "interop.h"

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

    SPEXPORT void sp_set_string(sp_data_t p_data, const char *name, const char *value);

    /** Assigns value of type @a SSC_ARRAY */
    SPEXPORT void sp_set_array(sp_data_t p_data, const char *name, sp_number_t *pvalues, int length);

    /** Assigns value of type @a SSC_MATRIX . Matrices are specified as a continuous array, in row-major order.  Example: the matrix [[5,2,3],[9,1,4]] is stored as [5,2,3,9,1,4]. */
    SPEXPORT void sp_set_matrix(sp_data_t p_data, const char *name, sp_number_t *pvalues, int nrows, int ncols);




    SPEXPORT sp_number_t sp_get_number(sp_data_t p_data, const char* name);

    /** Returns the value of a @a SSC_STRING variable with the given name. */
    SPEXPORT const char *ssc_data_get_string(sp_data_t p_data, const char *name);

    /** Returns the value of a @a SSC_ARRAY variable with the given name. */
    SPEXPORT void sp_get_array(sp_data_t p_data, const char *name, sp_number_t* values, int *length);




#ifdef __cplusplus
}
#endif


#endif  // _COPILOT_API_

SPEXPORT int sp_update_geometry(sp_data_t p_data);
