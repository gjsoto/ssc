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

    SPEXPORT void sp_cancel_simulation(sp_data_t p_data);

    /** Creates a new data object in memory.  A data object stores a table of named values, where each value can be of any SolarPILOT datatype. */
    SPEXPORT sp_data_t sp_data_create();

    /** Frees the memory associated with a data object, where p_data is the data container to free. */
    SPEXPORT void sp_data_free(sp_data_t p_data);

    SPEXPORT void var_free_memory(sp_number_t* varptr);

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

    SPEXPORT void sp_get_matrix(sp_data_t p_data, const char* name, sp_number_t* values, int* ncols, int* nrows);

    SPEXPORT void sp_reset_geometry(sp_data_t p_data);

    SPEXPORT int sp_add_receiver(sp_data_t p_data, const char* receiver_name);

    SPEXPORT int sp_drop_receiver(sp_data_t p_data, const char* receiver_name);

    SPEXPORT int sp_add_heliostat_template(sp_data_t p_data, const char* heliostat_name);

    SPEXPORT int sp_drop_heliostat_template(sp_data_t p_data, const char* heliostat_name);

    SPEXPORT int sp_update_geometry(sp_data_t p_data);

    SPEXPORT bool sp_generate_layout(sp_data_t p_data, int nthreads);

    SPEXPORT bool sp_assign_layout(sp_data_t p_data, sp_number_t* pvalues, int nrows, int ncols, int nthreads);

    SPEXPORT bool sp_get_layout_info(sp_data_t p_data, sp_number_t* layoutinfo, int* nhelio, int* ncol);

    SPEXPORT bool sp_simulate(sp_data_t p_data, int nthreads, bool save_detail, bool update_aimpoints);

    SPEXPORT const char* sp_summary_results(sp_data_t p_data);

    SPEXPORT bool sp_detail_results(sp_data_t p_data, sp_number_t* ret, int* nrows, int* ncols, const char* header, sp_number_t* selhel, int nselhel);

    SPEXPORT bool sp_get_fluxmap(sp_data_t p_data, sp_number_t* fluxmap, int* nrows, int* ncols, int rec_id);

    SPEXPORT void sp_optimize(sp_data_t p_data, sp_number_t* pvalues, int nvar);

    SPEXPORT void sp_clear_land(sp_data_t p_data, const char* type);

    SPEXPORT bool sp_add_land(sp_data_t p_data, const char* type, sp_number_t* polygon_points, int* npts, int* ndim, bool is_append);

    SPEXPORT bool sp_heliostats_by_region(sp_data_t p_data, sp_number_t* retvec, int* lenret, const char* coor_sys,
                                            sp_number_t* arguments, int* len_arg, const char* svgfname_data, sp_number_t* svg_opt_tab);

    SPEXPORT bool sp_modify_heliostats(sp_data_t p_data, sp_number_t* helio_data, int* nhel, int* ncols, const char* table_hdr);

    SPEXPORT bool sp_save_from_script(sp_data_t p_data, const char* sp_fname);

    SPEXPORT bool sp_dump_varmap(sp_data_t p_data, const char* sp_fname);

#ifdef __cplusplus
}
#endif


#endif  // _COPILOT_API_

SPEXPORT int sp_update_geometry(sp_data_t p_data);
