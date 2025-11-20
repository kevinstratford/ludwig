/*****************************************************************************
 *
 *  colloid_options.c
 *
 *  Container for options required to instantiate the colloid_info_t
 *  management object.
 *
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2025 The University of Edinburgh
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#include <assert.h>
#include <string.h>

#include "colloid_options.h"
#include "util_json.h"
#include "util.h"

static const double RHO_ZERO_DEFAULT      = 1.0; /* unit density lb units */
static const int    HAVE_COLLOIDS_DEFAULT = 0;   /* no colloids */

/*****************************************************************************
 *
 *  colloid_options_default
 *
 *****************************************************************************/

colloid_options_t colloid_options_default(void) {

  /* Default is no colloids, but other values are required. */

  colloid_options_t options = {
      .ncell  = {2, 2, 2},                   /* min. local list dimensions */
      .nvel   = 19,                          /* default bbl model */
      .rho0   = RHO_ZERO_DEFAULT,
      .input  = colloid_io_options_default(),
      .output = colloid_io_options_default(),
      .nfreq  = 0,                           /* No output */

      .have_colloids = HAVE_COLLOIDS_DEFAULT /* colloids expected or not */
  };

  return options;
}

/*****************************************************************************
 *
 *  colloid_options_have_colloids
 *
 *****************************************************************************/

colloid_options_t colloid_options_have_colloids(int have_colloids) {

  colloid_options_t options = colloid_options_default();

  options.have_colloids = have_colloids;

  return options;
}

/*****************************************************************************
 *
 *  colloid_options_ncell
 *
 *****************************************************************************/

colloid_options_t colloid_options_ncell(const int ncell[3]) {

  colloid_options_t options = colloid_options_default();

  /* Request for ncell => have_colloids */
  options.have_colloids = 1;
  options.ncell[0]      = ncell[0];
  options.ncell[1]      = ncell[1];
  options.ncell[2]      = ncell[2];

  return options;
}

/*****************************************************************************
 *
 *  colloid_options_ncell_valid
 *
 *  A minimal requirement.
 *
 *****************************************************************************/

int colloid_options_ncell_valid(const int ncell[3]) {

  int isvalid = (ncell[0] > 1 && ncell[1] > 1 && ncell[2] > 1);

  return isvalid;
}

/*****************************************************************************
 *
 *  colloid_options_valid
 *
 *****************************************************************************/

int colloid_options_valid(const colloid_options_t * options) {

  int isvalid = 1;

  assert(options);

  if (0 == colloid_options_ncell_valid(options->ncell)) isvalid = 0;
  if (options->rho0 <= 0.0)                             isvalid = 0;
  if (0 == colloid_io_options_valid(&options->input))   isvalid = 0;
  if (0 == colloid_io_options_valid(&options->output))  isvalid = 0;
  if (options->nfreq < 0)                               isvalid = 0;

  return isvalid;
}

/*****************************************************************************
 *
 *  colloid_options_to_json
 *
 *****************************************************************************/

int colloid_options_to_json(const colloid_options_t * opts, cJSON ** json) {

  int ifail = 0;

  if (json == NULL || opts == NULL) {
    ifail = -1;
  }
  else {
    /* Create new JSON object ... */
    cJSON * myjson = cJSON_CreateObject();
    cJSON * ncell  = cJSON_CreateIntArray(opts->ncell, 3);

    cJSON_AddBoolToObject(myjson, "have_colloids", opts->have_colloids);
    cJSON_AddItemToObject(myjson, "ncell", ncell);
    cJSON_AddNumberToObject(myjson, "nvel", opts->nvel);
    cJSON_AddNumberToObject(myjson, "rho0", opts->rho0);

    {
      cJSON * input  = NULL;
      cJSON * output = NULL;

      /* input */
      ifail = colloid_io_options_to_json(&opts->input, &input);
      if (ifail == 0) cJSON_AddItemToObject(myjson, "Input options", input);

      /* output */
      ifail = colloid_io_options_to_json(&opts->output, &output);
      if (ifail == 0) cJSON_AddItemToObject(myjson, "Output options", output);
    }

    cJSON_AddNumberToObject(myjson, "nfreq", opts->nfreq);

    *json = myjson;
  }

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_options_from_json
 *
 *  FIXME: decide what is mandatory, and what is optional.
 *
 *****************************************************************************/

int colloid_options_from_json(const cJSON * json, colloid_options_t * opts) {

  int ifail = 0;

  if (opts == NULL) {
    ifail = -1;
  }
  else {
    /* all elements currently required to be present */
    cJSON * have   = cJSON_GetObjectItemCaseSensitive(json, "have_colloids");
    cJSON * ncell  = cJSON_GetObjectItemCaseSensitive(json, "ncell");
    cJSON * nvel   = cJSON_GetObjectItemCaseSensitive(json, "nvel");
    cJSON * rho0   = cJSON_GetObjectItemCaseSensitive(json, "rho0");
    cJSON * input  = cJSON_GetObjectItemCaseSensitive(json, "Input options");
    cJSON * output = cJSON_GetObjectItemCaseSensitive(json, "Output options");
    cJSON * nfreq  = cJSON_GetObjectItemCaseSensitive(json, "nfreq");

    if (0 == cJSON_IsBool(have))                                  ifail +=  1;
    if (3 != util_json_to_int_array(ncell, opts->ncell, 3))       ifail +=  2;
    if (0 == cJSON_IsNumber(nvel))                                ifail +=  4;
    if (0 == cJSON_IsNumber(rho0))                                ifail +=  8;

    if (have) opts->have_colloids = cJSON_IsTrue(have);
    if (nvel) opts->nvel = cJSON_GetNumberValue(nvel);
    if (rho0) opts->rho0 = cJSON_GetNumberValue(rho0);

    if (0 != colloid_io_options_from_json(input, &opts->input))   ifail += 16;
    if (0 != colloid_io_options_from_json(output, &opts->output)) ifail += 32;
    if (0 == cJSON_IsNumber(nfreq))                               ifail += 64;

    if (nfreq) opts->nfreq = cJSON_GetNumberValue(nfreq);
  }

  return ifail;
}
