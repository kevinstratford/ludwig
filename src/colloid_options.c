/*****************************************************************************
 *
 *  colloid_options.c
 *
 *  Container for options required to instantiate the colloid_info_t
 *  management object.
 *
 *****************************************************************************/

#include <assert.h>
#include <string.h>

#include "colloid_options.h"
#include "util_json.h"
#include "util.h"

static const double RHO_ZERO_DEFAULT = 1.0;

/*****************************************************************************
 *
 *  colloid_io_mode_valid
 *
 *****************************************************************************/

int colloid_io_mode_valid(colloid_io_mode_enum_t mode) {

  int valid = 0;

  switch (mode) {
  case COLLOID_IO_MODE_ANSI:
  case COLLOID_IO_MODE_MPIIO:
    valid = 1;
    break;
  default:
    valid = 0;
  }

  return valid;
}

/*****************************************************************************
 *
 *  colloid_io_mode_to_string
 *
 *****************************************************************************/

const char * colloid_io_mode_to_string(colloid_io_mode_enum_t mode) {

  const char * str = NULL;

  switch (mode) {
  case COLLOID_IO_MODE_ANSI:
    str = "ansi";
    break;
  case COLLOID_IO_MODE_MPIIO:
    str = "mpiio";
    break;
  default:
    str = "invalid";
  }

  return str;
}

/*****************************************************************************
 *
 *  colloid_io_mode_from_string
 *
 *****************************************************************************/

colloid_io_mode_enum_t colloid_io_mode_from_string(const char * str) {

  colloid_io_mode_enum_t mode = COLLOID_IO_MODE_INVALID;
  char value[BUFSIZ] = {0};

  strncpy(value, str, BUFSIZ-1);
  util_str_tolower(value, strlen(value));

  if (strcmp(value, "ansi")  == 0) mode = COLLOID_IO_MODE_ANSI;
  if (strcmp(value, "mpiio") == 0) mode = COLLOID_IO_MODE_MPIIO; 

  return mode;
}

/*****************************************************************************
 *
 *  colloid_io_options_default
 *
 *****************************************************************************/

colloid_io_options_t colloid_io_options_default(void) {

  colloid_io_options_t opts = {
    .mode      = COLLOID_IO_MODE_ANSI,
    .iorformat = IO_RECORD_ASCII,
    .report    = 0,
    .iogrid    = {1, 1, 1}
  };

  return opts;
}

/*****************************************************************************
 *
 *  colloid_io_options_valid
 *
 *****************************************************************************/

int colloid_io_options_valid(const colloid_io_options_t * opts) {

  int valid = 1;

  if (0 == colloid_io_mode_valid(opts->mode)) valid = 0;
  if (0 == io_options_record_format_valid(opts->iorformat)) valid = 0;

  return valid;
}

/*****************************************************************************
 *
 *  colloid_io_options_to_json
 *
 *****************************************************************************/

int colloid_io_options_to_json(const colloid_io_options_t * opts,
			       cJSON ** json) {
  int ifail = 0;

  if (opts == NULL || json == NULL) {
    ifail = -1;
  }
  else {
    cJSON * myjson = cJSON_CreateObject();
    cJSON * iogrid = cJSON_CreateIntArray(opts->iogrid, 3);

    cJSON_AddStringToObject(myjson, "I/O mode",
			    colloid_io_mode_to_string(opts->mode));
    cJSON_AddStringToObject(myjson, "Format",
			    io_record_format_to_string(opts->iorformat));
    cJSON_AddBoolToObject(myjson, "Report", opts->report);
    cJSON_AddItemToObject(myjson, "I/O grid", iogrid);

    *json = myjson;
  }

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_io_options_from_json
 *
 *****************************************************************************/

int colloid_io_options_from_json(const cJSON * json,
				 colloid_io_options_t * opts) {
  int ifail = 0;

  if (json == NULL || opts == NULL) {
    ifail = -1;
  }
  else {
    /* We expect four key/value pairs */
    cJSON * mode   = cJSON_GetObjectItemCaseSensitive(json, "I/O mode");
    cJSON * format = cJSON_GetObjectItemCaseSensitive(json, "Format");
    cJSON * report = cJSON_GetObjectItemCaseSensitive(json, "Report");
    cJSON * iogrid = cJSON_GetObjectItemCaseSensitive(json, "I/O grid");

    if (mode   == NULL) ifail += 1;
    if (format == NULL) ifail += 2;
    if (report == NULL) ifail += 4;
    if (iogrid == NULL) ifail += 8;

    if (mode) {
      char * str = cJSON_GetStringValue(mode);
      opts->mode = colloid_io_mode_from_string(str);
    }
    if (format) {
      char * str = cJSON_GetStringValue(format);
      opts->iorformat = io_record_format_from_string(str);
    }

    if (report) opts->report = cJSON_IsTrue(report);

    if (iogrid) {
      if (3 != util_json_to_int_array(iogrid, opts->iogrid, 3)) ifail += 16;
    }
  }

  return ifail;
}

/*****************************************************************************
 *
 *  colloid_options_default
 *
 *****************************************************************************/

colloid_options_t colloid_options_default(void) {
  
  colloid_options_t options = {
    .ncell  = {2, 2, 2},                      /* min. local list dimensions */
    .nvel   = 19,                             /* default bbl model */
    .rho0   = RHO_ZERO_DEFAULT,
    .input  = colloid_io_options_default(),
    .output = colloid_io_options_default(),
    .nfreq  = 0                               /* No output */
  };

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
    cJSON * myjson = cJSON_CreateObject();
    cJSON * ncell  = cJSON_CreateIntArray(opts->ncell, 3);

    cJSON_AddItemToObject(myjson,   "ncell", ncell);
    cJSON_AddNumberToObject(myjson, "nvel", opts->nvel);
    cJSON_AddNumberToObject(myjson, "rho0", opts->rho0);

    {
      cJSON * input = NULL;
      cJSON * output = NULL;
      ifail = colloid_io_options_to_json(&opts->input, &input);
      if (ifail == 0) cJSON_AddItemToObject(myjson, "Input options", input);

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
 *  FIXME: all mandatory. optional elements wanted?
 *
 *****************************************************************************/

int colloid_options_from_json(const cJSON * json, colloid_options_t * opts) {

  int ifail = 0;

  if (opts == NULL) {
    ifail = -1;
  }
  else {
    /* all elements currently required to be present */
    cJSON * ncell  = cJSON_GetObjectItemCaseSensitive(json, "ncell");
    cJSON * nvel   = cJSON_GetObjectItemCaseSensitive(json, "nvel");
    cJSON * rho0   = cJSON_GetObjectItemCaseSensitive(json, "rho0");
    cJSON * input  = cJSON_GetObjectItemCaseSensitive(json, "Input options");
    cJSON * output = cJSON_GetObjectItemCaseSensitive(json, "Output options");
    cJSON * nfreq  = cJSON_GetObjectItemCaseSensitive(json, "nfreq");

    if (3 != util_json_to_int_array(ncell, opts->ncell, 3))       ifail +=  1;
    if (0 == cJSON_IsNumber(nvel))                                ifail +=  2;
    if (0 == cJSON_IsNumber(rho0))                                ifail +=  4;

    if (nvel) opts->nvel = cJSON_GetNumberValue(nvel);
    if (rho0) opts->rho0 = cJSON_GetNumberValue(rho0);

    if (0 != colloid_io_options_from_json(input, &opts->input))   ifail +=  8;
    if (0 != colloid_io_options_from_json(output, &opts->output)) ifail += 16;
    if (0 == cJSON_IsNumber(nfreq))                               ifail += 32;

    if (nfreq) opts->nfreq = cJSON_GetNumberValue(nfreq);
  }

  return ifail;
}
