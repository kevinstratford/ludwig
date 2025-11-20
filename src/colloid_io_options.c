/*****************************************************************************
 *
 *  colloid_io_options.c
 *
 *  Container for options related to colloid i/o (input, output).
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

#include "colloid_io_options.h"
#include "util_json.h"
#include "util.h"

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

  char                   value[BUFSIZ] = {0};
  colloid_io_mode_enum_t mode          = COLLOID_IO_MODE_INVALID;

  strncpy(value, str, BUFSIZ - 1);
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

  if (0 == colloid_io_mode_valid(opts->mode))               valid = 0;
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
      char * str      = cJSON_GetStringValue(format);
      opts->iorformat = io_record_format_from_string(str);
    }

    if (report) {
      opts->report = cJSON_IsTrue(report);
    }

    if (iogrid) {
      if (3 != util_json_to_int_array(iogrid, opts->iogrid, 3)) ifail += 16;
    }
  }

  return ifail;
}
