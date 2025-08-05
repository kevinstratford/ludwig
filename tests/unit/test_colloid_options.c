/*****************************************************************************
 *
 *  test_colloid_options.c
 *
 *  FIMME: io_options and options
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "colloid_options.h"

int test_colloid_io_mode_valid(void);
int test_colloid_io_mode_to_string(void);
int test_colloid_io_mode_from_string(void);

int test_colloid_io_options_default(void);
int test_colloid_io_options_valid(void);
int test_colloid_io_options_to_json(void);
int test_colloid_io_options_from_json(void);

int test_colloid_options_default(void);
int test_colloid_options_ncell_valid(void);
int test_colloid_options_valid(void);

int test_colloid_options_to_json(void);
int test_colloid_options_from_json(void);

/*****************************************************************************
 *
 *  test_colloid_options_suite
 *
 *****************************************************************************/

int test_colloid_options_suite(void) {

  pe_t * pe = NULL;

  pe_create(MPI_COMM_WORLD, PE_QUIET, &pe);

  test_colloid_io_mode_valid();
  test_colloid_io_mode_to_string();
  test_colloid_io_mode_from_string();

  test_colloid_io_options_default();
  test_colloid_io_options_valid();
  test_colloid_io_options_to_json();
  test_colloid_io_options_from_json();
  
  test_colloid_options_default();
  test_colloid_options_ncell_valid();
  test_colloid_options_valid();

  test_colloid_options_to_json();
  test_colloid_options_from_json();

  pe_info(pe, "%-9s %s\n", "PASS", __FILE__);
  pe_free(pe);

  return 0;
}

/*****************************************************************************
 *
 *  test_colloid_io_mode_valid
 *
 *****************************************************************************/

int test_colloid_io_mode_valid(void) {

  int ifail = 0;

  /* valid */
  {
    colloid_io_mode_enum_t mode = COLLOID_IO_MODE_ANSI;
    ifail = colloid_io_mode_valid(mode);
    assert(ifail != 0);
  }

  /* invalid */
  {
    colloid_io_mode_enum_t mode = COLLOID_IO_MODE_INVALID;
    ifail = colloid_io_mode_valid(mode);
    assert(ifail == 0);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_mode_to_string
 *
 *****************************************************************************/

int test_colloid_io_mode_to_string(void) {

  int ifail = 0;

  /* invalid */
  {
    colloid_io_mode_enum_t mode = COLLOID_IO_MODE_INVALID;
    const char * str = colloid_io_mode_to_string(mode);
    if (0 != strncmp("invalid", str, strlen(str))) ifail = -1;
    assert(ifail == 0);
  }

  /* ansi */
  {
    colloid_io_mode_enum_t mode = COLLOID_IO_MODE_ANSI;
    const char * str = colloid_io_mode_to_string(mode);
    if (0 != strncmp("ansi", str, strlen(str))) ifail = -1;
    assert(ifail == 0);
  }

  /* mpiio */
  {
    colloid_io_mode_enum_t mode = COLLOID_IO_MODE_MPIIO;
    const char * str = colloid_io_mode_to_string(mode);
    if (0 != strncmp("mpiio", str, strlen(str))) ifail = -1;
    assert(ifail == 0);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_mode_from_string
 *
 *****************************************************************************/

int test_colloid_io_mode_from_string(void) {

  int ifail = 0;

  /* invalid */
  {
    colloid_io_mode_enum_t mode = colloid_io_mode_from_string("rubiish");
    if (mode != COLLOID_IO_MODE_INVALID) ifail = -1;
    assert(ifail == 0);
  }

  /* ansi */
  {
    colloid_io_mode_enum_t mode = colloid_io_mode_from_string("ansi");
    if (mode != COLLOID_IO_MODE_ANSI) ifail = -1;
    assert(ifail == 0);
  }

  /* mpiio */
  {
    colloid_io_mode_enum_t mode = colloid_io_mode_from_string("MPIIO");
    if (mode != COLLOID_IO_MODE_MPIIO) ifail = -1;
    assert(ifail == 0);
  }
  
  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_options_default
 *
 *****************************************************************************/

int test_colloid_io_options_default(void) {

  colloid_io_options_t opts = colloid_io_options_default();

  assert(opts.mode      == COLLOID_IO_MODE_ANSI);
  assert(opts.iorformat == IO_RECORD_ASCII);
  assert(opts.report    == 0);

  assert(opts.iogrid[0] == 1);
  assert(opts.iogrid[1] == 1);
  assert(opts.iogrid[2] == 1);

  return 0;
}

/*****************************************************************************
 *
 *  test_colloid_io_options_valid
 *
 *****************************************************************************/

int test_colloid_io_options_valid(void) {

  int ifail = 0;

  /* valid */
  {
    colloid_io_options_t opts = colloid_io_options_default();
    ifail = colloid_io_options_valid(&opts);
    assert(ifail);
  }

  /* invalid */
  {
    colloid_io_options_t opts = {0};
    ifail = colloid_io_options_valid(&opts);
    assert(ifail == 0);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_options_to_json
 *
 *****************************************************************************/

int test_colloid_io_options_to_json(void) {

  int ifail = 0;

  {
    colloid_io_options_t opts = colloid_io_options_default();
    cJSON * json = NULL;
    ifail = colloid_io_options_to_json(&opts, &json);
    assert(ifail == 0);
    /* mode */
    {
      cJSON * key = cJSON_GetObjectItemCaseSensitive(json, "I/O mode");
      char  * str = cJSON_GetStringValue(key);
      colloid_io_mode_enum_t mode = colloid_io_mode_from_string(str);
      if (mode != COLLOID_IO_MODE_ANSI) ifail = -1;
      assert(ifail == 0);
    }
    /* ior format */
    {
      cJSON * key = cJSON_GetObjectItemCaseSensitive(json, "Format");
      char * str  = cJSON_GetStringValue(key);
      io_record_format_enum_t ior = io_record_format_from_string(str);
      if (ior != IO_RECORD_ASCII) ifail = -1;
      assert(ifail == 0);
    }
    /* report */
    {
      int report = -1;
      cJSON * key = cJSON_GetObjectItemCaseSensitive(json, "Report");
      if (key) report = cJSON_IsTrue(key);
      if (report != 0) ifail = -1;
    }
    /* iogrid (check key exsits only) */
    {
      cJSON * key = cJSON_GetObjectItemCaseSensitive(json, "I/O grid");
      if (key == NULL) ifail = -1;
      assert(ifail == 0);
    }
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_io_options_from_json
 *
 *****************************************************************************/

int test_colloid_io_options_from_json(void) {

  int ifail = 0;

  /* Form JSON from a string */
  {
    const char * str =
      "{\"I/O mode\": \"mpiio\","
      "\"Format\": \"binary\","
      "\"Report\": true,"
      "\"I/O grid\": [2,3,4] }";
    cJSON * json = cJSON_Parse(str);
    assert(json);

    colloid_io_options_t opts = {0};
    ifail = colloid_io_options_from_json(json, &opts);
    assert(ifail == 0);
    assert(opts.mode == COLLOID_IO_MODE_MPIIO);
    assert(opts.iorformat == IO_RECORD_BINARY);
    assert(opts.report);
    assert(opts.iogrid[0] == 2);
    assert(opts.iogrid[1] == 3);
    assert(opts.iogrid[2] == 4);
  }

  return 0;
}

/*****************************************************************************
 *
 *  test_colloid_options_default
 *
 *****************************************************************************/

int test_colloid_options_default(void) {

  int ifail = 0;
  colloid_options_t options = colloid_options_default();

  assert(options.ncell[0] == 2);
  assert(options.ncell[1] == 2);
  assert(options.ncell[2] == 2);
  assert(options.nvel     == 19);

  if (fabs(options.rho0 - 1.0) > DBL_EPSILON) ifail = -1;
  assert(ifail == 0);
  
  assert( colloid_io_options_valid(&options.input) );
  assert( colloid_io_options_valid(&options.output) );
  assert(options.nfreq     == 0);

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_options_ncell_valid
 *
 *****************************************************************************/

int test_colloid_options_ncell_valid(void) {

  int ifail= 0;

  /* valid */
  {
    int ncell[3] = {3, 4, 5};
    ifail = colloid_options_ncell_valid(ncell);
    assert(ifail);
  }

  /* invalid */
  {
    int ncell[3] = {2, 2, 1};
    ifail = colloid_options_ncell_valid(ncell);
    assert(ifail == 0);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_options_valid
 *
 *****************************************************************************/

int test_colloid_options_valid(void) {

  int ifail = 0;

  /* valid */
  {
    colloid_options_t options = colloid_options_default();
    ifail = colloid_options_valid(&options);
    assert(ifail);
  }

  /* invalid */
  {
    colloid_options_t options = {.rho0 = -1.0};
    ifail = colloid_options_valid(&options);
    assert(ifail == 0);
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_options_to_json
 *
 *****************************************************************************/

int test_colloid_options_to_json(void) {

  int ifail = 0;

  {
    /* Default */
    colloid_options_t opts = colloid_options_default();
    cJSON * json = NULL;

    ifail = colloid_options_to_json(&opts, &json);
    assert(ifail == 0);

    if (json) {
      cJSON * jnvel = cJSON_GetObjectItemCaseSensitive(json, "nvel");
      cJSON * jrho0 = cJSON_GetObjectItemCaseSensitive(json, "rho0");

      assert(cJSON_IsNumber(jnvel));
      assert(cJSON_IsNumber(jrho0));

      /* Test a subset of components */
      if (jnvel) {
	int nvel = cJSON_GetNumberValue(jnvel);
	if (nvel != opts.nvel) ifail = -1;
	assert(ifail == 0);
      }

      if (jrho0) {
	double rho0 = cJSON_GetNumberValue(jrho0);
	if (fabs(rho0 - opts.rho0) > DBL_EPSILON) ifail = -1;
	assert(ifail == 0);
      }

      cJSON_Delete(json);
    }
  }

  return ifail;
}

/*****************************************************************************
 *
 *  test_colloid_options_from_json
 *
 *****************************************************************************/

int test_colloid_options_from_json(void) {

  int ifail = 0;

  {
    /* All components are required. Slightly long-winded ... */
    const char * str =
      "{\"ncell\": [3, 4, 5], "
      "\"nvel\": 27, "
      "\"rho0\": 2.00, "
      "\"Input options\": { "
        "\"I/O mode\": \"ansi\", "
        "\"Format\":   \"binary\", "
        "\"Report\": false, "
        "\"I/O grid\": [6, 7, 8] }, "
      "\"Output options\": { "
        "\"I/O mode\": \"mpiio\", "
        "\"Format\": \"ascii\", "
        "\"Report\": true, "
        "\"I/O grid\": [0, 1, 2] }, "
      "\"nfreq\": 1000 }";

    cJSON * json = cJSON_Parse(str);
    assert(json);

    if (json) {
      colloid_options_t opts = {0};
      ifail = colloid_options_from_json(json, &opts);
      assert(ifail == 0);
      assert(opts.ncell[0] == 3);
      assert(opts.ncell[1] == 4);
      assert(opts.ncell[2] == 5);
      assert(opts.nvel     == 27);
      assert(fabs(opts.rho0 - 2.0) < DBL_EPSILON);
      /* Input */
      assert(opts.input.mode      == COLLOID_IO_MODE_ANSI);
      assert(opts.input.iorformat == IO_RECORD_BINARY);
      assert(opts.input.report    == 0);
      assert(opts.input.iogrid[0] == 6);
      assert(opts.input.iogrid[1] == 7);
      assert(opts.input.iogrid[2] == 8);
      /* Output */
      assert(opts.output.mode      == COLLOID_IO_MODE_MPIIO);
      assert(opts.output.iorformat == IO_RECORD_ASCII);
      assert(opts.output.report    != 0);
      assert(opts.output.iogrid[0] == 0);
      assert(opts.output.iogrid[1] == 1);
      assert(opts.output.iogrid[2] == 2);
      /* Frequency */
      assert(opts.nfreq == 1000);
    }

    cJSON_Delete(json);
  }

  return 0;
}
