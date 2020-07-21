/*****************************************************************************
 *
 *  runtime.h
 *
 *  Runtime input interface.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2010-2020 The University of Edinburgh
 *
 *  Contributing authors:
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef LUDWIG_RUNTIME_H
#define LUDWIG_RUNTIME_H

#include "pe.h"

typedef struct rt_s rt_t;

int rt_create(pe_t * pe, rt_t ** prt);
int rt_free(rt_t * rt);
int rt_read_input_file(rt_t * rt, const char * filename);
int rt_info(rt_t * rt);
int rt_int_parameter(rt_t * rt, const char * key, int * ivalue);
int rt_int_parameter_vector(rt_t * rt, const char * key, int ivalue[3]);
int rt_double_parameter(rt_t * rt, const char * key, double * value);
int rt_double_parameter_vector(rt_t * rt, const char * key, double value[3]);
int rt_string_parameter(rt_t * rt, const char * key, char * s, unsigned  int len);
int rt_switch(rt_t * rt, const char * key);
int rt_active_keys(rt_t * rt, int * nactive);
int rt_add_key_value(rt_t * rt, const char * key, const char * value);

#endif
