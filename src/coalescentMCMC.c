/* coalescentMCMC.c (2012-10-01) */

/* Copyright 2012 Emmanuel Paradis

/* This file is part of the R-package `coalescentMCMC'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

void get_single_index_integer(int *x, int *val, int *index)
{
	while (x[*index] != *val) (*index)++;
	*index += 1;
}

void get_two_index_integer(int *x, int *val, int *index)
{
	while (x[index[0]] != *val) index[0]++;
	index[1] = index[0] + 1;
	while (x[index[1]] != *val) index[1]++;
	index[0] += 1;
	index[1] += 1;
}
