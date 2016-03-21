//FastQTL: Fast and efficient QTL mapper for molecular phenotypes
//Copyright (C) 2015 Olivier DELANEAU, Alfonso BUIL, Emmanouil DERMITZAKIS & Olivier DELANEAU
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "data.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#define BETA_SHAPE1_MIN 0.1
#define BETA_SHAPE2_MIN 1
#define BETA_SHAPE1_MAX 10
#define BETA_SHAPE2_MAX 1000000		//to be changed for trans!

double betaLogLikelihood(const gsl_vector *v, void *params) {
	double * p = (double *) params;
	double beta_shape1 = gsl_vector_get(v, 0);
	double beta_shape2 = gsl_vector_get(v, 1);
	return -1.0 * ((beta_shape1 - 1) * p[0] + (beta_shape2 - 1) * p[1] - p[2] * gsl_sf_lnbeta(beta_shape1, beta_shape2));
}

int data::mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2) {

	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, beta_shape1);
	gsl_vector_set (x, 1, beta_shape2);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (2);
	gsl_vector_set (ss, 0, beta_shape1/10);
	gsl_vector_set (ss, 1, beta_shape2/10);

	//Initialize method and iterate
	double par [3];
	par[0] = 0.0;
	par[1] = 0.0;
	for (int e = 0 ; e < pval.size(); e ++) {
		if (pval[e] == 1.0) pval[e] = 0.99999999;
		par[0] += log (pval[e]);
		par[1] += log (1 - pval[e]);
	}
	par[2] = pval.size();
	gsl_multimin_function minex_func;
	minex_func.n = 2;
	minex_func.f = betaLogLikelihood;
	minex_func.params = par;

	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	//Optimization iteration
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
	} while (status == GSL_CONTINUE && iter < 1000);

	//Output new beta shape values
	beta_shape1 = gsl_vector_get (s->x, 0);
	beta_shape2 = gsl_vector_get (s->x, 1);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}
