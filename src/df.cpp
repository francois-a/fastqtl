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

struct data_to_function {
	data * D;
	int n;
	double * C;

	data_to_function (data * _D, int _n, double * _C) {
		D = _D;
		n = _n;
		C = _C;
	}

};

double degreeOfFreedom(const gsl_vector *v, void *params) {
	data_to_function * d = (data_to_function *) params;
	vector < double > pval = vector < double >(d->n, 0.0);
	double mean = 0.0;
	for (int c = 0 ; c < d->n ; c++) {
		pval[c] = d->D->getPvalue(d->C[c], gsl_vector_get(v, 0));
		mean += pval[c];
	}
	mean /= pval.size();
	double variance = 0.0;
	for (int c = 0 ; c < pval.size() ; c++) variance += (pval[c] - mean) * (pval[c] - mean);
	variance /= (pval.size() - 1);

	double shape2 = abs((mean * (mean * (1 - mean ) / variance - 1)) - 1.0);
	//cout << "O = " << mean << " " << shape2 << endl;
	return shape2;
}

int data::learnDF(vector < double > & corr, double & df) {

	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, df);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (1);
	gsl_vector_set (ss, 0, df * 0.1);

	//Initialize method and iterate
	data_to_function * par  = new data_to_function (this, corr.size(), &corr[0]);

	gsl_multimin_function minex_func;
	minex_func.n = 1;
	minex_func.f = degreeOfFreedom;
	minex_func.params = (void*)par;

	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 1);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	//Optimization iteration
	//cout << "\n ========================" << endl;
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
		//printf ("%d %10.3e f() = %7.10f size = %.10f\n", iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < 20);

	//Output new beta shape values
	df = gsl_vector_get (s->x, 0);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}
