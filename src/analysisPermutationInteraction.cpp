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


void data::runPermutationInteraction(string fout, int nPermutations) {

	//0. Loop over phenotypes
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		LOG.println("\nProcessing gene [" + phenotype_id[p] + "]");

		//1.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		//1.2. Loop over genotypes
		double nominal_correlation = 0.0;
		int nominal_variant = -1;
		int nominal_distance = -1;
		vector < double > permuted_correlations = vector < double > (nPermutations, 0.0);
		vector < vector < float > > permuted_phenotypes = vector < vector < float > > (nPermutations, vector < float > (sample_count, 0));

		for (int perm = 0 ; perm < nPermutations ; perm ++) {
			permuted_phenotypes[perm] = phenotype_orig[p];
			random_shuffle(permuted_phenotypes[perm].begin(), permuted_phenotypes[perm].end());
		}

		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {

			//RESIDUALIZER
			covariate_engine->clearSoft();
			covariate_engine->pushSoft(genotype_orig[targetGenotypes[g]]);

			//INTERACTION TERM
			vector < float > interaction_curr = genotype_orig[targetGenotypes[g]];
			for (int i = 0 ; i < sample_count ; i++) interaction_curr[i] *= interaction_val[i];
			covariate_engine->residualize(interaction_curr);
			normalize(interaction_curr);

			//NOMINAL STEP
			vector < float > phenotype_curr = phenotype_orig[p];
			covariate_engine->residualize(phenotype_curr);
			normalize(phenotype_curr);
			double corr = getCorrelation(interaction_curr, phenotype_curr);
			if (abs(corr) > abs(nominal_correlation)) {
				nominal_correlation = corr;
				nominal_variant = targetGenotypes[g];
				nominal_distance = targetDistances[g];
			}

			//PERMUTATION PASS
			for (int perm = 0 ; perm < nPermutations ; perm ++) {
				random_shuffle(phenotype_curr.begin(), phenotype_curr.end());
				phenotype_curr = permuted_phenotypes[perm];
				covariate_engine->residualize(phenotype_curr);
				normalize(phenotype_curr);
				corr = getCorrelation(interaction_curr, phenotype_curr);
				if (abs(corr) > abs(permuted_correlations[perm])) permuted_correlations[perm] = corr;
			}
		}

		//1.5. Calculate effective DFs & Beta distribution parameters
		vector < double > permuted_pvalues = permuted_correlations;
		double true_df = sample_count - 2 - covariate_engine->nCovariates();
		double mean = 0.0, variance = 0.0, beta_shape1 = 1.0, beta_shape2 = 1.0;
		fdo << phenotype_id[p] << " " << targetGenotypes.size();
		if (targetGenotypes.size() > 0) {
			//Estimate number of degrees of freedom
			if (putils::variance(permuted_correlations, putils::mean(permuted_correlations)) != 0.0) {
				learnDF(permuted_correlations, true_df);
				//LOG.println("  * Effective degree of freedom = " + sutils::double2str(true_df, 4));
			}
			//Compute mean and variance of p-values
			for (int c = 0 ; c < permuted_correlations.size() ; c ++) permuted_pvalues[c] = getPvalue(permuted_correlations[c], true_df);
			mean = putils::mean(permuted_pvalues);
			for (int pv = 0 ; pv < permuted_pvalues.size() ; pv++) variance += (permuted_pvalues[pv] - mean) * (permuted_pvalues[pv] - mean);
			variance /= (permuted_pvalues.size() - 1);
			//Estimate shape1 & shape2
			if (targetGenotypes.size() > 1 && mean != 1.0) {
				beta_shape1 = mean * (mean * (1 - mean ) / variance - 1);
				beta_shape2 = beta_shape1 * (1 / mean - 1);
				if (targetGenotypes.size() > 10) mleBeta(permuted_pvalues, beta_shape1, beta_shape2);	//ML estimate if more than 10 variant in cis
			}
			LOG.println("  * Beta distribution parameters = " + sutils::double2str(beta_shape1, 4) + " " + sutils::double2str(beta_shape2, 4));
			fdo << " " << beta_shape1 << " " << beta_shape2 << " " << true_df;
		} else fdo << " NA NA NA";

		//1.6. Writing results
		if (targetGenotypes.size() > 0 && nominal_variant >= 0) {
			double pval_fdo = getPvalue(nominal_correlation, true_df);
			double pval_nom = getPvalue(nominal_correlation, sample_count - 2 - covariate_engine->nCovariates());
			fdo << " " << genotype_id[nominal_variant];
			fdo << " " << nominal_distance;
			fdo << " " << pval_nom;
			fdo << " " << getPvalue(nominal_correlation, permuted_correlations);
			fdo << " " << pbeta(pval_fdo, beta_shape1, beta_shape2, 1, 0);
		} else fdo << " NA NA NA NA NA";
		fdo << endl;

		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
