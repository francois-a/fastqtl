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


void data::runPermutation(string fout, string fperm) {
	string buffer;
	vector < string > tokens;

	//Read perm sequence
	ifile fdp (fperm);
	vector < vector < unsigned int > > P;
	while(getline(fdp, buffer)) {
		P.push_back(vector < unsigned int > ());
		sutils::tokenize(buffer, tokens);
		for (int t = 0 ; t < tokens.size() ; t++ ) P.back().push_back(atoi(tokens[t].c_str()) - 1);
	}

	//0. Prepare genotypes
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes for covariates");
		covariate_engine->residualize(genotype_orig);
	}
	normalize(genotype_orig);

	//1. Loop over phenotypes
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

		//1.2. Copy original data
		vector < float > phenotype_curr = phenotype_orig[p];
		if (covariate_count > 0) covariate_engine->residualize(phenotype_curr);
		normalize(phenotype_curr);

		//1.3. Nominal pass: scan cis-window & compute statistics
		double bestCorr = 0.0;
		int bestDistance = ___LI___, bestIndex = -1;
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_curr);
			if (abs(corr) > abs(bestCorr) || (abs(corr) == abs(bestCorr) && abs(targetDistances[g]) < bestDistance)) {
				bestCorr = corr;
				bestDistance = targetDistances[g];
				bestIndex = targetGenotypes[g];
			}
		}
		if (targetGenotypes.size() > 0) LOG.println("  * Best correlation = " + sutils::double2str(bestCorr, 4));

		//1.4. Permutation pass:
		int nBetterCorrelation = 0;
		vector < double > permCorr;
		for (int perm = 0 ; perm < P.size() ; perm ++) {
			double bestCperm = 0.0;
			autils::reorder(phenotype_orig[p], P[perm]);
			phenotype_curr = phenotype_orig[p];
			if (covariate_count > 0) covariate_engine->residualize(phenotype_curr);
			normalize(phenotype_curr);
			for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
				double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_curr);
				if (abs(corr) > abs(bestCperm)) bestCperm = corr;
			}
			if (abs(bestCperm) >= abs(bestCorr)) nBetterCorrelation++;
			permCorr.push_back(bestCperm);
		}
		if (targetGenotypes.size() > 0) LOG.println("  * Number of permutations = " + sutils::int2str(nBetterCorrelation) + " / " + sutils::int2str(P.size()));

		//1.5. Calculate effective DFs & Beta distribution parameters
		vector < double > permPvalues;
		double true_df = sample_count - 2 - ((covariate_count>0)?covariate_engine->nCovariates():0);
		double mean = 0.0, variance = 0.0, beta_shape1 = 1.0, beta_shape2 = 1.0;
		fdo << phenotype_id[p] << " " << targetGenotypes.size();
		if (targetGenotypes.size() > 0) {
			//Estimate number of degrees of freedom
			if (putils::variance(permCorr, putils::mean(permCorr)) != 0.0) {
				learnDF(permCorr, true_df);
				//LOG.println("  * Effective degree of freedom = " + sutils::double2str(true_df, 4));
			}
			//Compute mean and variance of p-values
			for (int c = 0 ; c < permCorr.size() ; c ++) permPvalues.push_back(getPvalue(permCorr[c], true_df));
			for (int pv = 0 ; pv < permPvalues.size() ; pv++) mean += permPvalues[pv];
			mean /= permPvalues.size();
			for (int pv = 0 ; pv < permPvalues.size() ; pv++) variance += (permPvalues[pv] - mean) * (permPvalues[pv] - mean);
			variance /= (permPvalues.size() - 1);
			//Estimate shape1 & shape2
			if (targetGenotypes.size() > 1 && mean != 1.0) {
				beta_shape1 = mean * (mean * (1 - mean ) / variance - 1);
				beta_shape2 = beta_shape1 * (1 / mean - 1);
				if (targetGenotypes.size() > 10) mleBeta(permPvalues, beta_shape1, beta_shape2);	//ML estimate if more than 10 variant in cis
			}
			LOG.println("  * Beta distribution parameters = " + sutils::double2str(beta_shape1, 4) + " " + sutils::double2str(beta_shape2, 4));
			fdo << " " << beta_shape1 << " " << beta_shape2 << " " << true_df;
		} else fdo << " NA NA NA";

		//1.6. Writing results
		if (targetGenotypes.size() > 0 && bestIndex >= 0) {
			double pval_fdo = getPvalue(bestCorr, true_df);
			double pval_nom = getPvalue(bestCorr, sample_count - 2 - ((covariate_count>0)?covariate_engine->nCovariates():0));
			fdo << " " << genotype_id[bestIndex];
			fdo << " " << bestDistance;
			fdo << " " << pval_nom;
			fdo << " " << (nBetterCorrelation + 1) * 1.0 / (P.size() + 1.0);
			fdo << " " << pbeta(pval_fdo, beta_shape1, beta_shape2, 1, 0);
		} else fdo << " NA NA NA NA NA";
		fdo << endl;

		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
