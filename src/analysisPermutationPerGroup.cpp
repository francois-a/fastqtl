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


void data::runPermutationPerGroup(string fout, vector < int > nPermutations) {

	//0. Prepare genotypes
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes for covariates");
		covariate_engine->residualize(genotype_orig);
	}
	normalize(genotype_orig);

	//1. Prepare phenotype groups
	LOG.println("\nBuilding groups");
	vector < vector < int > > PG;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		if (p == 0 || phenotype_grp [p-1] != phenotype_grp [p]) PG.push_back(vector < int > (1, p));
		else PG.back().push_back(p);
	}
	LOG	.println("  * Number of groups = " + sutils::int2str(PG.size()));

	//2. Loop over groups
	ofile fdo (fout);
	for (int g = 0 ; g < PG.size() ; g ++) {

		LOG.println("\nProcessing group [" + phenotype_grp[PG[g][0]] + "]");
		LOG.println("  * Number of MP in group = " + sutils::int2str(PG[g].size()) + " / " + sutils::int2str(PG[g]));

		//1.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int l = 0 ; l < genotype_count ; l ++) {
			int cisdistance = genotype_pos[l] - phenotype_start[PG[g][0]];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(l);
				targetDistances.push_back(cisdistance);
			}
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		//1.2. Copy original data
		vector < vector < float > > phenotype_curr = vector < vector < float > > (phenotype_orig.begin() + PG[g][0], phenotype_orig.begin() + PG[g].back() + 1);
		if (covariate_count > 0) covariate_engine->residualize(phenotype_curr);
		normalize(phenotype_curr);

		//1.3. Nominal pass: scan cis-window & compute statistics
		double bestCorr = 0.0;
		int bestDistance = ___LI___, bestIndex = -1;
        string bestExon = "NA";
		for (int p = 0 ; p < phenotype_curr.size() ; p ++) {
			for (int l = 0 ; l < targetGenotypes.size() ; l ++) {
				double corr = getCorrelation(genotype_orig[targetGenotypes[l]], phenotype_curr[p]);
				if (corr > 1) {
					cerr << endl;
					cerr << "************************************************************************************************" << endl;
					cerr << g << "/" << PG.size() << " " << p << "/" << phenotype_curr.size() << " " << l << "/" << targetGenotypes.size() << endl;
					cerr << "-----Genotypes--------" << endl;
					cerr << sutils::double2str(genotype_orig[targetGenotypes[l]]) << endl;
					cerr << "-----Phenotypes-------" << endl;
					cerr << sutils::double2str(phenotype_curr[p]) << endl;
				}
				if (abs(corr) > abs(bestCorr) || (abs(corr) == abs(bestCorr) && abs(targetDistances[l]) < bestDistance)) {
					bestCorr = corr;
					bestDistance = targetDistances[l];
					bestIndex = targetGenotypes[l];
                    bestExon = phenotype_id[PG[g][p]];
				}
			}
		}
		if (targetGenotypes.size() > 0) {

			LOG.println("  * Best correlation = " + sutils::double2str(bestCorr, 4));
		}

		//1.4. Permutation pass:
		bool done = false;
		int countPermutations = 0, nBetterCorrelation = 0;
		vector < double > permCorr;
		do {
			double bestCperm = 0.0;
			copy(phenotype_orig.begin() + PG[g][0], phenotype_orig.begin() + PG[g].back() + 1, phenotype_curr.begin());
			putils::random_shuffle(phenotype_curr);
			if (covariate_count > 0) covariate_engine->residualize(phenotype_curr);
			normalize(phenotype_curr);
			for (int p = 0 ; p < phenotype_curr.size() ; p ++) {
				for (int l = 0 ; l < targetGenotypes.size() ; l ++) {
					double corr = getCorrelation(genotype_orig[targetGenotypes[l]], phenotype_curr[p]);
					if (abs(corr) > abs(bestCperm)) bestCperm = corr;
				}
			}
			if (abs(bestCperm) >= abs(bestCorr)) nBetterCorrelation++;
			permCorr.push_back(bestCperm);
			countPermutations++;

			if (nPermutations.size() == 1 && countPermutations >= nPermutations[0]) done = true;
			if (nPermutations.size() == 2 && (nBetterCorrelation >= nPermutations[0] || countPermutations >= nPermutations[1])) done = true;
			if (nPermutations.size() == 3 && (countPermutations >= nPermutations[0]) && (nBetterCorrelation >= nPermutations[1] || countPermutations >= nPermutations[2])) done = true;
		} while (!done);
		if (targetGenotypes.size() > 0) LOG.println("  * Number of permutations = " + sutils::int2str(nBetterCorrelation) + " / " + sutils::int2str(countPermutations));

		//1.5. Calculate effective DFs & Beta distribution parameters
		vector < double > permPvalues;
		double true_df = sample_count - 2 - ((covariate_count>0)?covariate_engine->nCovariates():0);
		double mean = 0.0, variance = 0.0, beta_shape1 = 1.0, beta_shape2 = 1.0;
		fdo << bestExon << " " << targetGenotypes.size();
		if (targetGenotypes.size() > 0) {
			//Estimate number of degrees of freedom
			if (putils::variance(permCorr, putils::mean(permCorr)) != 0.0) {
				learnDF(permCorr, true_df);
				LOG.println("  * Effective degree of freedom = " + sutils::double2str(true_df, 4));
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
			fdo << " " << (nBetterCorrelation + 1) * 1.0 / (countPermutations + 1.0);
			fdo << " " << pbeta(pval_fdo, beta_shape1, beta_shape2, 1, 0);
            fdo << " " << phenotype_grp[PG[g][0]];
            fdo << " " << PG[g].size();
		} else fdo << " NA NA NA NA NA NA NA";
		fdo << endl;

		LOG.println("  * Progress = " + sutils::double2str((g+1) * 100.0 / PG.size(), 1) + "%");
	}
	fdo.close();
}
