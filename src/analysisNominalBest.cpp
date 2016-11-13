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


void data::runNominalBest(string fout) {

	//0. Prepare genotypes
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes & phenotypes for covariates");
		covariate_engine->residualize(genotype_orig);
		covariate_engine->residualize(phenotype_orig);
	}
	for (int g = 0 ; g < genotype_count ; g ++) genotype_sd[g] = RunningStat(genotype_orig[g]).StandardDeviation();
	for (int p = 0 ; p < phenotype_count ; p ++) phenotype_sd[p] = RunningStat(phenotype_orig[p]).StandardDeviation();
	normalize(genotype_orig);
	normalize(phenotype_orig);

	double df = sample_count - 2 - covariate_count;

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

		//1.2. Nominal pass: scan cis-window for best variant & compute statistics
		int min_idx = 0;
		double min_pval = 1.0;
		double opt_corr = 0.0;
		double opt_tstat2 = 0.0;
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_orig[p]);
			double tstat2 = getTstat2(corr, df);
			double pval = getPvalueFromTstat2(tstat2, df);
			if (pval<min_pval) {
				min_idx = g;
				min_pval = pval;
				opt_corr = corr;
				opt_tstat2 = tstat2;
			}
		}
		if (targetGenotypes.size()>0) {
			double slope = getSlope(opt_corr, phenotype_sd[p], genotype_sd[targetGenotypes[min_idx]]);
			double slope_se = abs(slope) / sqrt(opt_tstat2);
			fdo << phenotype_id[p];
			fdo << "\t" << genotype_id[targetGenotypes[min_idx]];
			fdo << "\t" << targetDistances[min_idx];
			fdo << "\t" << min_pval;
			fdo << "\t" << slope;
			fdo << "\t" << slope_se;
			fdo << endl;
		}
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
