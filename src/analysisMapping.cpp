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

#define MAX_ANALYSIS_DEPTH	10

void data::runMapping(string fout, bool full) {
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		LOG.println("\nProcessing gene [" + phenotype_id[p] + "]");

		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		if (targetGenotypes.size() > 0) {

			LOG.println("  * Nominal significance threshold = " + sutils::double2str(phenotype_threshold[p]));

			// 1. Forward pass: Learn number of independent signals and Map the best candidates
			vector < double > bestCorr = vector < double > (MAX_ANALYSIS_DEPTH, 0.0);
			vector < double > uncorrected_pvalues = vector < double > (targetGenotypes.size(), 2.0);
			vector < int > bestIndex;
			bool done = false;

			for (int i = 0 ; i < MAX_ANALYSIS_DEPTH && !done; i ++) {
				int n_significant = 0;
				vector < float > phenotype_curr = phenotype_orig[p];
				copy(genotype_orig.begin() + targetGenotypes[0], genotype_orig.begin() + targetGenotypes.back() + 1, genotype_curr.begin() + targetGenotypes[0]);
				vector < int > bestIndex_tmp = bestIndex;
				bestIndex_tmp.push_back(-1);

				//Covariates + Best signals
				covariate_engine->clearSoft();
				for (int h = 0 ; h < bestIndex.size() ; h ++) covariate_engine->pushSoft(genotype_orig[bestIndex[h]]);
				covariate_engine->residualize(phenotype_curr);
				normalize(phenotype_curr);

				for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
					if (i == 0 || full || (!full && uncorrected_pvalues[g] <= phenotype_threshold[p])) {
						covariate_engine->residualize(genotype_curr[targetGenotypes[g]]);
						normalize(genotype_curr[targetGenotypes[g]]);
						double corr = getCorrelation(genotype_curr[targetGenotypes[g]], phenotype_curr);
						double pvalue = getPvalue(corr, sample_count - 2 - covariate_engine->nCovariates());
						if (abs(corr) > abs(bestCorr[i])) {
							bestCorr[i] = corr;
							bestIndex_tmp[i] = targetGenotypes[g];
						}
						if (i == 0) uncorrected_pvalues[g] = pvalue;
						if (pvalue <= phenotype_threshold[p]) n_significant ++;
					}
				}
				if (n_significant == 0) done = true;
				else bestIndex = bestIndex_tmp;
			}
			LOG.println("  * Number of independent signals found = " + sutils::int2str(bestIndex.size()));

			if (bestIndex.size() == 1) {
				int n_signals = 0;
				for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
					if (uncorrected_pvalues[g] <= phenotype_threshold[p]) {
						fdo << phenotype_id[p] << " 0 " << genotype_id[targetGenotypes[g]] << " " << targetDistances[g] << " " << uncorrected_pvalues[g] << " " << uncorrected_pvalues[g] << " " << (bestIndex[0] == targetGenotypes[g]) << endl;
						n_signals ++;
					}
				}
				LOG.println("  * Number of candidate QTLs reported for rank 1 = " + sutils::int2str(n_signals));
			} else if (bestIndex.size() > 1) {
				//2. Backward pass: Determine candidate variants and classify them
				vector < vector < double > > corrected_pvalues = vector < vector < double > > (bestIndex.size(), vector < double > (targetGenotypes.size(), 2.0));
				for (int i = bestIndex.size() - 1 ; i >= 0 ; i --) {
					vector < float > phenotype_curr = phenotype_orig[p];
					copy(genotype_orig.begin() + targetGenotypes[0], genotype_orig.begin() + targetGenotypes.back() + 1, genotype_curr.begin() + targetGenotypes[0]);

					vector < int > bestIndexOthers = bestIndex;
					bestIndexOthers.erase(bestIndexOthers.begin() + i);

					covariate_engine->clearSoft();
					for (int h = 0 ; h < bestIndexOthers.size() ; h ++) covariate_engine->pushSoft(genotype_orig[bestIndexOthers[h]]);
					covariate_engine->residualize(phenotype_curr);
					normalize(phenotype_curr);

					for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
						if (full || (!full && uncorrected_pvalues[g] <= phenotype_threshold[p])) {
							covariate_engine->residualize(genotype_curr[targetGenotypes[g]]);
							normalize(genotype_curr[targetGenotypes[g]]);
							double corr = getCorrelation(genotype_curr[targetGenotypes[g]], phenotype_curr);
							corrected_pvalues[i][g] = getPvalue(corr, sample_count - 2 - covariate_engine->nCovariates());
						}
					}
				}

				//
				for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
					double min_pvalue = 1.1;
					for (int i = 0; i < bestIndex.size() ; i ++) if (corrected_pvalues[i][g] < min_pvalue) min_pvalue = corrected_pvalues[i][g];
					for (int i = 0; i < bestIndex.size() ; i ++) if (corrected_pvalues[i][g] != min_pvalue) {
						corrected_pvalues[i][g] = 2.0;
					}
				}

				//
				vector < int > n_signals = vector < int > (bestIndex.size(), 0);
				for (int i = 0; i < bestIndex.size() ; i ++) {
					for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
						if ((bestIndex[i] == targetGenotypes[g]) || (corrected_pvalues[i][g] <= phenotype_threshold[p])) {
							fdo << phenotype_id[p] << " " << i << " " << genotype_id[targetGenotypes[g]] << " " << targetDistances[g] << " " << uncorrected_pvalues[g] << " " << corrected_pvalues[i][g] << " " << (bestIndex[i] == targetGenotypes[g]) << endl;
							n_signals [i] ++;
						}
					}
				}

				//
				for (int i = 0; i < bestIndex.size() ; i ++)
					LOG.println("  * Number of candidate QTLs reported for rank "+ sutils::int2str(i) + " = " + sutils::int2str(n_signals[i]));
			}
		}

		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
