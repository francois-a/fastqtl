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


void data::runNominalInteractionBest(string fout) {

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
        int min_idx = 0;
        double min_pval = 1.0;
        double opt_corr = 0.0;
        double opt_tstat2 = 0.0;
        double opt_phenotype_sd = 0.0;
        double opt_interaction_sd = 0.0;
        for (int g = 0 ; g < targetGenotypes.size() ; g ++) {

            //RESIDUALIZER
            covariate_engine->clearSoft();
            covariate_engine->pushSoft(genotype_orig[targetGenotypes[g]]);

            //INTERACTION TERM
            vector < float > interaction_curr = genotype_orig[targetGenotypes[g]];
            for (int i = 0 ; i < sample_count ; i++) interaction_curr[i] *= interaction_val[i];
            covariate_engine->residualize(interaction_curr);
            double interaction_sd = RunningStat(interaction_curr).StandardDeviation();
            normalize(interaction_curr);

            //NOMINAL STEP: scan cis-window for best variant & compute statistics
            vector < float > phenotype_curr = phenotype_orig[p];
            covariate_engine->residualize(phenotype_curr);
            double phenotype_sd = RunningStat(phenotype_curr).StandardDeviation();
            normalize(phenotype_curr);
            double corr = getCorrelation(interaction_curr, phenotype_curr);
            double df = sample_count - 2 - covariate_engine->nCovariates();
            double tstat2 = getTstat2(corr, df);
            double pval = getPvalueFromTstat2(tstat2, df);
            if (pval < min_pval) {
                min_idx = g;
                min_pval = pval;
                opt_corr = corr;
                opt_tstat2 = tstat2;
                opt_phenotype_sd = phenotype_sd;
                opt_interaction_sd = interaction_sd;
            }
        }
        if (targetGenotypes.size()>0) {
            int gi = targetGenotypes[min_idx];
            double slope = getSlope(opt_corr, opt_phenotype_sd, opt_interaction_sd);
            double slope_se = abs(slope) / sqrt(opt_tstat2);
            fdo << phenotype_id[p];
            fdo << "\t" << genotype_id[gi];
            fdo << "\t" << targetDistances[min_idx];
            fdo << "\t" << genotype_ma_samples[gi];
            fdo << "\t" << genotype_ma_count[gi];
            fdo << "\t" << genotype_maf[gi];
            fdo << "\t" << min_pval;
            fdo << "\t" << slope;
            fdo << "\t" << slope_se;
            fdo << endl;
        }
        LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
    }
    fdo.close();
}
