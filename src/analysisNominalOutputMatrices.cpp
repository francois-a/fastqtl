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
#include "cnpy.h"

void data::runNominalOutputMatrices(string dout, string fout, double threshold) {

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

    // -------------------------------------------------------------------------------------------------------------- 
    // Save the data to file 
    // pass 1: count how many variants there are per gene, save the variants to file as well
    LOG.println("\nSaving genes and variant information");
    string g_out = dout + "/Genes.txt";
	if (!futils::createFile(g_out)) LOG.error(g_out + " is impossible to create, check permissions");
    ofile g_file (g_out);
    int max_target_genotypes = 0;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		// enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
        int tg_size = (int) targetGenotypes.size();
		LOG.println("  * Gene: " + phenotype_id[p] + "\t Number of variants in cis = " + sutils::int2str(tg_size));
        max_target_genotypes = max(max_target_genotypes, tg_size);
        // save the gene name and cis var counts 
        g_file << phenotype_id[p] << "\t" << tg_size <<  "\n";
        // save the variant names to file
        string v_out = dout+"/Vars_" + sutils::int2str(p) + ".txt";
        ofile v_file (v_out);
        for (int g = 0 ; g < targetGenotypes.size() ; g ++)    
            v_file << genotype_id[targetGenotypes[g]] << "\n";
        v_file.close();
    }
    g_file.close();
    LOG.println("Maximum number of variants in cis = " + sutils::int2str(max_target_genotypes));
    LOG.println("Written genes to: " + g_out); 

    // pass 2: save the data to npy files (used info from pass 1 to allocate memory)
    LOG.println("\nSaving data (X, y)");
    float* x_data = new float[sample_count*max_target_genotypes];
    float* y_data = new float[sample_count];
	for (int p = 0 ; p < phenotype_count ; p ++) {
		//0.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
			int cisdistance = genotype_pos[g] - phenotype_start[p];
			if (abs(cisdistance) <= cis_window) {
				targetGenotypes.push_back(g);
				targetDistances.push_back(cisdistance);
			}
		}
        int tg_size = (int) targetGenotypes.size();
        //0.2 Save the responses for this gene
        string y_out = dout+"/y_" + sutils::int2str(p) + ".npy";
        string x_out = dout+"/X_" + sutils::int2str(p) + ".npy";
        const unsigned int x_shape[] = {sample_count,tg_size};
        const unsigned int y_shape[] = {sample_count};
        for(unsigned i=0;i<x_shape[0];i++){
            for(unsigned j=0;j<x_shape[1];j++){
                x_data[i*x_shape[1]+j] = genotype_orig[targetGenotypes[j]][i];
            }
            y_data[i] = phenotype_orig[p][i];
        }
        cnpy::npy_save(x_out,x_data,x_shape,2,"w");
		LOG.println("  * Saved eQTL X matirx to [" + x_out + "]");
        cnpy::npy_save(y_out,y_data,y_shape,1,"w");
		LOG.println("  * Saved eQTL y vector to [" + y_out + "]");
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
    }
    delete[] x_data;
    delete[] y_data;
    // -------------------------------------------------------------------------------------------------------------- 


	normalize(genotype_orig);
	normalize(phenotype_orig);

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

		//1.2. Nominal pass: scan cis-window & compute statistics
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_orig[p]);
			double df = sample_count - 2 - covariate_count;
			double tstat2 = getTstat2(corr, df);
			double pval = getPvalueFromTstat2(tstat2, df);
			double slope = getSlope(corr, phenotype_sd[p], genotype_sd[targetGenotypes[g]]);
			double slope_se = abs(slope) / sqrt(tstat2);
			if (pval <= threshold ) {
				fdo << phenotype_id[p];
				fdo << "\t" << genotype_id[targetGenotypes[g]];
				fdo << "\t" << targetDistances[g];
				fdo << "\t" << pval;
				fdo << "\t" << slope;
				fdo << "\t" << slope_se;
				fdo << endl;
			}
		}
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}

