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

#ifndef _DATA_H
#define _DATA_H

#define ___NA___ (0.0/0.0)
#define ___LI___ 1000000000

#define MATHLIB_STANDALONE

#include "utils/utils.h"
#include "region.h"
#include "residualizer.h"
#include <Rmath.h>

class data {
public:
	//INCLUDE/EXCLUDE LISTS
	set < string > sample_inclusion;
	set < string > sample_exclusion;
	set < string > genotype_inclusion;
	set < string > genotype_exclusion;
	set < string > phenotype_inclusion;
	set < string > phenotype_exclusion;
	set < string > covariate_inclusion;
	set < string > covariate_exclusion;

	//REGIONS
	region regionPhenotype;
	region regionGenotype;
	float cis_window;
	float maf_threshold;								// minor allele frequency threshold
	int ma_sample_threshold;							// minor allele sample threshold
	float global_af_threshold;

	//SAMPLES
	int sample_count;									//sample number
	vector < string > sample_id;						//sample IDs

	//GENOTYPES
	int genotype_count;									//variant site number
	vector < vector < float > > genotype_orig;			//original genotype dosages
	vector < vector < float > > genotype_curr;			//current genotype dosages
	vector < string > genotype_chr;						//variant site chromosome
	vector < string > genotype_id;						//variant site IDs
	vector < int > genotype_pos;						//variant site positions
	vector < float > genotype_maf;						//variant minor allele frequency
	vector < int > genotype_ma_count;					//variant minor allele count
	vector < int > genotype_ma_samples;					//variant minor allele samples
	vector < int > genotype_ref_factor;					//variant minor allele reference factor: 1 if ALT is MA; -1 if REF is MA

	//PHENOTYPES
	int phenotype_count;								//phenotype number
	vector < vector < float > > phenotype_orig;			//original phenotype values
	vector < string > phenotype_id;						//phenotype ids
	vector < string > phenotype_grp;					//phenotype groups
	vector < string > phenotype_chr;					//phenotype chromosomes
	vector < int > phenotype_start;						//phenotype start positions
	vector < int > phenotype_end;						//phenotype end positions
	vector < vector < int > > phenotype_cluster;		//phenotype cluster (parallel jobs)
	vector < double > phenotype_threshold;				//phenotype genome-wide significance thresholds

	//COVARIATES
	int covariate_count;								//covariate number
	vector < vector < string > > covariate_val;			//covariate values
	vector < string > covariate_id;						//covariate ids
	residualizer * covariate_engine;

	//INTERACTION
	vector < float > interaction_val;					//interaction values

	//CONSTRUCTOR / DESTRUCTOR
	data();
	void clear();

	//READ EXCLUSION/INCLUSION LISTS (IO/readInclusionsExcusions.cpp)
	void readSamplesToExclude(string);
	void readSamplesToInclude(string);
	void readGenotypesToExclude(string);
	void readGenotypesToInclude(string);
	void readPhenotypesToExclude(string);
	void readPhenotypesToInclude(string);
	void readCovariatesToExclude(string);
	void readCovariatesToInclude(string);
	bool checkSample(string &, bool checkDuplicates = true);
	bool checkGenotype(string &);
	bool checkPhenotype(string &, bool include = true);
	bool checkCovariate(string &);

	//REGIONS
	bool setPhenotypeRegion(string);
	void setPhenotypeRegion(int);
	string getPhenotypeRegion(int);
	void deduceGenotypeRegion(int);

	//READ DATA
	void readGenotypesVCF(string);
	void readPhenotypes(string);
	void scanPhenotypes(string);
	void readCovariates(string);
	void readThresholds(string);
	void readGroups(string);
	void readInteractions(string);

	//GENERAL MANAGMENT
	void clusterizePhenotypes(int);
	void imputeGenotypes();
	void imputePhenotypes();
	void normalTranformPhenotypes();
	void initResidualizer();
	void rankTranformPhenotypes();
	void rankTranformGenotypes();
	void normalize(vector < float > &);
	void normalize(vector < vector < float > > &);
	void correct(vector < float > &, vector < float > &);
	int mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2);
	int learnDF(vector < double > & corr, double &);

	//COMPUTATION METHODS [ALL INLINES FOR SPEED]
	double getCorrelation(vector < float > &, vector < float > &);
	double getCorrelation(vector < float > &, vector < float > &, vector < int > &);
	double getTstat2(double, double);
	double getPvalueFromTstat2(double, double);
	double getPvalue(double, double);
	double getPvalue(double, vector < double > &);
	double getSlope(double, double, double);

	//ANALYSIS
	void runNominal(string, double);
	void runNominalBest(string);
	void runNominalOutputMatrices(string, string, double);
	void runPermutation(string, vector < int >);
	void runPermutation(string, string);
	void runPermutationPerGroup(string, vector < int >);
	void runMapping(string, bool);
	void runPermutationInteraction(string, int);

	//COMMANDS
	void writeCommands(string, int, int, char **);
};


//***************************************************************//
//******************** INLINE FUNCTIONS *************************//
//***************************************************************//

#ifdef _FAST_CORRELATION

/*
 * This implementation of inner_product is optimized for 64 bytes cache lines.
 */
inline double data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
	int i = 0;
	int repeat = (sample_count / 4);
	int left = (sample_count % 4);
	double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

	while (repeat --) {
		sum0 += vec1[i] * vec2[i];
		sum1 += vec1[i+1] * vec2[i+1];
		sum2 += vec1[i+2] * vec2[i+2];
		sum3 += vec1[i+3] * vec2[i+3];
		i += 4;
	}

	switch (left) {
	case 3:	sum0 += vec1[i+2] * vec2[i+2];
	case 2:	sum0 += vec1[i+1] * vec2[i+1];
	case 1:	sum0 += vec1[i+0] * vec2[i+0];
	case 0: ;
	}

	return sum0 + sum1 + sum2 + sum3;
}

#else

inline double data::getCorrelation(vector < float > & vec1, vector < float > & vec2) {
	double corr = 0.0;
	for (int s = 0 ; s < sample_count ; s ++) corr += vec1[s] * vec2[s];
	return corr;
}

#endif

inline double data::getCorrelation(vector < float > & vec1, vector < float > & vec2, vector < int > & indexes) {
	double corr = 0.0;
	for (int s = 0 ; s < sample_count ; s ++) corr += vec1[indexes[s]] * vec2[indexes[s]];
	return corr;
}

inline double data::getTstat2(double corr, double df) {
	return df * corr * corr / (1 - corr * corr);
}

inline double data::getPvalueFromTstat2(double tstat2, double df) {
	return pf(tstat2, 1, df, 0, 0);
}

inline double data::getPvalue(double corr, double df) {
	return pf(df * corr * corr / (1 - corr * corr), 1, df,0,0);
}

inline double data::getPvalue(double ncorr, vector < double > & pcorr) {
	unsigned int n_hits = 0;
	for (int p = 0 ; p < pcorr.size() ; p++) if (abs(pcorr[p]) >= abs(ncorr)) n_hits++;
	return ((n_hits + 1) * 1.0 / (pcorr.size() + 1.0));
}

inline double data::getSlope(double nominal_correlation, double psd, double gsd) {
	if (gsd < 1e-16 || psd < 1e-16) return 0;
	else return nominal_correlation * psd / gsd;
}

#endif
