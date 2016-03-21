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


#define _GLOBAL

#include "data.h"
#include "utils/ranker.h"


data::data() {
	sample_count = 0;
	genotype_count = 0;
	phenotype_count = 0;
	covariate_count = 0;
	covariate_engine = NULL;
}

void data::clear() {
	sample_count = 0;
	sample_id.clear();
	genotype_count = 0;
	genotype_orig.clear();
	genotype_curr.clear();
	genotype_chr.clear();
	genotype_id.clear();
	genotype_pos.clear();
	phenotype_count = 0;
	phenotype_orig.clear();
	phenotype_id.clear();
	phenotype_chr.clear();
	phenotype_start.clear();
	phenotype_end.clear();
	covariate_count = 0;
	covariate_val.clear();
	covariate_id.clear();
}

bool data::checkSample(string & id, bool checkDuplicate) {
	bool included = ((sample_inclusion.size() == 0)?true:sample_inclusion.count(id));
	bool excluded = ((sample_exclusion.size() == 0)?false:sample_exclusion.count(id));
	if (!included || excluded) return false;
	if (checkDuplicate) for (int s = 0 ; s < sample_id.size() ; s++) if (id == sample_id[s]) LOG.error("Duplicate sample id [" + id +"]");
	return true;
}

bool data::checkGenotype(string & id) {
	bool included = ((genotype_inclusion.size() == 0)?true:genotype_inclusion.count(id));
	bool excluded = ((genotype_exclusion.size() == 0)?false:genotype_exclusion.count(id));
	if (!included || excluded) return false;
	//for (int g = 0 ; g < genotype_id.size() ; g++) if (id == genotype_id[g]) LOG.error("Duplicate variant site id [" + id + "]");
	return true;
}

bool data::checkPhenotype(string & id, bool include) {
	if (include) {
		bool included = ((phenotype_inclusion.size() == 0)?true:phenotype_inclusion.count(id));
		bool excluded = ((phenotype_exclusion.size() == 0)?false:phenotype_exclusion.count(id));
		if (!included || excluded) return false;
	}
	for (int p = 0 ; p < phenotype_id.size() ; p++) if (id == phenotype_id[p]) LOG.error("Duplicate phenotype id [" + id + "]");
	return true;
}

bool data::checkCovariate(string & id) {
	bool included = ((covariate_inclusion.size() == 0)?true:covariate_inclusion.count(id));
	bool excluded = ((covariate_exclusion.size() == 0)?false:covariate_exclusion.count(id));
	if (!included || excluded) return false;
	for (int c = 0 ; c < covariate_id.size() ; c++) if (id == covariate_id[c]) LOG.error("Duplicate covariate id [" + id + "]");
	return true;
}

void data::clusterizePhenotypes(int K) {
	//check K
	if (K < 1) LOG.error("Number of chunks needs to be > 0");
	if (K > phenotype_count) LOG.error("Number of chunks (" + sutils::int2str(K) + ") is greater than the number of phenotypes (" + sutils::int2str(phenotype_count) + ")");

	//initialise (1 cluster = 1 chromosome)
	phenotype_cluster = vector < vector < int > > (1, vector < int > (1, 0));
	for (int p = 1 ; p < phenotype_count ; p ++) {
		if (phenotype_chr[p] != phenotype_chr[p-1]) phenotype_cluster.push_back(vector < int > (1, p));
		else phenotype_cluster.back().push_back(p);
	}

	//iterate (split cluster in the middle until K clusters are built)
	if (phenotype_cluster.size() > K) LOG.error("Number of chunks needs to be > #chr");
	if (phenotype_cluster.size() < K) {
		int max_idx, max_val, max_mid;
		do {
			max_idx = -1;
			max_val = +1;
			max_mid = -1;
			for (int k = 0 ; k < phenotype_cluster.size() ; k ++) {
				if (phenotype_cluster[k].size() > max_val) {
					max_val = phenotype_cluster[k].size();
					max_idx = k;
				}
			}
			if (max_idx >= 0) {
				max_mid = max_val / 2;
				while (max_mid > 1 && phenotype_end[phenotype_cluster[max_idx][max_mid-1]] >= phenotype_start[phenotype_cluster[max_idx][max_mid]]) max_mid --;
				phenotype_cluster.push_back(vector < int > (phenotype_cluster[max_idx].begin() + max_mid, phenotype_cluster[max_idx].end()));
				phenotype_cluster[max_idx].erase(phenotype_cluster[max_idx].begin() + max_mid, phenotype_cluster[max_idx].end());
			}
		} while (phenotype_cluster.size() < K && max_idx >= 0);
	}
}

bool data::setPhenotypeRegion(string reg) {
	return regionPhenotype.set(reg);
}

void data::setPhenotypeRegion(int k) {
	regionPhenotype.chr = phenotype_chr[phenotype_cluster[k][0]];
	regionPhenotype.start = phenotype_start[phenotype_cluster[k][0]];
	regionPhenotype.end = phenotype_end[phenotype_cluster[k].back()];
}

string data::getPhenotypeRegion(int k) {
	return string (phenotype_chr[phenotype_cluster[k][0]] + ":" + sutils::int2str(phenotype_start[phenotype_cluster[k][0]]) + "-" + sutils::int2str(phenotype_end[phenotype_cluster[k].back()]));
}

void data::deduceGenotypeRegion(int W) {
	regionGenotype.chr = regionPhenotype.chr;
	regionGenotype.start = regionPhenotype.start - W;
	if (regionGenotype.start < 0) regionGenotype.start = 0;
	regionGenotype.end = regionPhenotype.end + W;
}

void data::imputeGenotypes() {
	LOG.println("\nImputing missing genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		double mean = 0.0;
		int c_mean = 0;
		for (int s = 0; s < sample_count ; s ++) {
			if (genotype_orig[g][s] >= 0) {
				mean += genotype_orig[g][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (genotype_orig[g][s] < 0) genotype_orig[g][s] = mean;
	}
}

void data::imputePhenotypes() {
	LOG.println("\nImputing missing phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		double mean = 0.0;
		int c_mean = 0;
		for (int s = 0; s < sample_count; s ++) {
			if (!isnan(phenotype_orig[p][s])) {
				mean += phenotype_orig [p][s];
				c_mean ++;
			}
		}
		mean /= c_mean;
		for (int s = 0; s < sample_count ; s ++) if (isnan(phenotype_orig[p][s])) phenotype_orig[p][s] = mean;
	}
}

void data::normalTranformPhenotypes() {
	string method = "average";
	LOG.println("\nQuantile normalise phenotypes to Normal distribution");
	for (int p = 0; p < phenotype_count ; p ++) {
		vector < float > R;
        myranker::rank(phenotype_orig[p], R, method);
		double max = 0;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] = R[s] - 0.5;
			if (R[s] > max) max = R[s];
		}
		max = max + 0.5;
		for (int s = 0 ; s < sample_count ; s ++) {
			R[s] /= max;
			phenotype_orig[p][s] = qnorm(R[s], 0.0, 1.0, 1, 0);
		}
	}
}

void data::initResidualizer() {
	LOG.println("\nInitialize covariate ");
	covariate_engine = new residualizer (sample_count);
	for (int c = 0 ; c < covariate_count ; c ++) covariate_engine->pushHard(covariate_val[c]);
	if (interaction_val.size() > 0) covariate_engine->pushHard(interaction_val);
}

void data::rankTranformPhenotypes() {
	string method = "average";
	LOG.println("\nRanking phenotypes");
	for (int p = 0; p < phenotype_count ; p ++) {
		vector < float > R;
        myranker::rank(phenotype_orig[p], R, method);
		phenotype_orig[p] = R;
	}
}

void data::rankTranformGenotypes() {
	string method = "average";
	LOG.println("\nRanking genotypes");
	for (int g = 0; g < genotype_count ; g ++) {
		vector < float > R;
        myranker::rank(genotype_orig[g], R, method);
		genotype_orig[g] = R;
	}
}

void data::normalize(vector < float > & X) {
	double mean = 0.0, sum = 0.0;
	for (int s = 0; s < sample_count ; s ++) mean += X[s];
	mean /= sample_count;
	for (int s = 0; s < sample_count ; s ++) {
		X[s] -= mean;
		sum += X[s] * X[s];
	}
	sum = sqrt(sum);
	if (sum == 0) sum = 1;
	for (int s = 0; s < sample_count ; s ++) X[s] /= sum;
}

void data::normalize(vector < vector < float > > & X) {
	for (int x = 0 ; x < X.size() ; x++) {
		double mean = 0.0, sum = 0.0;
		for (int s = 0; s < sample_count ; s ++) mean += X[x][s];
		mean /= sample_count;
		for (int s = 0; s < sample_count ; s ++) {
			X[x][s] -= mean;
			sum += X[x][s] * X[x][s];
		}
		sum = sqrt(sum);
		if (sum == 0) sum = 1;
		for (int s = 0; s < sample_count ; s ++) X[x][s] /= sum;
	}
}

void data::correct(vector < float > & X, vector < float > & Y) {
	double corr = getCorrelation(X, Y);
	for (int s = 0 ; s < sample_count ; s ++) Y[s] = Y[s] - corr * X[s];
}
