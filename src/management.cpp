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
#include <map>

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
    genotype_maf.clear();
    genotype_ma_count.clear();
    genotype_ma_samples.clear();
    genotype_ref_factor.clear();
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

    if (K < 1) LOG.error("Number of chunks needs to be > 0");

    std::map<std::string, vector<int> > chunk_sizes;
    // only used if phenotype groups are provided:
    std::map<std::string, vector<int> > group_sizes;
    int group_count = 0;
    std::map<std::string,int> num_pheno;

    if (phenotype_grp.empty()) {
        if (K > phenotype_count) LOG.error("Number of chunks (" + sutils::int2str(K) + ") is greater than the number of phenotypes (" + sutils::int2str(phenotype_count) + ")");

        // count number of phenotypes per chromosome
        for (int p=0; p<phenotype_count; ++p) {
            num_pheno[phenotype_chr[p]]++;
        }
        if (num_pheno.size() > K) LOG.error("Number of chunks needs to be >= #chr");

    } else {
        LOG.println("\nCalculating chunks based on group sizes");

        // count number of groups and group sizes for each chr
        int group_size = 1;
        // phenotype_grp: mapping of BED position to group_id, parsed in readGroups.cpp
        for (int p=1; p<phenotype_count; ++p) {
            if (phenotype_grp[p]!=phenotype_grp[p-1]) {  // new group
                group_sizes[phenotype_chr[p-1]].push_back(group_size);  // push previous group
                group_size = 1;
            } else {
                group_size += 1;
            }
        }
        // push last
        group_sizes[phenotype_chr[phenotype_count-1]].push_back(group_size);
        // count totals
        for(const auto &kv : group_sizes) {
            num_pheno[kv.first] += kv.second.size();
            group_count += kv.second.size();
        }
        if (num_pheno.size() > K) LOG.error("Number of chunks needs to be >= #chr");
        if (K > group_count) LOG.error("Number of chunks cannot exceed number of phenotype groups");
    }

    // initialize chunk size to all phenotypes
    std::map<std::string,int> chunk_size(num_pheno);

    // each chr must have at least 1 chunk
    std::map<std::string,int> num_chunks;
    for (const auto &kv : num_pheno) {
        num_chunks[kv.first] = 1;
    }

    // determine number of chunks for each chr
    for (auto i=0;i<K-num_pheno.size(); ++i) {
        // find largest chunk size
        int max_chunk = 0;
        std::string max_key;
        for(const auto &kv : chunk_size) {
            if (kv.second>max_chunk) {
                max_chunk = kv.second;
                max_key = kv.first;
            }
        }
        num_chunks[max_key] += 1;  // add chunk for this chr
        chunk_size[max_key] = ceil((double)num_pheno[max_key]/(double)num_chunks[max_key]);  // adjust chunk size
    }

    // evenly distributed chunk sizes
    for (const auto &kv : num_pheno) {  // loop chrs
        chunk_sizes[kv.first] = vector <int> (num_chunks[kv.first], (int)floor(num_pheno[kv.first]/num_chunks[kv.first]));
        int d = num_pheno[kv.first] - std::accumulate(chunk_sizes[kv.first].begin(), chunk_sizes[kv.first].end(), 0);
        for (int i=0;i<d;++i) {
            chunk_sizes[kv.first][i] += 1;
        }
    }

    // now, split according to chr and chunk size
    phenotype_cluster = vector < vector < int > > (1, vector < int > (1, 0));  // first phenotype already stored (index 0)
    int chunk = 0;
    int chunk_index = 1;  // next position in chunk

    if (phenotype_grp.empty()) {
        for (int p=1; p<phenotype_count; ++p) {
            if (phenotype_chr[p]!=phenotype_chr[p-1]) {  // new chr --> new chunk
                chunk = 0;
                chunk_index = 1;
                phenotype_cluster.push_back(vector < int > (1, p));
            } else if (chunk_index==chunk_sizes[phenotype_chr[p]][chunk]) {  // new chunk
                chunk += 1;
                chunk_index = 1;
                phenotype_cluster.push_back(vector < int > (1, p));
            } else {  // add to current chunk
                phenotype_cluster.back().push_back(p);
                chunk_index += 1; // next index in current chunk
            }
        }
    } else {
        int current_group = 0;  // within chunk
        int current_chunk = 0;  // within chr
        int groups_parsed = 1;
        for (int p=1; p<phenotype_count; ++p) {
            if (phenotype_grp[p]!=phenotype_grp[p-1]) {  // new group
                current_group += 1;
                groups_parsed += 1;
                if (phenotype_chr[p]!=phenotype_chr[p-1]) {  // new chr -> new chunk, new group
                    phenotype_cluster.push_back(vector < int > (1, p));
                    current_group = 0;
                    current_chunk = 0;
                } else if (current_group==chunk_sizes[phenotype_chr[p]][current_chunk]) {  // previous group was last in chunk
                    phenotype_cluster.push_back(vector < int > (1, p));
                    current_group = 0;
                    current_chunk += 1;
                } else {  // add to current chunk
                    phenotype_cluster.back().push_back(p);
                }
            } else {  // same group
                phenotype_cluster.back().push_back(p);
            }
        }
        LOG.println("  * groups parsed = "+ sutils::int2str(groups_parsed) );
    }

    if (phenotype_cluster.size() != K)
        LOG.error("Chunks ("+sutils::int2str(phenotype_cluster.size())+") do not match input ("+sutils::int2str(K)+")");
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
