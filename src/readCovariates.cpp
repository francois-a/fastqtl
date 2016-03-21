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


void data::readCovariates(string fcov) {
	string buffer;
	vector < string > str;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_includedC = 0;
	int n_excludedC = 0;
	int n_missingS = 0;
	vector < int > mappingS;

	LOG.println("\nReading covariates in [" + fcov + "]");
	ifile fd (fcov);

	//Read samples
	getline(fd, buffer);
	if (buffer.size() == 0) LOG.error("No header line detected!");
	sutils::tokenize(buffer, str, "\t");
	for (int t = 1 ; t < str.size() ; t ++) {
		if (checkSample(str[t], false)) {
			int idx_sample = -1;
			for (int i = 0 ; i < sample_count && idx_sample < 0 ; i++) if (sample_id[i] == sutils::remove_spaces(str[t])) idx_sample = i;
			mappingS.push_back(idx_sample);
			if (idx_sample >= 0) n_includedS ++;
			else n_missingS ++;
        } else {
            mappingS.push_back(-1);
            n_excludedS ++;
        }
	}
	if (n_includedS != sample_count) LOG.error("Number of overlapping samples in covariate file is " + sutils::int2str(n_includedS) + " and should be " + sutils::int2str(sample_count));

	//Read covariates
	//int n_type0 = 0, n_type1 = 0;
	while(getline(fd, buffer)) {
        if (buffer.size() == 0) continue;
		sutils::tokenize(buffer, str, "\t");
		if (str.size() < 2) LOG.error("Wrong genotype covariate file format");
		if (checkCovariate(str[0])) {
			covariate_val.push_back(vector < string > (sample_count, "0"));
			for (int t = 1 ; t < str.size() ; t ++) {
				if (mappingS[t-1] >= 0) {
                    covariate_val.back()[mappingS[t-1]] = str[t];
				}
			}
            n_includedC ++;
		} else n_excludedC ++;
	}

	//Finalise
	covariate_count = n_includedC;
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	if (n_missingS > 0) LOG.println("  * " + sutils::int2str(n_missingS) + " samples excluded without phenotype data");
	LOG.println("  * " + sutils::int2str(n_includedC) + " covariate(s) included");
	if (n_excludedC > 0) LOG.println("  * " + sutils::int2str(n_excludedC) + " covariate(s) excluded");
	fd.close();
}

