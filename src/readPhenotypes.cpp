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
#include "utils/tabix.hpp"


void data::readPhenotypes(string fbed) {
	string buffer;
	vector < string > str;
	int n_includedP = 0;
	int n_excludedP = 0;
	int n_includedS = 0;
	int n_excludedS = 0;
	int n_parsed = 0;
	vector < int > mappingS;

	//Initialise
	LOG.println("\nReading phenotype data in [" + fbed + "]");
	if (!futils::isFile(fbed + ".tbi")) LOG.error("index file missing [" + fbed + ".tbi]");
	Tabix fd (fbed);

	//Read samples ids
	fd.getLastHeader(buffer);
	if (buffer.size() == 0) LOG.error("No header line detected");
	sutils::tokenize(buffer, str);
    if (str.size() < 5) LOG.error("Wrong BED header format for sample ids");
	for (int t = 4 ; t < str.size() ; t ++) {
		if (checkSample(str[t])) {
			sample_id.push_back(str[t]);
			mappingS.push_back(n_includedS);
			n_includedS ++;
		} else {
			mappingS.push_back(-1);
			n_excludedS ++;
		}
	}
	sample_count = sample_id.size();

	//Read phenotypes
	if (!fd.setRegion(regionPhenotype.str())) LOG.error("Failed to get region " + regionPhenotype.str() + " in [" + fbed + "]");
	while (fd.getNextLine(buffer)) {
        if (buffer.size() == 0) continue;
		sutils::tokenize(buffer, str);
		if (str.size() < 5) LOG.error("Wrong BED file format 2 = " +sutils::int2str(str.size()));
		if (checkPhenotype(str[3])) {
			phenotype_id.push_back(str[3]);
			phenotype_chr.push_back(str[0]);
			phenotype_start.push_back(atoi(str[1].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_end.push_back(atoi(str[2].c_str()));
			phenotype_orig.push_back(vector < float > (sample_count, 0.0));
			for (int t = 4 ; t < str.size() ; t ++) {
				if (mappingS[t-4] >= 0) {
					if (str[t] == "NA") phenotype_orig.back()[mappingS[t-4]] = ___NA___;
					else if (!sutils::isNumeric(str[t])) LOG.error("Phenotype encountered is not a number, check: [" + str[t] + "]");
					else phenotype_orig.back()[mappingS[t-4]] = atof(str[t].c_str());
				}
			}
			n_includedP++;
		} else n_excludedP ++;
        n_parsed++;
		if (n_parsed % 100000 == 0) LOG.println("  * " + sutils::int2str(n_parsed) + " lines parsed");
	}

	//Finalise
	phenotype_count = phenotype_id.size();
	LOG.println("  * region = " + regionPhenotype.str());
	LOG.println("  * " + sutils::int2str(n_includedS) + " samples included");
	if (n_excludedS > 0) LOG.println("  * " + sutils::int2str(n_excludedS) + " samples excluded");
	LOG.println("  * " + sutils::int2str(n_includedP) + " phenotypes included");
	if (n_excludedP > 0) LOG.println("  * " + sutils::int2str(n_excludedP) + " phenotypes excluded");
    if (phenotype_count == 0) LOG.error("No phenotypes in this region");
}

void data::scanPhenotypes(string fbed) {
	string buffer;
	vector < string > str;
	int n_includedP = 0;
	int n_duplicateP = 0;

	//Initialise
	LOG.println("\nScanning phenotype data in [" + fbed + "]");
	if (!futils::isFile(fbed + ".tbi")) LOG.error("index file missing [" + fbed + ".tbi]");
	Tabix fd (fbed);

	//Read lines
	fd.getLastHeader(buffer);
	while (fd.getNextLine(buffer)) {
        if (buffer.size() == 0) continue;
		sutils::tokenize(buffer, str, 4);
		if (str.size() != 4) LOG.error("Wrong BED file format");
		if (checkPhenotype(str[3], false)) {
			phenotype_id.push_back(str[3]);
			phenotype_chr.push_back(str[0]);
			phenotype_start.push_back(atoi(str[1].c_str()) + 1); //convert to 1-based, tabix works on 1-based coordinates and all other files are 1-based
			phenotype_end.push_back(atoi(str[2].c_str()));
			n_includedP++;
		} else n_duplicateP ++;
	}

	//Finalise
	phenotype_count = phenotype_id.size();
	LOG.println("  * " + sutils::int2str(n_includedP) + " phenotypes");
	if (n_duplicateP > 0) LOG.error("Duplicated phenotype IDs found n=" + sutils::int2str(n_duplicateP));
    if (phenotype_count == 0) LOG.error("No phenotypes in this region");
}
