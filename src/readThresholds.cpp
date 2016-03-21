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


void data::readThresholds(string fres) {
	string buffer;
	vector < string > str;

	//1.0 Allocation
	phenotype_threshold = vector < double > (phenotype_count, -1);
	vector < bool > phenotype_mask = vector < bool > (phenotype_count, false);

	//2.0 Read results
	LOG.println("\nReading Qvalues in [" + fres + "]");
	ifile fdr(fres);
	while (getline(fdr, buffer)) {
		sutils::tokenize(buffer, str);
		if (str.size() < 2) LOG.error("Incorrect number of columns in [" + fres + "], observed = " + sutils::int2str(str.size())  + " expected > 2");

		int phenotype_idx = -1;
		for (int p = 0 ; p < phenotype_count && phenotype_idx < 0 ; p ++) if (phenotype_id[p] == str[0]) phenotype_idx = p;

		if (phenotype_idx >= 0) {
			phenotype_threshold[phenotype_idx] = atof(str[1].c_str());
			phenotype_mask[phenotype_idx] = true;
		}
	}
	fdr.close();

	//3.0 Make sure that each MP has a qvalue
	int n_set= 0, n_unset = 0;
	for (int p = 0 ; p < phenotype_count ; p ++) {
		if (phenotype_mask[p]) n_set ++;
		else n_unset ++;
	}
	LOG.println("  * #phenotypes set = " + sutils::int2str(n_set));
	if (n_unset > 0) LOG.error(sutils::int2str(n_unset) + " phenotypes without qvalues!");
}
