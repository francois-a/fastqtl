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


void data::readInteractions(string fcov) {
	string buffer;
	vector < string > str;
	int n_sample_found = 0;

	LOG.println("\nReading interaction term in [" + fcov + "]");
	ifile fd (fcov);

	//Read interactions
	interaction_val = vector < float >(sample_count, 0.0);
	while(getline(fd, buffer)) {
        sutils::tokenize(buffer, str, "\t");
        if (str.size() != 2) LOG.error("Wrong interaction file format");

        //
        int idx_sample = -1;
        for (int i = 0 ; i < sample_count && idx_sample < 0 ; i++) if (sample_id[i] == sutils::remove_spaces(str[0])) idx_sample = i;

        //
        if (idx_sample >=0) {
        	interaction_val[idx_sample] = atof(str[1].c_str());
        	n_sample_found ++;
        }
	}
	if (n_sample_found != sample_count)
		LOG.error("Number of overlapping samples in interaction file is " + sutils::int2str(n_sample_found) + " and should be " + sutils::int2str(sample_count));

	//Finalise
	LOG.println("  * " + sutils::int2str(n_sample_found) + " samples included");
	fd.close();
}

