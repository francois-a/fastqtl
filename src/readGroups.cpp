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


void data::readGroups(string fgrp) {
    string buffer;
    vector < string > str;

    phenotype_grp = vector < string > (phenotype_count, "");
    vector < bool > phenotype_mask = vector < bool > (phenotype_count, false);

    // map phenotype IDs to index
    std::map<std::string,int> phenotype_index_map;
    for (int p=0; p<phenotype_count; ++p) {
        phenotype_index_map[phenotype_id[p]] = p;
    }

    //Initialise
    LOG.println("\nReading phenotype groups in [" + fgrp + "]");
    ifile fdg (fgrp);
    while (getline(fdg, buffer)) {
        sutils::tokenize(buffer, str);

        if (str.size() < 2) LOG.error("Incorrect number of columns in [" + fgrp + "], observed = " + sutils::int2str(str.size())  + " expected == 2");

        // if phenotype was present in input BED, add group ID
        if (phenotype_index_map.count(str[0]) > 0) {
            int phenotype_idx = phenotype_index_map[str[0]];
            phenotype_grp[phenotype_idx] = str[1];
            phenotype_mask[phenotype_idx] = true;
        }
    }
    fdg.close();

    //3.0 Make sure that each MP has a value
    int n_set= 0, n_unset = 0;
    for (int p = 0 ; p < phenotype_count ; p ++) {
        if (phenotype_mask[p]) n_set ++;
        else n_unset ++;
    }
    LOG.println("  * #phenotypes set = " + sutils::int2str(n_set));
    if (n_unset > 0) LOG.error(sutils::int2str(n_unset) + " phenotypes without groups!");
}
