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

void data::readSamplesToExclude(string file) {
	string buffer;
	LOG.println("\nReading list of samples to exclude in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) sample_exclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(sample_exclusion.size()) + " sample(s) found");
	fd.close();
}

void data::readSamplesToInclude(string file) {
	string buffer;
	LOG.println("\nReading list of samples to include in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) sample_inclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(sample_inclusion.size()) + " sample(s) found");
	fd.close();
}

void data::readGenotypesToExclude(string file) {
	string buffer;
	LOG.println("\nReading list of sites to exclude in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) genotype_exclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(genotype_exclusion.size()) + " site(s) found");
	fd.close();
}

void data::readGenotypesToInclude(string file) {
	string buffer;
	LOG.println("\nReading list of sites to include in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) genotype_inclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(genotype_inclusion.size()) + " site(s) found");
	fd.close();
}

void data::readPhenotypesToExclude(string file) {
	string buffer;
	LOG.println("\nReading list of phenotypes to exclude in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) phenotype_exclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(phenotype_exclusion.size()) + " phenotype(s) found");
	fd.close();
}

void data::readPhenotypesToInclude(string file) {
	string buffer;
	LOG.println("\nReading list of phenotypes to include in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) phenotype_inclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(phenotype_inclusion.size()) + " phenotype(s) found");
	fd.close();
}

void data::readCovariatesToExclude(string file) {
	string buffer;
	LOG.println("\nReading list of covariates to exclude in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) covariate_exclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(covariate_exclusion.size()) + " covariate(s) found");
	fd.close();
}

void data::readCovariatesToInclude(string file) {
	string buffer;
	LOG.println("\nReading list of covariates to include in [" + file + "]");
	ifile fd(file);
	while(getline(fd, buffer, '\n')) covariate_inclusion.insert(buffer);
	LOG.println("  * " + sutils::int2str(covariate_inclusion.size()) + " covariate(s) found");
	fd.close();
}
