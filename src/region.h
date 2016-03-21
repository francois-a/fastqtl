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

#ifndef _REGION_H
#define _REGION_H

#define POS_MIN 0
#define POS_MAX 1000000000

class region {
public:
	string chr;
	int start;
	int end;

	region() {
		chr = "chrNA";
		start = POS_MIN;
		end = POS_MAX;
	}

	string str() {
		ostringstream s2( stringstream::out );
		s2 << chr << ":" << start << "-" << end;
		return s2.str();
	}

	bool check (string _chr, int _pos) {
		if (_chr == chr && _pos >= start && _pos < end) return true;
		else return false;
	}

	bool set(string r) {
		vector < string > chr_split, pos_split;
		sutils::tokenize(r, chr_split, ":");
		if (chr_split.size() == 2) {
			chr = chr_split[0];
			sutils::tokenize(chr_split[1], pos_split, "-");
			if (pos_split.size() == 2) {
				start = atoi(pos_split[0].c_str());
				end = atoi(pos_split[1].c_str());
				return true;
			} else return false;
		} else if (chr_split.size() == 1) {
			chr = chr_split[0];
			start = POS_MIN;
			end = POS_MAX;
			return true;
		} else return false;
	}
};

#endif
