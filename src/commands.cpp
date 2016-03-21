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

void data::writeCommands(string fout, int nChunks, int argc, char ** argv) {
	LOG.println("\nGenerate " + sutils::int2str(nChunks) + " commands in [" + fout + "]");

	vector < string > args;
	for (int a = 0 ; a < argc ; a ++) args.push_back(argv[a]);

	//map commands argument
	int idx_commands_arg = -1;
	for (int a = 0 ; a < argc ; a ++) if (args[a] == "--commands") idx_commands_arg = a;
	args.erase (args.begin() + idx_commands_arg, args.begin() + idx_commands_arg + 3);

	//map output argument
	int idx_output_arg = -1;
	string str_output_arg = "";
	for (int a = 0 ; a < argc ; a ++) {
		if (args[a] == "--out") {
			idx_output_arg = a;
			str_output_arg = string(args[a+1]);
		}
	}

	//write commands [loop over regions]
	ofile fd(fout);
	for (int c = 0 ; c < nChunks ; c ++) {
		setPhenotypeRegion(c);
		string region = getPhenotypeRegion(c);

		for (int a = 0 ; a < args.size() ; a ++ )
			if (a == idx_output_arg + 1) fd << " " << str_output_arg + "." + region;
			else fd << " " << args[a];

		fd << " --region " << region << endl;
	}
}
