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

int main(int argc, char ** argv) {
	data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------

	bpo::options_description opt_basic ("\33[33mBasic options\33[0m");
	opt_basic.add_options()
		("help", "Produces this help")
		("silent", "Silent mode on terminal")
		("seed", bpo::value< int >()->default_value(time(NULL)), "Random number seed. Useful to replicate runs.");

	bpo::options_description opt_files ("\33[33mInput/Output files\33[0m");
	opt_files.add_options()
		("log,L", bpo::value< string >()->default_value("fastQTL_date_time_UUID.log"), "Screen output is copied in this file.")
		("vcf,V", bpo::value< string >(), "Genotypes in VCF format.")
		("bed,B", bpo::value< string >(), "Phenotypes in BED format.")
		("cov,C", bpo::value< string >(), "Covariates in TXT format.")
		("grp,G", bpo::value< string >(), "Phenotype groups in TXT format.")
		("mtx,M", bpo::value< string >(), "Output folder for matrices.")
		("out,O", bpo::value< string >(), "Output file.");

	bpo::options_description opt_exclusion ("\33[33mExclusion/Inclusion files\33[0m");
	opt_exclusion.add_options()
		("exclude-samples", bpo::value< string >(), "List of samples to exclude.")
		("include-samples", bpo::value< string >(), "List of samples to include.")
		("exclude-sites", bpo::value< string >(), "List of sites to exclude.")
		("include-sites", bpo::value< string >(), "List of sites to include.")
		("exclude-phenotypes", bpo::value< string >(), "List of phenotypes to exclude.")
		("include-phenotypes", bpo::value< string >(), "List of phenotypes to include.")
		("exclude-covariates", bpo::value< string >(), "List of covariates to exclude.")
		("include-covariates", bpo::value< string >(), "List of covariates to include.");

	bpo::options_description opt_parameters ("\33[33mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform the phenotypes.")
		("window,W", bpo::value< double >()->default_value(1e6), "Cis-window size.")
		("threshold,T", bpo::value< double >()->default_value(1.0), "P-value threshold used in nominal pass of association")
		("maf-threshold", bpo::value< double >()->default_value(0.0), "Minor allele frequency threshold used when parsing genotypes")
		("ma-sample-threshold", bpo::value< int >()->default_value(0), "Minimum number of samples carrying the minor allele; used when parsing genotypes")
		("global-af-threshold", bpo::value< double >()->default_value(0.0), "AF threshold for all samples in VCF (used to filter AF in INFO field)");

	bpo::options_description opt_modes ("\33[33mModes\33[0m");
	opt_modes.add_options()
		("permute,P", bpo::value< vector < int > >()->multitoken(), "Permutation pass to calculate corrected p-values for molecular phenotypes.")
		("psequence", bpo::value< string >(), "Permutation sequence.")
		("map", bpo::value< string >(), "Map best QTL candidates per molecular phenotype.")
		("map-full", "Scan full cis-window to discover independent signals.")
		("interaction", bpo::value< string >(), "Test for interactions with variable specified in file.")
		("report-best-only", bpo::bool_switch()->default_value(false), "Report best variant only (nominal mode)");

	bpo::options_description opt_parallel ("\33[33mParallelization\33[0m");
	opt_parallel.add_options()
		("chunk,K", bpo::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed")
		("commands", bpo::value< vector < string > >()->multitoken(), "Generates all commands")
		("region,R", bpo::value< string >(), "Region of interest.");

	bpo::options_description descriptions;
	descriptions.add(opt_basic).add(opt_files).add(opt_exclusion).add(opt_parameters).add(opt_modes).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	bpo::variables_map options;
	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing command line :" << string(e.what()) << endl;
		exit(0);
	}

	//-----------------------
	// 3. PRINT HEADER/HELP
	//-----------------------
	if (! options.count("silent")) {
		cout << endl;
		cout << "\33[33mF\33[0mast \33[33mQTL\33[0m" << endl;
		cout << "  * Authors : Olivier DELANEAU, Halit ONGEN, Alfonso BUIL & Manolis DERMITZAKIS" << endl;
		cout << "  * Contact : olivier.delaneau@gmail.com" << endl;
		cout << "  * Webpage : http://fastqtl.sourceforge.net/" << endl;
		cout << "  * Version : v2.184_gtex" << endl;
		if (options.count("help")) { cout << descriptions<< endl; exit(1); }
	}

	//--------------
	// 4. LOG FILE
	//--------------
	struct timeval start_time, stop_time;
	gettimeofday(&start_time, 0);
	START_DATE = time(0);
	//localtime(&START_DATE);
	//string logfile = "fastQTL_" + sutils::date2str(&START_DATE, "%d%m%Y_%Hh%Mm%Ss") + "_" + putils::getRandomID() + ".log";
	if (!options["log"].defaulted()) {
		if (!LOG.open(options["log"].as < string > ())) {
			cerr << "Impossible to open log file[" << options["log"].as < string > () << "] check writing permissions!" << endl;
			exit(1);
		}
	} else LOG.muteL();
	if (options.count("silent")) LOG.muteC();

	//------------------------
	// 5. OPTIONS COMBINATIONS
	//------------------------
	if (!options.count("vcf")) LOG.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!options.count("bed")) LOG.error("Phenotype data needs to be specified with --bed [file.bed]");
	if (!options.count("out")) LOG.error("Output needs to be specified with --out [file.out]");

	int nParallel = options.count("chunk") + options.count("commands") + options.count("region");
	if (nParallel != 1) LOG.error("Please, specify one of these options [--region, --chunk, --commands]");

	int nMode = options.count("permute") + options.count("map");
	if (nMode > 1) LOG.error("Please, specify only one of these options [--permute, --map]");

	if (options.count("mtx"))  LOG.println("\n *** cis-eQTL matrices will be written to folder ["+options["mtx"].as <string> () +"]\n");

	//---------------
	// 6. CHECK FILES
	//---------------
	if (!futils::isFile(options["vcf"].as < string > ())) LOG.error(options["vcf"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (!futils::isFile(options["bed"].as < string > ())) LOG.error(options["bed"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("cov") && !futils::isFile(options["cov"].as < string > ())) LOG.error(options["cov"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("interaction") && !futils::isFile(options["interaction"].as < string > ())) LOG.error(options["interaction"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("grp") && !futils::isFile(options["grp"].as < string > ())) LOG.error(options["grp"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("map") && !futils::isFile(options["map"].as < string > ())) LOG.error(options["map"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (!futils::createFile(options["out"].as < string > ())) LOG.error(options["out"].as < string > () + " is impossible to create, check writing permissions");

	//-----------------------------------
	// 6. CHECK INCLUSION/EXCLUSION FILES
	//-----------------------------------
	if (options.count("exclude-samples") && !futils::isFile(options["exclude-samples"].as < string > ())) LOG.error(options["exclude-samples"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-samples") && !futils::isFile(options["include-samples"].as < string > ())) LOG.error(options["include-samples"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-sites") && !futils::isFile(options["exclude-sites"].as < string > ())) LOG.error(options["exclude-sites"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-sites") && !futils::isFile(options["include-sites"].as < string > ())) LOG.error(options["include-sites"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-phenotypes") && !futils::isFile(options["exclude-phenotypes"].as < string > ())) LOG.error(options["exclude-phenotypes"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-phenotypes") && !futils::isFile(options["include-phenotypes"].as < string > ())) LOG.error(options["include-phenotypes"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("exclude-covariates") && !futils::isFile(options["exclude-covariates"].as < string > ())) LOG.error(options["exclude-covariates"].as < string > () + " is impossible to open, check file existence or reading permissions");
	if (options.count("include-covariates") && !futils::isFile(options["include-covariates"].as < string > ())) LOG.error(options["include-covariates"].as < string > () + " is impossible to open, check file existence or reading permissions");

	//----------------------------
	// 7. CHECK METHODS/PARAMETERS
	//----------------------------
	if (options.count("interaction")) {
		if (options.count("permute")) {
			LOG.println("\nPerform permutation based interaction analysis (used to calculate corrected p-values for MPs)");
			vector < int > nPerm = options["permute"].as < vector < int > > ();
			if (nPerm.size() != 1) LOG.error("Interactions only work with a fixed number of permutations!");
			else {
				if (nPerm[0] < 50) LOG.warning("Permutation number seems to be low, check parameters");
				LOG.println("  * Perform " + sutils::int2str(nPerm[0]) + " permutations");
			}
		} else LOG.error("Interactions only work with permutations!");
		LOG.println("  * Test interaction with term from [" + options["interaction"].as < string > () + "]");
	} else if (options.count("permute")) {
		LOG.println("\nPerform permutation based analysis (used to calculate corrected p-values for MPs)");
		vector < int > nPerm = options["permute"].as < vector < int > > ();
		if (nPerm.size() > 3 || nPerm.size() < 1) LOG.error ("Option --permute takes 1, 2 or 3 arguments");
		if (nPerm.size() == 1) {
			if (nPerm[0] <= 0) LOG.error("Permutation number needs to be positive integer");
			if (nPerm[0] < 50) LOG.warning("Permutation number seems to be low, check parameters");
			LOG.println("  * Perform " + sutils::int2str(nPerm[0]) + " permutations");
		} else if (nPerm.size() == 2) {
			if (nPerm[0] <= 0 || nPerm[1] <= 0) LOG.error("Permutation number needs to be positive");
			if (nPerm[1] <= nPerm[0]) LOG.error("For adaptive permutation scheme, arg1 needs to be smaller than arg2!");
			LOG.println("  * Perform between " + sutils::int2str(nPerm[0]) + " and " + sutils::int2str(nPerm[1]) + " permutations");
		} else {
			if (nPerm[0] <= 0 || nPerm[1] <= 0 || nPerm[2] <= 0) LOG.error("Permutation number needs to be positive");
			if (nPerm[2] <= nPerm[0]) LOG.error("For adaptive permutation scheme, arg1 needs to be smaller than arg3!");
			if (nPerm[0] <= nPerm[1]) LOG.error("For adaptive permutation scheme, arg2 needs to be smaller than arg1!");
			LOG.println("  * Perform between " + sutils::int2str(nPerm[0]) + " and " + sutils::int2str(nPerm[2]) + " permutations and stop when " + sutils::int2str(nPerm[1]) + " best associations are found");
		}
		if (options.count("grp")) LOG.println("  * Using MP groups from [" + options["grp"].as < string > () + "]");
	} else if (options.count("map")) {
		LOG.println("\nPerform conditional based analysis (used to map significant QTLs for MPs");
		LOG.println("  * Using per MP p-value threshold from [" + options["map"].as < string > () + "]");
		if (options.count("map-full")) LOG.println("  * Scanning all variants in cis and not only nominally significant ones");
	} else {
		LOG.println("\nPerform nominal analysis (used to get raw p-values of association)");
		double threshold = options["threshold"].as < double > ();
		if (threshold <= 0.0 || threshold > 1.0) LOG.error("Incorrect --threshold value  :  0 < X <= 1");
		LOG.println("  * Using p-value threshold = " + sutils::double2str(threshold, 10));
	}

	if (options["seed"].as < int > () < 0) LOG.error("Random number generator needs a positive seed value");
	else srand(options["seed"].as < int > ());
	LOG.println("  * Random number generator is seeded with " + sutils::int2str(options["seed"].as < int > ()));

	if (options["window"].as < double > () <= 0) LOG.error ("Incorrect value for option --window (null or negative value)");
	if (options["window"].as < double > () > 1e9) LOG.error ("Cis-window cannot be larger than 1e9bp");
	LOG.println("  * Considering variants within " + sutils::double2str(options["window"].as < double > ()) + " bp of the MPs");
	D.cis_window = options["window"].as < double > ();
	
	D.maf_threshold = options["maf-threshold"].as < double > ();
	if (D.maf_threshold < 0.0 || D.maf_threshold >= 0.5) LOG.error("Incorrect --maf-threshold value  :  0 <= X < 0.5");
	LOG.println("  * Using minor allele frequency threshold = " + sutils::double2str(D.maf_threshold, 4));
	
	D.ma_sample_threshold = options["ma-sample-threshold"].as < int > ();
	if (D.ma_sample_threshold < 0) LOG.error("Incorrect --ma-sample-threshold  :  0 <= X");
	LOG.println("  * Using minor allele sample threshold = " + sutils::int2str(D.ma_sample_threshold));
	
	D.global_af_threshold = options["global-af-threshold"].as < double > ();
	if (D.global_af_threshold < 0.0 || D.global_af_threshold >= 0.5) LOG.error("Incorrect --global-af-threshold value  :  0 <= X < 0.5");
	LOG.println("  * Using INFO field AF threshold = " + sutils::double2str(D.global_af_threshold, 4));

	if (options.count("chunk")) {
		vector < int > nChunk = options["chunk"].as < vector < int > > ();
		if (nChunk.size() != 2) LOG.error ("--chunk needs 2 integer arguments");
		if (nChunk[0] > nChunk[1]) LOG.error ("arg1 for --chunk needs to be smaller or equal to arg2");
		LOG.println ("  * Chunk processed " + sutils::int2str(nChunk[0]) + " / " + sutils::int2str(nChunk[1]));
	} else if (options.count("commands")) {
		vector < string > nCommands = options["commands"].as < vector < string > > ();
		if (nCommands.size() != 2) LOG.error ("--commands needs 2 arguments");
		LOG.println ("  * " + nCommands[0] + " commands output in [" + nCommands[1] +"]");
	} else LOG.println ("  * Focus on region [" + options["region"].as < string > () +"]");

	//--------------------------------
	// 7. READ EXCLUDE / INCLUDE FILES
	//--------------------------------
	if (options.count("exclude-samples")) D.readSamplesToExclude(options["exclude-samples"].as < string > ());
	if (options.count("include-samples")) D.readSamplesToInclude(options["include-samples"].as < string > ());
	if (options.count("exclude-sites")) D.readGenotypesToExclude(options["exclude-sites"].as < string > ());
	if (options.count("include-sites")) D.readGenotypesToInclude(options["include-sites"].as < string > ());
	if (options.count("exclude-phenotypes")) D.readPhenotypesToExclude(options["exclude-phenotypes"].as < string > ());
	if (options.count("include-phenotypes")) D.readPhenotypesToInclude(options["include-phenotypes"].as < string > ());
	if (options.count("exclude-covariates")) D.readCovariatesToExclude(options["exclude-covariates"].as < string > ());
	if (options.count("include-covariates")) D.readCovariatesToInclude(options["include-covariates"].as < string > ());

	if (options.count("commands")) {
		//---------------------
		// 8. GENERATE COMMANDS
		//---------------------
		int nChunks = atoi(options["commands"].as < vector < string > > ()[0].c_str());
		D.scanPhenotypes(options["bed"].as < string > ());
		D.clusterizePhenotypes(nChunks);
		D.writeCommands(options["commands"].as < vector < string > > ()[1], nChunks, argc, argv);
	} else {
		//--------------
		// 9. SET REGION
		//--------------
		if (options.count("chunk")) {
			D.scanPhenotypes(options["bed"].as < string > ());
			D.clusterizePhenotypes(options["chunk"].as < vector < int > > ()[1]);
			D.setPhenotypeRegion(options["chunk"].as < vector < int > > ()[0] - 1);
			D.clear();
		} else 	if (!D.setPhenotypeRegion(options["region"].as < string > ())) LOG.error("Impossible to interpret region [" + options["region"].as < string > () + "]");
		D.deduceGenotypeRegion(options["window"].as < double > ());

		//---------------
		// 10. READ FILES
		//---------------
		D.readPhenotypes(options["bed"].as < string > ());
		D.readGenotypesVCF(options["vcf"].as < string > ());
		if (options.count("cov")) D.readCovariates(options["cov"].as < string > ());
		if (options.count("map")) D.readThresholds(options["map"].as < string > ());
		if (options.count("grp")) D.readGroups(options["grp"].as < string > ());
		if (options.count("interaction")) D.readInteractions(options["interaction"].as < string > ());

		//------------------------
		// 11. INITIALIZE ANALYSIS
		//------------------------
		D.imputeGenotypes();
		D.imputePhenotypes();
		if (options.count("normal")) D.normalTranformPhenotypes();
		D.initResidualizer();

		//-----------------
		// 12. RUN ANALYSIS
		//-----------------
		if (options.count("interaction"))
			D.runPermutationInteraction(options["out"].as < string > (), options["permute"].as < vector < int > > ()[0]);
		else if (options.count("permute") && options.count("grp"))
			D.runPermutationPerGroup(options["out"].as < string > (), options["permute"].as < vector < int > > ());
		else if (options.count("permute")) {
			if (options.count("psequence")) D.runPermutation(options["out"].as < string > (), options["psequence"].as < string > ());
			else D.runPermutation(options["out"].as < string > (), options["permute"].as < vector < int > > ());
		} else if (options.count("map")) {
			D.runMapping(options["out"].as < string > (), options.count("map-full"));
		} else if (options["report-best-only"].as<bool>()) {
			D.runNominalBest(options["out"].as < string > ());
		} else if (options.count("mtx")) {
			D.runNominalOutputMatrices(options["mtx"].as < string > (), options["out"].as < string > (), options["threshold"].as < double > ());
		} else {
			D.runNominal(options["out"].as < string > (), options["threshold"].as < double > ());
		}
	}

	//----------------
	// 13. TERMINATION
	//----------------
	D.clear();
	gettimeofday(&stop_time, 0);
	int n_seconds = (int)floor(stop_time.tv_sec - start_time.tv_sec);
	LOG.println("\nRunning time: " + sutils::int2str(n_seconds) + " seconds");
	if (!options["log"].defaulted()) LOG.close();
}
