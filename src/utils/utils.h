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

#ifndef _UTILS_H
#define _UTILS_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.14159265358979323846

#define LOW_POS_DOUBLE 1e-300
#define BIG_POS_DOUBLE 1e300
#define LOW_NEG_DOUBLE -1e-300
#define BIG_NEG_DOUBLE -1e300
#define LOW_POS_FLOAT 1e-8
#define BIG_POS_FLOAT 1e8
#define LOW_NEG_FLOAT -1e-8
#define BIG_NEG_FLOAT -1e8
#define BIG_POS_INT 1000000000
#define BIG_NEG_INT -1000000000

#define ___NA___ (0.0/0.0)

#include <string>
#include <complex>
#include <vector>
#include <queue>
#include <map>
#include <numeric>
#include <bitset>
#include <list>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;


/******************************************************/
/*                  UTILS RUNNING STATS               */
/******************************************************/

class RunningStat {
public:
	RunningStat() : m_n(0) {}

	RunningStat(vector < double > & X): m_n(0) {
		for (unsigned int e = 0 ; e < X.size() ; e ++) Push(X[e]);
	}

	RunningStat(vector < float > & X): m_n(0) {
		for (unsigned int e = 0 ; e < X.size() ; e ++) Push((double)X[e]);
	}

	void Clear() { m_n = 0; }

	void Push(double x) {
		m_n++;

		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
            m_oldM = m_newM;
            m_oldS = m_newS;
		}
	}

	int NumDataValues() const {
		return m_n;
	}

	double Mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	double Variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	double StandardDeviation() const {
		return sqrt( Variance() );
	}

	void MeanStandardDeviation (double & mean, double & sd) {
		mean = Mean();
		sd = StandardDeviation();
	}

	void MeanVariance (double & mean, double & variance) {
		mean = Mean();
		variance = Variance();
	}

private:
	int m_n;
	double m_oldM, m_newM, m_oldS, m_newS;
};


/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
namespace putils {
	double mean(vector < double > &);
	double variance(vector < double > &, double );
	double mean(vector < float > &);
	double variance(vector < float > &, double );
	bool isVariable(vector < float > &);
	bool isVariable(vector < double > &);
	void initRandom(long s);
	double getRandom();
	string getRandomID();
	int getRandom(int);
	void bootstrap(int, vector < int > &);
	long getSeed();
	void normalise(vector < double > & v);
	void random_shuffle(vector < vector < float > > &);
	int sample(vector< double > & v, double sum);
	double entropy(vector < double > & v);
	double KLdistance(vector < double > & P, vector < double > & Q);
	double qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
	int max(vector < double > & v);
	int max(vector < int > & v);
	void findUniqueSet(vector < bool > & B, vector < int > & U);
	void reorder(vector < float > &, vector < unsigned int > const & ) ;
	void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB);
	int checkDuo (int pa1, int pa2, int ca1, int ca2);
	int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2);
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
	int tokenize(string &, vector < string > &);
	int tokenize(string &, vector < string > &, int);
	int tokenize(string &, vector < string > &, string);
	string int2str(int n);
	string int2str(vector < int > & v);
	string long2str(long int n);
	string double2str(double n);
	string double2str(double n, int prc);
	string double2str(vector < double > &v);
	string double2str(vector < double > &v, int prc);
	string double2str(vector < float > &v);
	string double2str(vector < float > &v, int prc);
	string bool2str(vector<bool> & v);
	string date2str(time_t * t, string format);
	bool isNumeric(string &);
	string remove_spaces(const string &);
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
	bool isFile(string f);
	bool createFile(string f);
	void checkFile(string f, bool);
	string extensionFile(string & filename);
	void bool2binary(vector < bool > & V, ostream &fd);
	bool binary2bool(vector < bool > & V, istream & fd);
};


/******************************************************/
/*                  EXCEPTIONS                        */
/******************************************************/
class myException : public exception {
public:
   explicit myException(std::string msg) : msg_(msg) {}

   virtual ~myException() throw() {}

   virtual const char* what() const throw() {
      return msg_.c_str();
   }

private:
   std::string msg_;
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
class ifile : public bio::filtering_istream {
private:
	string file;
	ifstream fd;

public:
	ifile();
	ifile(string filename , bool binary = false);
	~ifile();
	string name();
	bool open(string filename, bool binary = false);
	bool readString(string &);
	void close();
};

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
class ofile : public bio::filtering_ostream {
private:
	string file;
	ofstream fd;

public:
	ofile();
	ofile(string filename , bool binary = false);
	~ofile();
	string name();
	bool open(string filename, bool binary = false);
	void writeString(string &);
	void close();
};

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
class lfile {
private:
	string file;
	ofstream fd;
	bool verboseC;
	bool verboseL;

public:
	lfile();
	~lfile();
	string name();
	bool open(string filename = "file.log");
	void close();
	string getPrefix();
	void muteL();
	void unmuteL();
	void muteC();
	void unmuteC();
	void print(string s);
	void printC(string s);
	void printL(string s);
	void println(string s);
	void printlnC(string s);
	void printlnL(string s);
	void warning(string s);
	void error(string s);
};

#ifdef	_GLOBAL
#define _EXTERN
#else
#define _EXTERN	 extern
#endif

_EXTERN lfile LOG;
_EXTERN time_t START_DATE;

#endif
