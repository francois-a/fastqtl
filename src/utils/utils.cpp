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

#include "utils.h"

long seed = -123456789;
pthread_mutex_t mutex_rng;

namespace putils {
	void initRandom(long s) {
		seed = - s;
		pthread_mutex_init(&mutex_rng, NULL);
	}

	double mean(vector < double > & X) {
		double mean = 0.0;
		for (int x = 0 ; x < X.size() ; x ++) mean += X[x];
		mean /= X.size();
		return mean;
	}

	double variance(vector < double > & X, double mean) {
		double variance = 0.0;
		for (int x = 0 ; x < X.size() ; x++) variance += (X[x] - mean) * (X[x] - mean);
		variance /= (X.size() - 1);
		return variance;
	}

	double mean(vector < float > & X) {
		double mean = 0.0;
		for (int x = 0 ; x < X.size() ; x ++) mean += X[x];
		mean /= X.size();
		return mean;
	}

	double variance(vector < float > & X, double mean) {
		double variance = 0.0;
		for (int x = 0 ; x < X.size() ; x++) variance += (X[x] - mean) * (X[x] - mean);
		variance /= (X.size() - 1);
		return variance;
	}

	bool isVariable(vector < float > & v) {
		for (int i = 1 ; i < v.size(); i++) if (v[i] != v[i-1]) return true;
		return false;
	}

	bool isVariable(vector < double > & v) {
		for (int i = 1 ; i < v.size(); i++) if (v[i] != v[i-1]) return true;
		return false;
	}

	double getRandom() {
		pthread_mutex_lock(&mutex_rng);
		int j;
		long k;
		static long iy=0;
		static long iv[NTAB];
		double temp;
		if (seed <= 0 || !iy) {
			if (-(seed) < 1) seed=1;
			else seed = -(seed);
			for (j=NTAB+7;j>=0;j--) {
				k=(seed)/IQ;
				seed=IA*(seed-k*IQ)-IR*k;
				if (seed < 0) seed += IM;
				if (j < NTAB) iv[j] = seed;
			}
			iy=iv[0];
		}
		k=(seed)/IQ;
		seed=IA*(seed-k*IQ)-IR*k;
		if (seed < 0) seed += IM;
		j=iy/NDIV;
		iy=iv[j];
		iv[j] = seed;
		temp=AM*iy;
		pthread_mutex_unlock(&mutex_rng);
		if (temp > RNMX) return RNMX;
		else return temp;
	}

	string getRandomID(){
		return boost::lexical_cast<string>(boost::uuids::random_generator()());
	}

	long getSeed() {
		return seed;
	}

	int getRandom(int n) {
		return (int)floor(getRandom() * n);
	}

	void bootstrap(int N, vector < int > & sample) {
		sample = vector < int > (N, 0);
		for (int e = 0 ; e < N ; e ++) sample[e] = getRandom(N);
	}

	void normalise(vector < double > & v) {
		double sum = 0.0;
		for (int i = 0 ; i < v.size() ; i++) sum += v[i];
		if (sum != 0.0) for (int i = 0 ; i < v.size() ; i++) v[i] /= sum;
	}

	int sample(vector< double > & v, double sum) {
		double csum = v[0];
		double u = getRandom() * sum;
		for (int i = 0; i < v.size() - 1; ++i) {
			if ( u < csum ) return i;
			csum += v[i+1];
		}
		return v.size() - 1;
	}

	void random_shuffle(vector < vector < float > > & X) {
		int n = X[0].size();
		vector < int > O;
		for (int i = 0; i < n; i++) O.push_back(i);
		random_shuffle(O.begin(), O.end());
		vector < vector < float > > Xtmp = X;
		for (int c = 0 ; c < X.size() ; c++) for (int e = 0 ; e < n ; e++) X[c][e] = Xtmp[c][O[e]];
	}

	double entropy(vector < double > & v) {
		double e = 0.0;
		for (int i = 0 ; i < v.size() ; i ++) {
			if (v[i] > 0.0) e += v[i] * log(v[i]);
		}
		return e;
	}

	double KLdistance(vector < double > & P, vector < double > & Q) {
		assert(P.size() == Q.size());
		double d = 0.0;
		for (int i = 0 ; i < Q.size() ; i ++) {
			if (Q[i] > 0.0 && P[i] > 0.0) d += P[i] * (log(P[i]) - log(Q[i]));
		}
		return d;
	}

	double qnorm(double p, double mu, double sigma, int lower_tail, int log_p) {
		double q, r, val;
		q = p - 0.5;
		if (fabs(q) <= .425) {
			r = .180625 - q * q;
			val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r + 45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r + 133.14166789178437745) * r + 3.387132872796366608) / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r + 39307.89580009271061) * r + 21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
		} else { /* closer than 0.075 from {0,1} boundary */
			if (q > 0) r = 1 - p;
			else r = p;

			r = sqrt(- ((log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : log(r)));

			if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
				r += -1.6;
	            val = (((((((r * 7.7454501427834140764e-4 + .0227238449892691845833) * r + .24178072517745061177) * r + 1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r + 4.6303378461565452959) * r + 1.42343711074968357734) / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + .0151986665636164571966) * r + .14810397642748007459) * r + .68976733498510000455) * r + 1.6763848301838038494) * r + 2.05319162663775882187) * r + 1.);
			} else { /* very close to  0 or 1 */
				r += -5.;
				val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + .0012426609473880784386) * r + .026532189526576123093) * r + .29656057182850489123) * r + 1.7848265399172913358) * r + 5.4637849111641143699) * r + 6.6579046435011037772) / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7)* r + 1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r + .0148753612908506148525) * r + .13692988092273580531) * r + .59983220655588793769) * r + 1.);
			}
			if (q < 0.0) val = -val;
		}
		return mu + sigma * val;
	}
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
	int max(vector < double > & v) {
		double max = -1e300;
		int index_max = 0;
		for (int i = 0 ; i < v.size() ; i ++)
			if (v[i] > max) {
				max = v[i];
				index_max = i;
				}
		return index_max;
	}

	int max(vector < int > & v) {
		int max = -1000000000;
		int index_max = 0;
		for (int i = 0 ; i < v.size() ; i ++)
		if (v[i] > max) {
			max = v[i];
			index_max = i;
		}
		return index_max;
	}

	void reorder(vector < float > & v, vector < unsigned int > const & order )  {
		vector < float > vtmp = v;
		for (int e = 0 ; e < order.size() ; e ++) v[e] = vtmp[order[e]];
	}

	void findUniqueSet(vector < bool > & B, vector < int > & U) {
		U.clear();
		for (int b = 0 ; b < B.size() ; b ++)
			if (B[b]) U.push_back(b);
		//sort( U.begin(), U.end() );
		//U.erase(unique( U.begin(), U.end() ), U.end());
	}

	void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB) {
		if (B.size() < 2 * min || B.size() == 2) {
			BB.push_back(B);
			return;
		}

		int l = putils::getRandom(B.size() - 2 * min) + min;
		vector < vector < int > > LB = vector < vector < int > > (B.begin(), B.begin() + l + 1);
		vector < vector < int > > RB = vector < vector < int > > (B.begin() + l - 1, B.end());
		vector < vector < vector < int > > > LBB;
		vector < vector < vector < int > > > RBB;
		decompose(min, LB, LBB);
		decompose(min, RB, RBB);
		BB = LBB;
		BB.insert(BB.end(), RBB.begin(), RBB.end());
	}

	int checkDuo (int pa1, int pa2, int ca1, int ca2) {
		if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return -1; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return 1; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 1; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 1; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return 1; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return -1; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
		if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 0; }
		return 0;
	}

	int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2) {
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return -1; }
		if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 0; }
		return 0;
	}
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
	int tokenize (string & str, vector < string > & tokens) {
		tokens.clear();
		string::size_type p_last = str.find_first_not_of(" 	", 0);
		string::size_type p_curr = str.find_first_of(" 	", p_last);
		while (string::npos != p_curr || string::npos != p_last) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of("	 ", p_curr);
			p_curr = str.find_first_of("	 ", p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r')
			tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	int tokenize(string & str, vector < string > & tokens, int n_max_tokens) {
		tokens.clear();
		string::size_type p_last = str.find_first_not_of(" 	", 0);
		string::size_type p_curr = str.find_first_of(" 	", p_last);
		while ((string::npos != p_curr || string::npos != p_last) && tokens.size() < n_max_tokens) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of("	 ", p_curr);
			p_curr = str.find_first_of("	 ", p_last);
		}
		return tokens.size();
	}

	int tokenize (string & str, vector < string > & tokens, string sep) {
		tokens.clear();
		string::size_type p_last = str.find_first_not_of(sep, 0);
		string::size_type p_curr = str.find_first_of(sep, p_last);
		while (string::npos != p_curr || string::npos != p_last) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of(sep, p_curr);
			p_curr = str.find_first_of(sep, p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r')
			tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	bool isNumeric(string & str) {
		float n;
		std::istringstream in(str);
		if (!(in >> n)) return false;
		return true;
	}

	string int2str(int n) {
		ostringstream s2( stringstream::out );
		s2 << n;
		return s2.str();
	}

	string int2str(vector < int > & v) {
		ostringstream s2( stringstream::out );
		for (int l = 0 ; l < v.size() ; l++) {
			s2 << " " << v[l] ;
		}
		return s2.str();
	}

	string long2str(long int n) {
		ostringstream s2( stringstream::out );
		s2 << n;
		return s2.str();
	}

	string double2str(double n) {
		ostringstream s2;
		s2 << n;
		return s2.str();
	}

	string double2str(double n, int prc) {
		ostringstream s2;
		s2 << setiosflags( ios::fixed );
		if ( prc > 0 ) s2.precision(prc);
		s2 << n;
		return s2.str();
	}

	string double2str(vector < double > &v) {
		ostringstream s2;
		for (int l = 0 ; l < v.size() ; l++) s2 << v[l] << ((l!=v.size()-1)?" ":"");
		return s2.str();
	}

	string double2str(vector < double > &v, int prc) {
		ostringstream s2;
		s2 << setiosflags( ios::fixed );
		if ( prc >= 0 ) s2.precision(prc);
		for (int l = 0 ; l < v.size() ; l++) s2 << v[l] << ((l!=v.size()-1)?" ":"");
		return s2.str();
	}

	string double2str(vector < float > &v) {
		ostringstream s2;
		for (int l = 0 ; l < v.size() ; l++) s2 << v[l] << ((l!=v.size()-1)?" ":"");
		return s2.str();
	}

	string double2str(vector < float > &v, int prc) {
		ostringstream s2;
		s2 << setiosflags( ios::fixed );
		if ( prc >= 0 ) s2.precision(prc);
		for (int l = 0 ; l < v.size() ; l++) s2 << v[l] << ((l!=v.size()-1)?" ":"");
		return s2.str();
	}

	string bool2str(vector<bool> & v) {
		ostringstream s2( stringstream::out );
		for (int l = 0 ; l < v.size() ; l++) {
			if (v[l]) s2 << "1";
			else s2 << "0";
		}
		return s2.str();
	}

	string date2str(time_t * t, string format) {
		struct tm * timeinfo = localtime(t);
		char buffer[128];
		strftime(buffer, 128, format.c_str(), timeinfo);
		ostringstream s2( stringstream::out );
		s2 << buffer;
		return s2.str();
	}

    string remove_spaces(const string &s) {
             int last = s.size() - 1;
             while (last >= 0 && s[last] == ' ') --last;
             return s.substr(0, last + 1);
     }
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
	bool isFile(string f) {
		ifile inp;
		inp.open(f.c_str(), ifstream::in);
		if(inp.fail()) {
			inp.clear(ios::failbit);
			inp.close();
			return false;
		}
		inp.close();
		return true;
	}

	bool createFile(string f) {
		ofile out;
		//out.open(f.c_str(), ofstream::out);
		out.open(f.c_str());
		if (out.fail()) {
			out.clear(ios::failbit);
			out.close();
			return false;
		} else out << "";
		out.close();
		return true;
	}

	void checkFile(string f, bool readmode) {
		if (readmode && !isFile(f)) LOG.error(f + " is impossible to open, check file existence or reading permissions");
		if (!readmode && !createFile(f)) LOG.error(f + " is impossible to open, check writing permissions");
	}

	string extensionFile(string & filename) {
		if (filename.find_last_of(".") != string::npos)
			return filename.substr(filename.find_last_of(".") + 1);
		return "";
	}


	void bool2binary(vector < bool > & V, ostream &fd) {
		int nb = V.size();
		fd.write((char*)&nb, sizeof(int));
		int cpt_byte = 0;
		int cpt_bin = 0;
		int nb_byte = (int)ceil( (V.size() * 1.0) / 8);
		while (cpt_byte < nb_byte) {
			bitset<8> byte_bitset;
			for (int l = 0; l < 8 && cpt_bin < V.size() ; l++) {
				byte_bitset[l] = V[cpt_bin];
				cpt_bin ++;
			}
			char byte_char = (char)byte_bitset.to_ulong();
			fd.write(&byte_char, 1);
			cpt_byte++;
		}
	}

	bool binary2bool(vector < bool > & V, istream & fd) {
		int nb;
		fd.read((char*)&nb, sizeof(int));
		if (!fd) return false;
		int cpt_byte = 0;
		int cpt_bin = 0;
		int nb_byte = (int)ceil( (nb * 1.0) / 8);
		V = vector < bool >(nb);
		while (cpt_byte < nb_byte) {
			char byte_char;
			fd.read(&byte_char, 1);
			if (!fd) return false;
			bitset<8> byte_bitset = byte_char;
			for (int l = 0; l < 8 && cpt_bin < nb ; l++) {
				V[cpt_bin] = byte_bitset[l];
				cpt_bin++;
			}
			cpt_byte++;
		}
		return true;
	}
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
ifile::ifile() {
}

ifile::ifile(string filename , bool binary) {
	open(filename, binary);
}

ifile::~ifile() {
	close();
}

string ifile::name() {
	return file;
}

bool ifile::open(string filename, bool binary) {
	file = filename;
	string ext = futils::extensionFile(filename);
	if (ext == "gz") {
		fd.open(file.c_str(), ios::in | ios::binary);
		push(bio::gzip_decompressor());
	} else if (ext == "bz2") {
		fd.open(file.c_str(), ios::in | ios::binary);
		push(bio::bzip2_decompressor());
	} else if (binary) {
		fd.open(file.c_str(), ios::in | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

bool ifile::readString(string & str) {
	int s;
	if (!read((char*)&s, sizeof(int))) return false;
	char  * buffer = new char [s + 1];
	if (!read(buffer, s)) return false;
	buffer[s] = '\0';
	str = string(buffer);
	delete [] buffer;
	return true;
}

void ifile::close() {
	 if (!empty()) reset();
	 fd.close();
}

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
ofile::ofile() {
}

ofile::ofile(string filename , bool binary) {
	open(filename, binary);
}

ofile::~ofile() {
	close();
}

string ofile::name() {
	return file;
}

bool ofile::open(string filename, bool binary) {
	file = filename;
	string ext = futils::extensionFile(filename);
	if (ext == "gz") {
		fd.open(file.c_str(), ios::out | ios::binary);
		push(bio::gzip_compressor());
	} else if (ext == "bz2") {
		fd.open(file.c_str(), ios::out | ios::binary);
		push(bio::bzip2_compressor());
	} else if (binary) {
		fd.open(file.c_str(), ios::out | ios::binary);
	} else  {
		fd.open(file.c_str());
	}
	if (fd.fail()) return false;
	push(fd);
	return true;
}

void ofile::writeString(string & str) {
	int s = str.size();
	write((char*)&s, sizeof(int));
	write(str.c_str(), s * sizeof(char));
}

void ofile::close() {
	 if (!empty()) reset();
	 fd.close();
}

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
lfile::lfile() {
	verboseC = true;
	verboseL = true;
}

lfile::~lfile() {
	close();
}

string lfile::name() {
	return file;
}

bool lfile::open(string filename) {
	file = filename;
	if (futils::extensionFile(file) != "log") file += ".log";
	fd.open(file.c_str());
	if (fd.fail()) return false;
	return true;
}

void lfile::close() {
	 fd.close();
}

string lfile::getPrefix() {
	return file.substr(0, file.find_last_of("."));
}

void lfile::muteL() {
	verboseL = false;
}

void lfile::unmuteL() {
	verboseL = true;
}

void lfile::muteC() {
	verboseC = false;
}

void lfile::unmuteC() {
	verboseC = true;
}


void lfile::print(string s) {
	if (verboseL) { fd << s; fd.flush(); }
	if (verboseC) { cout << s; cout.flush(); }
}

void lfile::printC(string s) {
	if (verboseC) { cout << s; cout.flush(); }
}

void lfile::printL(string s) {
	if (verboseL) { fd << s; fd.flush(); }
}

void lfile::println(string s) {
	if (verboseL) { fd << s << endl; }
	if (verboseC) { cout << s << endl; }
}

void lfile::printlnC(string s) {
	if (verboseC) { cout << s << endl; }
}

void lfile::printlnL(string s) {
	if (verboseL) { fd << s << endl; }
}

void lfile::warning(string s) {
	cout << endl << "\33[33mWARNING:\33[0m " << s << endl;
	if (verboseL) fd << endl << "WARNING: " << s << endl;
}

void lfile::error(string s) {
	cout << endl << "\33[33mERROR:\33[0m " << s << endl;
	if (verboseL) fd << endl << "ERROR: " << s << endl;
	close();
	exit(1);
}
