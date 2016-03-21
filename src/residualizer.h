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

#ifndef _RESIDUALIZER_H
#define _RESIDUALIZER_H

#include <Dense>
#include <LU>
#include "utils/utils.h"

using namespace Eigen;

class residualizer {
public:

	int n_samples;
	vector < vector < float > > hard_covariates;
	vector < vector < float > > soft_covariates;

	bool matrices_uptodate;
    int m_rank;
    MatrixXd covarM;
    MatrixXd PQR_Q;
    MatrixXd PQR_Q_A;
    ColPivHouseholderQR < MatrixXd > PQR;

    residualizer(int);
    ~residualizer();

    int nCovariates();
    void clearSoft();
    void pushHard(vector < string > &);
    void pushHard(vector < float > &);
    void pushSoft(vector < float > &);
    void update();
    void residualize(vector < float > &);
    void residualize(vector < vector < float > > &);
};

#endif
