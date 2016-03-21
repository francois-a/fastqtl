
#define FASTQTL_USE_QR
#define FASTQTL_QR_TOLERANCE 1e-7 //this is set to the same tolerance as R

#include "residualizer.h"


residualizer::residualizer(int _n_samples) : n_samples (_n_samples) {
	matrices_uptodate = false;
}

residualizer::~residualizer() {
	hard_covariates.clear();
	soft_covariates.clear();
	m_rank = 0;
	covarM.resize(0,0);
    PQR_Q.resize(0,0);
    PQR_Q_A.resize(0,0);
    matrices_uptodate = false;
}

int residualizer::nCovariates() {
	return soft_covariates.size() + hard_covariates.size();
}

void residualizer::clearSoft() {
	if (soft_covariates.size() > 0) {
		soft_covariates.clear();
		matrices_uptodate = false;
	}
}

void residualizer::pushHard(vector < string > & covariate) {
	//INIT
	set < int > i_yesmissing, i_nonmissing;
	set < string > factors;

	//TEST FOR NUMERIC & MISSING
	for (int i = 0 ; i < covariate.size() ; i++) {
		if (covariate[i] == "NA") i_yesmissing.insert(i);
		else {
			if (!sutils::isNumeric(covariate[i])) factors.insert(covariate[i]);
			i_nonmissing.insert(i);
		}
	}

	//FILL IN VALUES
	unsigned int starting_row = hard_covariates.size();
	if (factors.size() == 1) LOG.warning("Non-numeric covariate with only 1 factor, covariate dropped!");
	else if (factors.size() > 1) {
		set < string >::iterator itF = factors.begin(); itF ++;
		for (; itF != factors.end() ; ++itF) {
			hard_covariates.push_back(vector < float > (n_samples, 0.0));
			for (int i = 0 ; i < covariate.size() ; i++) if (covariate[i] == *itF) hard_covariates.back()[i] = 1.0;
		}
	} else {
		hard_covariates.push_back(vector < float > (n_samples, 0.0));
		for (int i = 0 ; i < covariate.size() ; i++) hard_covariates.back()[i] = atof(covariate[i].c_str());
	}

	//IMPUTE MISSING
	if (i_yesmissing.size() > 0) {
		for (int i_row = starting_row ; i_row < hard_covariates.size() ; i_row ++) {
			float mean_row = 0;
			for (set < int >::iterator itNM = i_nonmissing.begin(); itNM != i_nonmissing.end() ; ++itNM) mean_row += hard_covariates[i_row][*itNM];
			mean_row /= i_nonmissing.size();
			for (set < int >::iterator itYM = i_yesmissing.begin(); itYM != i_yesmissing.end() ; ++itYM) hard_covariates[i_row][*itYM] = mean_row;
		}
	}

	matrices_uptodate = false;
}

void residualizer::pushSoft(vector < float > & covariate) {
	for (int i = 1, same = 1 ; i < covariate.size(); i++) {
		if (covariate[i] != covariate[i-1]) same = 0;
		if (i == covariate.size() - 1 && same == 1) return;
	}
	soft_covariates.push_back(vector < float > (n_samples, 0.0));
	for (int i = 0 ; i < covariate.size() ; i++) soft_covariates.back()[i] = covariate[i];
	matrices_uptodate = false;
}

void residualizer::pushHard(vector < float > & covariate) {
	for (int i = 1, same = 1 ; i < covariate.size(); i++) {
		if (covariate[i] != covariate[i-1]) same = 0;
		if (i == covariate.size() - 1 && same == 1) return;
	}
	hard_covariates.push_back(vector < float > (n_samples, 0.0));
	for (int i = 0 ; i < covariate.size() ; i++) hard_covariates.back()[i] = covariate[i];
	matrices_uptodate = false;
}

void residualizer::update() {
	covarM.resize(n_samples,nCovariates() + 1);
	covarM.col(0) = VectorXd::Ones(n_samples);
	for(int h = 0; h < hard_covariates.size() ; h ++) for(int c = 0 ; c < n_samples ; c ++) covarM(c, h + 1) = hard_covariates[h][c];
	for(int s = 0; s < soft_covariates.size() ; s ++) for(int c = 0 ; c < n_samples ; c ++) covarM(c, s + hard_covariates.size() + 1) = soft_covariates[s][c];
	PQR = ColPivHouseholderQR<MatrixXd>(covarM);
	PQR.setThreshold(FASTQTL_QR_TOLERANCE);
	m_rank = PQR.rank();
	if (m_rank != nCovariates()+1) {
		PQR_Q = PQR.householderQ();
		PQR_Q_A = PQR.householderQ().adjoint();
		LOG.warning("Linear dependency between covariates [#dropped=" +sutils::int2str(nCovariates()+1-m_rank) + "]");
	}
	matrices_uptodate = true;
}

void residualizer::residualize(vector < float > & data) {
	if (nCovariates() == 0) return;
	//TEST IF DATA VARIABLE
	for (int i = 1, same = 1 ; i < data.size(); i++) {
		if (data[i] != data[i-1]) same = 0;
		if (i == data.size() - 1 && same == 1) return;
	}

	//FILL IN DATA
	VectorXd counts(n_samples);
	for(int i = 0; i < n_samples ; i ++) counts(i) = (double) data[i];

	//MATRICES UPDATE
	if (!matrices_uptodate) update();

	//CORRECTION
	if (m_rank == nCovariates()+1) {
		VectorXd m_coef = PQR.solve(counts);
		VectorXd fitted = covarM * m_coef;
		VectorXd e = counts - fitted;
		for (int i = 0; i < e.size(); i ++) data[i] = e(i);
	} else {
		VectorXd effects(PQR_Q_A * counts);
		effects.tail(n_samples - m_rank).setZero();
		VectorXd fitted = PQR_Q * effects;
		VectorXd e = counts - fitted;
		for (int i = 0; i < e.size(); i ++) data[i] = (float)e(i);
	}
}

void residualizer::residualize(vector < vector < float > > & data) {
	for (int d = 0; d < data.size() ; d ++) residualize(data[d]);
}
