/*
 * Histogram.cpp
 *
 *  Created on: Sep 15, 2016
 *      Author: konbob
 */

#include "Histogram.h"

#include <boost/math/special_functions/zeta.hpp> // for PowerLaw ML
#include <boost/math/tools/minima.hpp> // for PowerLaw ML
#include <boost/math/tools/roots.hpp> // for PowerLaw ML

#include <armadillo> // for nonlinear least squares

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <exception>
#include <fstream>
#include <utility>
#include <array>

Histogram::Histogram(std::vector<double>& rawdata, double binWidth) {

	if (rawdata.size() > 0 && binWidth >0) {
		rawData = rawdata;
		// sort rawdata
		std::sort(rawData.begin(), rawData.end());
		// get minimum,maximum and number of bins
		minimum = rawData.front();
		maximum = rawData.back();
		this->binsize = binWidth;
		// calculate bin number
		numberOfBins = ceil((maximum - minimum) / binsize)+1;

		// make sure data member variable is correctly initialized
		binCounts.assign(numberOfBins, 0);

		double bin = 0;
		for (auto value : rawData) {
			// calculate bin and increment counter
			bin = floor((value - minimum) / binsize);
			try {
				binCounts.at(bin)++;
			} catch (...) {
				std::cerr << "Forbidden index for value " << value << std::endl;
				std::cerr << " min " << minimum << " max " << maximum << " binWidth " << binWidth << " binnumber " << numberOfBins << std::endl;
			}

		}

		// set up flag for sums
		sumsCalculated = false;

		// initialize sums
		sumOfEntries = 0.0;
		sumOfSquaredEntries  = 0.0;
		guess = 0.0;

	} else {
		std::cerr << "Empty data set or forbidden bin width (passed value: "<< binWidth<<") at histogram constructor. Exit."<< std::endl;
		exit(1);
	}

}

/*
 * Prints all bins of the histogram
 */
void Histogram::printHistogram(const std::string& message) const {
	std::cerr<< message << std::endl;
	for (unsigned int var = 0; var < binCounts.size(); ++var) {
		std::cerr<< minimum +var*binsize<<" "<< binCounts.at(var) << std::endl;
	}
}

/*
 * Writes all bins of the histogram to file with given name in tabular separated format
 * Arguments:
 * filename - file with data
 * title - title of the plot in pdf
 * xLabel - Label for x axis of the plot
 * yLabel - Label for y axis of the plot
 * plotname - name of pdf
 *
 * For use with script Histpoplot.py
 */
void Histogram::writeHistogram(std::string filename,std::string title,std::string xLabel,std::string yLabel, std::string plotname) const  {

	std::ofstream file;
	try {
		// Open file
		file.open(filename,std::ofstream::out);

		// write legend to file
		file << title << std::endl;
		file << xLabel << std::endl;
		file << yLabel << std::endl;
		file << plotname << std::endl;

		//write data to file
		for (unsigned int var = 0; var < binCounts.size(); ++var) {
			file << minimum +var*binsize <<"  "<< binCounts.at(var) << std::endl;
		}

		// close file
		file.close();
		std::cerr << "#Finished writing data to " << filename << std::endl;

	} catch (...) {
		std::cerr << "Could not open or write to file " << filename << std::endl;
	}
}

/*
 * Returns the value x of the p percentile  ( sum_{i}^{x} = p * sum_{i})
 * NOTE: x will be linearly interpolated between integer values
 */
double Histogram::getPercentile(const double p) const  {
	if (0 <= p && p <= 1) {

		// cumulate the entries
		std::vector<double> cumulatedCounts(binCounts.size(), 0);
		cumulatedCounts.at(0) = binCounts.at(0);
		for (unsigned int bin = 1; bin < binCounts.size(); ++bin) {
			try {
				cumulatedCounts.at(bin) = binCounts.at(bin) + cumulatedCounts.at(bin - 1);
			} catch (...) {
				std::cout << "Out of Range error at getPercentile() " << std::endl;
			}
		}
		// calculate the p value
		double percentileValue = p * cumulatedCounts.back();

		// get the iterator of the first value greater than the percentile value
		auto upper = std::upper_bound(cumulatedCounts.begin(), cumulatedCounts.end(), percentileValue);

		// get the index of the iterator
		int index = std::distance(cumulatedCounts.begin(), upper);

		double upperValue, lowerValue, xi;
		try {
			if (index > 1) {
				// calculate the non integer part of the interpolation
				upperValue = cumulatedCounts.at(index);
				lowerValue = cumulatedCounts.at(index - 1);

			} else {
				// in that case most of the data is in the first bin
				return minimum;
			}
			xi = (percentileValue - lowerValue) / (upperValue - lowerValue);

		} catch (...) {
			std::cerr << "Out of Range error at getPercentile() - lower / upper part. " << std::endl;
			std::cerr << "upperValue: "<< upperValue << " lower Value: "<<lowerValue << " index "<< index<<std::endl;
			exit(1);
		}
		// calculate the interpolated bin value
		double binValue = index - 1 + xi;

		// return the percentile
		return minimum+ binsize*binValue;

	} else {
		std::cout << "p must be in range [0,1] " << std::endl;
		exit(1);
	}
}

/*
 * Calculates the sum of entries and the sum of squared entries according to Blobel and Lohrmann with high precision.
 */
void Histogram::calculateSums() {

	// use guess for better precision
	guess = minimum + binsize * (floor(binCounts.size() / 2.0));

	// loop over entries
	for (unsigned long int index = 0; index < binCounts.size(); ++index) {
		sumOfEntries += binCounts.at(index) * ((minimum + binsize * index)-guess);
		sumOfSquaredEntries += binCounts.at(index) * pow(minimum + binsize * index-guess, 2.0);
	}
}

/*
 * Returns an unbiased estimator for the mean of the distribution.
 * Numerical scheme according to Blobel & Lohrmann: "Statistische und numerische Methoden der Datenanalyse", page 94
 */
double Histogram::getMeanEstimator(){
	// check if sums have already been calculated
	if (sumsCalculated==false) {
		calculateSums();
		sumsCalculated = true;
	}

	// calculate estimate and return it
	return guess +(sumOfEntries)/(rawData.size());


}

/*
 * Returns the error of the mean estimator.
 * Numerical scheme according to Blobel & Lohrmann: "Statistische und numerische Methoden der Datenanalyse", page 94
 */
double Histogram::getErrorOfMeanEstimator() {
	// check if sums have already been calculated
	if (sumsCalculated == false) {
		calculateSums();
		sumsCalculated = true;
	}

	// calculate variance
	double n = rawData.size();
	if (n == 1.0) {
		std::cerr << "WARNING: Sample size is one, infinity will be produced. " << std::endl;
	}
	double variance = (sumOfSquaredEntries - (pow(sumOfEntries, 2.0)) / (n)) / (n - 1.0);

	// calculate error and return it
	return sqrt((variance) / (n));
}
/*
 * Returns an unbiased estimator for the variance of the distribution.
 * Numerical scheme according to Blobel & Lohrmann: "Statistische und numerische Methoden der Datenanalyse", page 94
 */
double Histogram::getVarianceEstimator() {
	// check if sums have already been calculated
	if (sumsCalculated == false) {
		calculateSums();
		sumsCalculated = true;
	}

	// calculate variance
	double n = rawData.size();
	if (n == 1.0) {
		std::cerr << "WARNING: Sample size is one, infinity will be produced. " << std::endl;
	}
	return (sumOfSquaredEntries - (pow(sumOfEntries, 2.0)) / (n)) / (n - 1.0);
}
/*
 * Returns an error of the variance estimator under the ASSUMPTION THAT THE SAMPLE IS NORMAL DISTRIBUTED.
 * Numerical scheme according to Blobel & Lohrmann: "Statistische und numerische Methoden der Datenanalyse", page 94
 */
double Histogram::getVarianceEstimatorErrorNormalDistributed() {
	// check if sums have already been calculated
	if (sumsCalculated == false) {
		calculateSums();
		sumsCalculated = true;
	}

	// calculate variance
	double n = rawData.size();
	if (n == 1.0) {
		std::cerr << "WARNING: Sample size is one, infinity will be produced. " << std::endl;
	}
	double variance = (sumOfSquaredEntries - (pow(sumOfEntries, 2.0)) / (n)) / (n - 1.0);

	// calculate error and return it
	return sqrt((variance) / (2.0*(n-1.0)));
}

void Histogram::writeRawData(std::string filename) const  {
	printf("writing raw data!\n");
	std::ofstream file;
		try {
			// Open file
			file.open(filename,std::ofstream::out);

			//write data to file
			for (unsigned int var = 0; var < rawData.size(); ++var) {
				file << rawData.at(var) << std::endl;
			}

			// close file
			file.close();
			std::cerr << "Finished writing data to " << filename << std::endl;

		} catch (...) {
			std::cerr << "Could not open or write to file " << filename << std::endl;
		}


}

/*
 * Does numerical maximum likelihood fit to power law p(x) = (a-1)/(xmin)*((x)/(xmin))**(-a) according to
 *  SIAM REVIEW Vol. 51, No. 4, pp. 661–703 and error estimation according to Blobel & Lohrmann p.131.
 *
 *  Return values in order:
 *  - estimate of exponent
 *  - left hand side error
 *  - right hand side error
 *
 */
std::vector<double> Histogram::doPowerLawMLFit() {

	// for normalized data
	std::vector<double> normalizedData(rawData.size());

	// check if sum of entries is already calculated
	if (!sumsCalculated) {
		calculateSums();
	}

	// normalize data
	for (unsigned long int index = 0; index < normalizedData.size(); ++index) {
		normalizedData.at(index) = rawData.at(index)/sumOfEntries;
	}

	// get sum of logarithms
	double sumLogs = 0;
	for (double& date : rawData) {
		if (date != 0.0) {
			sumLogs += log(date);
		}
	}

	// get the size of the data
	double n = rawData.size();

	// ================== Estimator ================================
	// calculate the best exponent by minimizing the negative loglikelihood numerically

	// struct to for the negative loglikelihood function
	struct NegLogLikelihood {

		// short notation of constructor, initializes  the parameters n and  sum of logs
		NegLogLikelihood(double n, double sum): n(n), sumLogs(sum) {}

		// return the negative loglikelihood for a given alpha
		double operator()(double const& alpha) {
			return n*boost::math::zeta(alpha)+alpha*sumLogs;
		}

	private:
		double n;
		double sumLogs;

	};

	// create instance of function struct
	NegLogLikelihood negLogLike = NegLogLikelihood(n,sumLogs);

	// print out for visual check
	bool printOut = false;
	if (printOut == true) {
		for (double alpha = 1.1; alpha< 4; alpha+=0.01) {
			printf("%f %f\n",alpha, negLogLike(alpha));
		}
	}

	//  do numerical search for minimum in range [a_min,a_max] for estimate
	double a_min = 1.1;
	double a_max = 4.0;
	int bits = std::numeric_limits<double>::digits;
	std::pair<double, double> estimatePair =boost::math::tools::brent_find_minima(negLogLike, a_min, a_max, bits);

	// results
	double exponent = estimatePair.first; // estimated value of exponent
	double minimumValue = estimatePair.second; // mimimum value of NegLogLikelihood

	// ======================== Error of estimator =============================
	// calculate the error by find the roots of negLogLike = min(negLogLike) + 0.5

	// struct to for the shifted negative loglikelihood function
	struct ShiftedNegLogLikelihood {

		// short notation of constructor, initializes  the parameters n, sum of logs and minimum
		ShiftedNegLogLikelihood(double n, double sum, double min) :	n(n), sumLogs(sum) , min(min){}

		// return the shifted negative loglikelihood for a given alpha
		double operator()(double const& alpha) {
			return n * boost::math::zeta(alpha) + alpha * sumLogs -min -0.5 ;
		}

	private:
		double n;
		double sumLogs;
		double min;

	};

	// create instance
	ShiftedNegLogLikelihood shiftNegLogLike = ShiftedNegLogLikelihood(n,sumLogs,minimumValue);

	// print out for visual check
	printOut = false;
	if (printOut == true) {
		for (double alpha = 2.09; alpha < 2.2; alpha += 0.0001) {
			printf("%f %f\n", alpha, shiftNegLogLike(alpha));
		}
	}

	// create tolerance function
	boost::math::tools::eps_tolerance<double> tolerance(bits);

	// search for left sided error
	double leftBound = 1.0001*exponent; // must be a value where the sign has already changed but not to close to second root
	std::pair<double,double> leftRootPair = boost::math::tools::bisect(shiftNegLogLike,a_min,leftBound,tolerance);
	double errorLeft  = exponent -0.5*(leftRootPair.first+ leftRootPair.second);

	// search for left sided error
	double rightBound = 0.9999 * exponent; // must be a value where the sign has already changed but not to close to second root
	std::pair<double, double> rightRootPair = boost::math::tools::bisect(shiftNegLogLike, rightBound,a_max ,tolerance);
	double errorRight = 0.5 * (rightRootPair.first + rightRootPair.second)-exponent;

	//return the results
	std::vector<double> result  = {exponent,errorLeft,errorRight};
	return result;

}

/*
 * Does a nonlinear least squares fit to a power law f(x) = a0 * (x**a1).
 * Method according to Blobel and Lohrmann p.166
 *
 * structure of return value
 * {a0,a1,StdError[a0],StdError[a1],Cov(a0,a1),Cov(a1,a0)}
 */
std::vector<double> Histogram::doNonlinearPowerFit(const std::vector<double>& xValues, const std::vector<double>& meanValues, const std::vector<double>& errors){

	// store size of data
	unsigned long int dataSize = xValues.size();

	// ======= check input ========
	if (meanValues.size() != dataSize || errors.size() != dataSize) {
		std::cerr << " Error at doNonlinearPowerFit: Input of different lengths."<< std::endl;
		exit(1);
	}
	for (unsigned int index = 0; index < dataSize; ++index) {
		if (errors.at(index) <= 0.0) {
			std::cerr << " Error at doNonlinearPowerFit: Forbidden error value (" << errors.at(index) << ") at index " << index << "." << std::endl;
			std::cerr << "Output of data: \n x 		y		 error" << std::endl;
			for (unsigned int index = 0; index < dataSize; ++index) {
				std::cerr << xValues.at(index) <<" "<< meanValues.at(index) <<" "<< errors.at(index)<<std::endl;
			}
			exit(1);
		}
	}

	// store number of fit parameters
	int numberFitParam = 2;

	// ======= set up variables =====

	// for estimating starting values (Blobel and Lohrmann eq (7.65))
	std::vector<double> logX(dataSize,1.0);
	std::vector<double> logY(dataSize,1.0);
	double sumX = 0.;
	double sumY = 0.;
	double sumXX = 0.;
	double sumXY = 0.;

	for (unsigned long int count = 0; count < dataSize; ++count) {
		// get the logarithms
		logX.at(count) = log(xValues.at(count));
		logY.at(count) = log(meanValues.at(count));

		// calculate the sums
		sumX += logX.at(count);
		sumY += logY.at(count);
		sumXX += pow(logX.at(count),2.0);
		sumX += logX.at(count)*logY.at(count);

	}

	// estimated starting values
	double D = dataSize*sumXX-pow(sumX,2.0);
	double aZeroGuess =exp((sumXX*sumY-sumX*sumXY)/(D));
	double aOneGuess = exp((dataSize*sumXY-sumX*sumY)/D);

	// parameter vectors
	arma::vec aGuess ={aZeroGuess,aOneGuess};
	arma::vec deltaA = {0.,0.};

	// weights matrix
	arma::vec weightVector(dataSize);
	for (unsigned int index = 0; index < dataSize; ++index) {
		weightVector(index) = 1.0/pow(errors.at(index),2.0);
	}
	arma::mat W = arma::diagmat(weightVector);

	// data vectors
	arma::vec y(meanValues);
	arma::vec x(xValues);

	// hessian matrix
	arma::mat A(dataSize,numberFitParam);

	// covariance matrix
	arma::mat V(numberFitParam,numberFitParam);

	// fit values
	arma::vec f(dataSize);

	// ===== do calculation =====

	for (int count = 0; count < 18; ++count) {
		// set up fit values and hessian matrix
		for (unsigned int index = 0; index < dataSize; ++index) {
			f(index) = aGuess(0) * pow(x(index), aGuess(1));
			A(index, 0) = pow(x(index), aGuess(1));
			A(index, 1) = f(index)* log(x(index));

		}
		// calculate the variance matrix
		V = arma::inv(A.t() * W * A);

		// calculate the increment for the parameters
		deltaA = V * A.t() * W * (y - f);

		// calculate new parameters
		aGuess += deltaA;

		// check if increment is small enough to finish calculation
		if (arma::norm(deltaA, 2) / arma::norm(aGuess, 2) < 1e-2 ) {
			break;
		}
	}

	//return the results
	std::vector<double> result(3*numberFitParam);
	result.at(0) = aGuess(0);
	result.at(1) = aGuess(1);
	result.at(2) = sqrt(V(0,0));
	result.at(3) = sqrt(V(1,1));
	result.at(4) = V(0,1);
	result.at(5) = V(1,0);
	return result;
}

/*
 * Loads file for input in Histogram::doNonlinearPowerFit.
 * Expected file format:
 * || 4 comment lines (symbol '#')
 * || x meanY error fit (white space separated values)
 * Fit values will be ignored.
 *
 *Return values in order
 * x, meanY, error
 */
std::array<std::vector<double>,3> Histogram::loadTxtForFit(std::string filename) {
	try {

		// file object for reading
		std::ifstream file;

		// open the file
		file.open(filename, std::ifstream::in);

		// skip the first four lines
		std::string line;
		for (int count = 0; count < 4; ++count) {
			std::getline(file, line);
		}

		// ===== parse the lines one to get the number of nodes =====
		// prepare vectors for entries
		std::vector<double> xValues;
		std::vector<double> meanValues;
		std::vector<double> errors;

		while (std::getline(file, line)) {

			// split it into parts at delimeter " "
			std::stringstream stringstream;
			stringstream.str(line);
			std::string substring;
			std::vector<std::string> parsedEntries;
			while (std::getline(stringstream, substring, ' ')) {
				parsedEntries.push_back(substring);
			}
			// cast entries and push them back
			xValues.push_back(std::stod(parsedEntries.at(0)));
			meanValues.push_back(std::stod(parsedEntries.at(1)));
			errors.push_back(std::stod(parsedEntries.at(2)));
		}

		// close file
		file.close();

		// return values
		std::array<std::vector<double>,3> Data;
		Data[0] = xValues;
		Data[1] = meanValues;
		Data[2] = errors;

		return Data;

	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		exit(1);
	}
}





































