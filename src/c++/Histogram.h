/*
 * Histogram.h
 *
 *  Created on: Sep 15, 2016
 *      Author: konbob
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_
#include <vector>
#include <string>
#include <array>


class Histogram {
public:
	Histogram(std::vector<double>& rawdata,double binWidth);

	// output methods
	void printHistogram(const std::string& message) const ;
	void writeHistogram(std::string filename,std::string title,std::string xLabel,std::string yLabel, std::string plotname) const ;
	void writeRawData(std::string filename) const ;

	// input methods
	std::array<std::vector<double>,3> loadTxtForFit(std::string filename);

	// analysis methods
	double getPercentile(const double value) const ;
	double getMeanEstimator();
	double getErrorOfMeanEstimator();
	double getVarianceEstimator();
	double getVarianceEstimatorErrorNormalDistributed();

	// fitting methods
	std::vector<double> doPowerLawMLFit();
	static std::vector<double> doNonlinearPowerFit(const std::vector<double>& xValues,const std::vector<double>& meanValues, const std::vector<double>& errors);

	// watch out if methods change the raw data! Then it could be necessary to recalculate the sums!
private:
	//methods
	void calculateSums();

	// member
	std::vector<double> binCounts;
	std::vector<double> rawData;
	double minimum ;
	double maximum;
	int numberOfBins;
	double binsize;
	bool sumsCalculated;
	double sumOfEntries;
	double sumOfSquaredEntries;
	double guess; // for higher numerical precision



};

#endif /* HISTOGRAM_H_ */
