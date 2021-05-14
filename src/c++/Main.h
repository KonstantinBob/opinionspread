/*
 * Main.h
 *
 *  Created on: Sep 16, 2016
 *      Author: konbob
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "Graph.h"
#include "LearningNetwork.hpp"

class Main {
public:
	template<typename T> static void effDiaProduction(T& graph);
	template<typename T> static void DPL_Production(T& graph, const bool& doFit);
	template<typename T> static void degreeDistribution_Production(T& graph, const bool& doFit);
	template<typename T> static void pairSimulation_Production(const T& graph, const double& p);
	template<typename T> static LearningNetwork<T> generateTimeSeries(T& graph,const double& p, const double& q);
	template<typename T> static std::vector<double> generateIsingTimeSeries(T& graph,const double& J, const double& h,const double& beta);
	template<typename T> static void generateCoarseGrainedView(T& graph,const double& p, const double& q);
	template<typename T> static void writeDegreeDifference(T& graph);
	template<typename T> static void writeCorrelationGraph(T& graph,const double& p, const double& q);
	template<typename T> static double computeAverageSkew(T& graph);
};



#endif /* MAIN_H_ */
