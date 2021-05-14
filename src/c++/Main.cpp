/*
 * Main.cpp
 *
 *  Created on: Sep 16, 2016
 *      Author: konbob
 */

#include "Main.h"
#include "Graph.h"
#include "Histogram.h"
#include "Scientist.h"
#include "RandomNumbers.h"
#include "ForestFire.h"
#include "Erdos.h"
#include "StochKronecker.h"
#include "Arxiv.h"
#include "BarabasiAlbert.h"
#include "WattsStrogatz.h"

#include "Clustering.hpp"
#include "LearningNetwork.hpp"

#include <armadillo>

#include <string>
#include <set>
#include <ctime>
#include <chrono>
#include <fstream>
#include <stdio.h>

int main(int argc, char *argv[]) {

	// for command line arguments
	double q, p;
	int nodes;
	int netType;

	// read the command line args
	if (argc == 5) {
		p = std::atof(argv[1]);
		q = std::atof(argv[2]);
		nodes = std::atof(argv[3]);
		netType = std::atoi(argv[4]);
	} else {
		printf("No arguments !\n");
		exit(1);
	}

	// ==== creation of networks
	int littlenodes = 111;

	//Arxiv a("arxivData.snap");

	double forward = 0.37;
	double backward = 0.32;
	ForestFire ff(forward, backward, littlenodes);

	Erdos er((1.0 * ff.getNumberOfEdges()) / (littlenodes * (littlenodes - 1)), littlenodes);

	//==== run simulations ====
	if (netType == 0) {
		printf("#ff with size %d p %f q %f \n",nodes,p,q);
		ff.resize(nodes);
		//Main::generateIsingTimeSeries(ff,1.0,p,q);
		Main::generateCoarseGrainedView(ff,p,q);
		//Main::generateTimeSeries(ff,p,q);
		//printf("ff %f\n",Main::computeAverageSkew(ff));

	} else if (netType == 1) {
		printf("#er with size %d p %f q %f \n",nodes,p,q);
		er.resize(nodes);
		//Main::generateIsingTimeSeries(er,1.0,p,q);
		Main::generateCoarseGrainedView(er,p,q);
		//Main::generateTimeSeries(er,p,q);
		//printf("er %f\n",Main::computeAverageSkew(er));

	} else if (netType == 2) {
		printf("#sk with size %d p %f q %f \n",nodes,p,q);
		arma::mat iterMat( { { 0.999, 0.437 }, { 0.437, 0.484 } });
		StochKronecker Sto(iterMat, nodes);
		Main::generateCoarseGrainedView(Sto,p,q);
		//Main::generateIsingTimeSeries(Sto,1.0,p,q);
		//Main::generateTimeSeries(Sto,p,q);
		//printf("Sto %f\n",Main::computeAverageSkew(Sto));

	} else if (netType == 3) {
		printf("#ws with size %d p %f q %f \n",nodes,p,q);
		WattsStrogatz w(13, 0.1, nodes);
		Main::generateCoarseGrainedView(w,p,q);
		//Main::generateIsingTimeSeries(w,1.0,p,q);
		//Main::generateCoarseGrainedView(w,p,q);
		//Main::generateTimeSeries(w, p, q);
		//printf("w %f\n",Main::computeAverageSkew(w));

	} else if (netType == 4) {
		printf("#ba with size %d p %f q %f \n",nodes,p,q);
		BarabasiAlbert ba(13, 13, nodes);
		Main::generateCoarseGrainedView(ba,p,q);
		//Main::generateIsingTimeSeries(ba,1.0,p,q);
		//Main::generateTimeSeries(ba, p,q);
		//Main::generateCoarseGrainedView(ba,p,q);
		//printf("ba %f\n",Main::computeAverageSkew(ba));


	}
}

template<typename T> inline double Main::computeAverageSkew(T& graph) {
	// prepare sizes list
	std::vector<long int> sizes = {100, 200,400,800,1600,3200};

	// prepare number of trails
	int trailMax = 20;

	// prepare accumulator
	double skew = 0.0;

	// loop over sizes
	for (const long int& size :sizes ){
		for (int trial = 0; trial < trailMax; ++trial) {
			skew+= graph.getSkewness();
		}
	}

	return skew/(sizes.size() * trailMax);
}

template<typename T> inline void Main::writeCorrelationGraph(T& graph, const double& p, const double& q) {

	// create Learning network
	LearningNetwork<T> ln(graph, p, q);

	// create meta graph
	MetaGraph meta = ln.createCorrelationGraph(5000);

	// write results
	meta.writeGephiFile("meta" + graph.getTypeIdentifier());
}

template<typename T> void Main::pairSimulation_Production(const T& graph, const double& p) {

	// list of qValues
	std::vector<double> qValues = { 0.0, 0.1, 0.2, 0.3, 0.4, .5, .6, .7, .8, .9, 1.0 };

	// loop over qValues
	for (double q : qValues) {

		// prepare the network
		LearningNetwork<T> ln(graph, p, q);

		// store distribution
		std::array<double, 4> opDistr = { 0., 0., 0., 0. };

		// number of samples
		double sampleMax = 1e5;

		for (double sample = 0; sample < sampleMax; ++sample) {
			// perform cycles
			for (int counter = 1; counter < 50; ++counter) {
				ln.performCycle();
			}

			// output
			std::array<double, 4> opinionVec = ln.calculatePairDistributions();

			// accumulate
			for (int index = 0; index < opinionVec.size(); ++index) {
				opDistr.at(index) += opinionVec.at(index);
			}
		}

		// final output
		for (int index = 0; index < opDistr.size(); ++index) {
			double ratio = opDistr.at(index) / sampleMax;
			printf("%f %f \n", ratio, sqrt((ratio * (1.0 - ratio)) / (sampleMax)));
		}
	}

}

template<typename T> LearningNetwork<T> Main::generateTimeSeries(T& graph, const double& p, const double& q) {

	// generate learning Network
	LearningNetwork<T> ln(graph, p, q);

	//  Perform the cycles
	for (int counter = 0; counter < 5e2; ++counter) {

		ln.performCycle();
		printf(" %f \n", ln.getFractionUp());

	}

	return ln;
}

/*
 * Write the difference between Out- and Indegree.
 */
template<typename T> inline void Main::writeDegreeDifference(T& graph) {

	// file object
	std::ofstream file;

	// create filename
	std::string filename = "DegreeDiff_" + graph.getTypeIdentifier();

	try {
		// Open file
		file.open(filename, std::ofstream::out);

		// write heade
		file << "# Node OutDegree-Indegree" << std::endl;

		// loop over the nodes
		for (unsigned long int node = 0; node < graph.getSize(); ++node) {
			file << node << " " << (double) graph.getAdjacentOutNodes(node).size() - (double) graph.getAdjacentInNodes(node).size() << std::endl;
		}

		// close file
		file.close();

	} catch (...) {
		std::cerr << "Could not open or write to file " << filename << "." << std::endl;
	}

}

template<typename T> void Main::effDiaProduction(T& graph) {

	std::printf("# size mean mean error var varError \n");

	for (double size = 0.02; size < 3.5; size += 0.5) {
		clock_t begin = clock();
		long int numberNodes = size * pow(10., 4.);

		double trials = 0.1 * numberNodes;

		if (trials < 100) {
			trials = numberNodes;
		}

		if (trials > 2000) {
			trials = 2000;
		}

		int countMax = 25;
		std::vector<double> diameters;
		diameters.reserve(countMax);

		for (int count = 0; count < countMax; ++count) {
			// resize the graph
			graph.resize(numberNodes);

			// measure effective diameter
			diameters.push_back(graph.getApproximateEffectiveDiameter(trials, 50).at(0));
		}

		Histogram histo(diameters, 0.01);
		std::printf("%ld %e %e %e %e \n", numberNodes, histo.getMeanEstimator(), histo.getErrorOfMeanEstimator(), histo.getVarianceEstimator(),
				histo.getVarianceEstimatorErrorNormalDistributed());
		clock_t end = clock();
		std::printf("#time for %d created networks of size %ld and bootstraped diameter estimation is %f hours. \n", countMax, numberNodes,
				double(end - begin) / (CLOCKS_PER_SEC * 3600.0));

	}
}

template<typename T> void Main::DPL_Production(T& graph, const bool& doFit) {

	// variables for experiment sizes
	int sizeMax = 15;
	int trailMax = 25;

	// vectors for storage
	std::vector<double> nodes(sizeMax, 0.);
	std::vector<double> edges(sizeMax, 0.);
	std::vector<double> error(sizeMax, 0.);

	// create different sizes
	for (int size = 0; size < sizeMax; ++size) {
		std::vector<double> trails(trailMax, 0);

		double graphSize = pow(1.5, size + 1) * 1e2;
		nodes.at(size) = graphSize;

		// sample several graphs
		for (int trail = 0; trail < trailMax; ++trail) {
			graph.resize(graphSize);
			trails.at(trail) = graph.getNumberOfEdges();
		}

		// get mean and error
		Histogram trailHisto(trails, 1.0);
		edges.at(size) = trailHisto.getMeanEstimator();
		error.at(size) = sqrt(trailHisto.getVarianceEstimator());

	}


	if (doFit == true) {
		// do fit
		std::vector<double> result = Histogram::doNonlinearPowerFit(nodes, edges, error);

		// print results with fit
		std::printf("# a0 = %f \\pm %f \n", result.at(0), result.at(2));
		std::printf("# a1 = %f \\pm %f \n", result.at(1), result.at(3));
		std::printf("# covariances %f %f \n", result.at(4), result.at(5));
		std::printf("# nodes edges error fit \n");
		for (unsigned int index = 0; index < nodes.size(); ++index) {
			std::printf("%f %f %f %f \n", nodes.at(index), edges.at(index), error.at(index), result.at(0) * pow(nodes.at(index), result.at(1)));
		}
	} else {
		// print results without fit
		for (unsigned int index = 0; index < nodes.size(); ++index) {
			std::printf("%f %f %f\n", nodes.at(index), edges.at(index), error.at(index));
		}
	}
}

template<typename T> void Main::degreeDistribution_Production(T& graph, const bool& doFit) {
	// create graph of right size
	unsigned long int desiredSize = 27770;
	if (graph.getNumberOfNodes() != desiredSize) {
		graph.resize(desiredSize);
	}

	// get the degree distribution
	Histogram histo = graph.getDegreeDistribution();

	if (doFit == true) {
		std::vector<double> fit = histo.doPowerLawMLFit();
		printf("#exponent %f error left %f error right %f \n", fit.at(0), fit.at(1), fit.at(2));
	}

	// write the data
	histo.writeHistogram("degreeDistr_" + graph.getTypeIdentifier(), "", "Degree", "Count", "degreeDistrPlot_" + graph.getTypeIdentifier());

}

template<typename T> std::vector<double> Main::generateIsingTimeSeries(T& graph, const double& J, const double& h, const double& beta) {
	printf("#Ising time series \n");

	// generate learning Network
	LearningNetwork<T> ln(graph, 0., 0.);
	printf("#ln size %ld  \n",ln.getSize());

	// prepare list of logarithmic ratios, i.e. ratios[i]= ln(P_{i+1}/P_0)
	std::vector<double> ratios(ln.getSize(),0.0);

	// prepare list of probabilities, i.e. prob[i] = P_{i}
	std::vector<double> prob(ln.getSize()+1,0.0);

	// sum for ratios
	double sumRatios = 0.0;

	//  do the cycles
	for (long int index = 0; index < ln.getSize(); ++index) {

		// get ratios from simulation
		ratios.at(index) = ln.doIsingCycleWindow(J, h, beta,index,1e3);

		// calculate cumulative sum, i.e. log of ratio
		prob.at(index+1) = ratios.at(index) + prob.at(index);
	}

	// === calculate probabilities ===

	// sum of p/p0
	double sumPdiv = 0.0;

	// remove logarithms and accumulate the sum pDiv
	for (long int index = 1; index <= ln.getSize(); ++index){
		prob.at(index) = exp(prob.at(index));
		sumPdiv += prob.at(index);
	}

	// calculate p_{0} by normalization
	double pZero = 1./(1. + sumPdiv);

	// fill pZero
	prob.at(0) = pZero;

	// multiply entries by pZero
	for (long int index = 1; index <= ln.getSize(); ++index){
		prob.at(index) = pZero * prob.at(index);

	}

	// output
	for (long int index = 0; index <= ln.getSize(); ++index){
		printf("%ld,%f\n",index,prob.at(index));
	}

	// return prob array
	return prob;


}

template<typename T> void Main::generateCoarseGrainedView(T& graph, const double& p, const double& q) {
	// generate learning Network
	LearningNetwork<T> ln(graph, p, q);

	Clustering<T> cl(ln);
	MetaGraph meta = cl.createMetaCluster();
	meta.writeGephiFile("ModularityOnly" + graph.getTypeIdentifier());

	/*graph.writeDotFile("Detail"+ graph.getTypeIdentifier()+".dot");

	//  Perform the cycles
	for (int counter = 0; counter < 5e2; ++counter) {

		ln.performCycle();
		printf(" %f \n", ln.getFractionUp());

	}

	cl.updateClustering();
	MetaGraph meta = cl.createMetaCluster();
	meta.writeGephiFile("CoarseGrained" + graph.getTypeIdentifier());*/
}
