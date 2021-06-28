/*
 * LearningNetwork.h
 *
 *  Created on: Oct 11, 2016
 *      Author: konbob
 */

#ifndef LEARNINGNETWORK_HPP_
#define LEARNINGNETWORK_HPP_

#include "Scientist.h"
#include "RandomNumbers.h"
#include "Graph.h"
#include "StochKronecker.h"
#include "ForestFire.h"
#include "Erdos.h"
#include "Arxiv.h"
#include "Histogram.h"
#include "MetaGraph.h"

#include <armadillo>

#include <string>
#include <vector>
#include <set>
#include <array>
#include <functional>
#include <iostream>
#include <algorithm>
#include <string>

template<typename GraphType> class LearningNetwork {
public:
	// constructor
	LearningNetwork(const GraphType &graphArgument, const double& expGap, const double& q);

	// simulation methods
	inline void performCycle();
	inline void doCommunication(Scientist& scientist);
	inline void doIsingCycle(const double& J, const double& h, const double& beta);
	inline double doIsingCycleWindow(const double& J, const double& h, const double& beta, const long int & lowerBound, const long int & samples);

	// utility methods
	inline long int getSize() const;
	inline double getSimilarity(const long int & scientistNumberI, const long int & scientistNumberJ) const;
	inline bool getOpinion(const long int & scientistNumber) const;

	// graph wrapper
	inline void writeDotFileWithColorData(const std::string& filename, const std::vector<std::string>& colorData, const std::string& appendix) const;
	inline std::vector<long int> getAdjacentAllNodes(const long int& nodeNumber) const;
	inline std::vector<long int> getAdjacentInNodes(const long int& nodeNumber) const;
	inline std::vector<long int> getAdjacentOutNodes(const long int& nodeNumber) const;
	inline unsigned long int getNumberOfEdges() const;

	// analysis methods
	inline std::array<double, 4> calculatePairDistributions() const;
	inline double getFractionUp() const;
	inline double getNumberUp() const;
	inline std::array<double, 2> trackWithDegreeThreshold(const long int& degreeThreshold) const;
	inline MetaGraph createCorrelationGraph(const long int& timeSteps);

private:
	std::vector<Scientist> listOfScientists;
	GraphType graph;
	double flipUp, flipDown, q, m, expGap;

	inline double getCorrelation(const arma::Mat<short>& matrix, const long int& i, const long int & j) const;
};
;

template<typename GraphType>
LearningNetwork<GraphType>::LearningNetwork(const GraphType& graphArgument, const double& expGap, const double& q) :
		graph(graphArgument) {
	// store the size
	long int size = graphArgument.getNumberOfNodes();

	// store set of scientists
	listOfScientists.reserve(size);

	// fill the list
	for (unsigned int index = 0; index < size; ++index) {
		listOfScientists.push_back(Scientist(index, 0.0));

	}

	// store the values
	this->expGap = expGap;
	this->q = q;
	this->m = 1.0 - q;

	// Calculate the flip rates
	flipUp = 1. / (expGap * (1. + expGap));
	flipDown = (expGap) / (1. + expGap);

}

template<typename GraphType>
inline void LearningNetwork<GraphType>::performCycle() {
	for (Scientist& scientist : listOfScientists) {

		//do mutation or communication
		if (RandomNumbers::getRandomReal() < m) {
			scientist.doMutation(expGap);
		} else {
			doCommunication(scientist);
		}
	}
}

template<typename GraphType>
inline void LearningNetwork<GraphType>::doCommunication(Scientist& scientist) {
	// get the id
	long int id = scientist.getId();

	// check if neighborhood is not empty
	if (graph.getNumberOfOutNodes(id) > 0) {

		//pick a random neighbor
		std::vector<long int> neighborhood = graph.getAdjacentOutNodes(id);
		long int randomID = neighborhood.at(RandomNumbers::getRandomLongInt(0, neighborhood.size()));

		// store the opinion
		bool neighborOpinion = listOfScientists.at(randomID).getOpinion();

		// see if the opinions differ
		if (scientist.getOpinion() != neighborOpinion) {
			// ==== gamble for a flip of opinions. ====
			// flip up
			if (neighborOpinion == true && RandomNumbers::getRandomReal() < flipUp) {
				scientist.setOpinion(true);
			}
			// flip down
			if (neighborOpinion == false && RandomNumbers::getRandomReal() < flipDown) {
				scientist.setOpinion(false);
			}
		}
	}
}

/*
 * Returns the number of scientists in the network.
 */
template<typename GraphType>
inline long int LearningNetwork<GraphType>::getSize() const {
	return listOfScientists.size();
}

/*
 * Wrapper for graph type.
 */
template<typename GraphType>
inline void LearningNetwork<GraphType>::writeDotFileWithColorData(const std::string& filename, const std::vector<std::string>& colorData, const std::string& appendix) const {
	graph.writeDotFileWithColorData(filename, colorData, appendix);
}

template<typename GraphType>
inline double LearningNetwork<GraphType>::getSimilarity(const long int& scientistNumberI, const long int& scientistNumberJ) const {
	return Scientist::getSimilarity(listOfScientists.at(scientistNumberI), listOfScientists.at(scientistNumberJ));
}

template<typename GraphType>
inline std::vector<long int> LearningNetwork<GraphType>::getAdjacentAllNodes(const long int& nodeNumber) const {
	return graph.getAdjacentAllNodes(nodeNumber);
}

template<typename GraphType>
inline std::vector<long int> LearningNetwork<GraphType>::getAdjacentInNodes(const long int& nodeNumber) const {
	return graph.getAdjacentInNodes(nodeNumber);
}

template<typename GraphType>
inline std::vector<long int> LearningNetwork<GraphType>::getAdjacentOutNodes(const long int& nodeNumber) const {
	return graph.getAdjacentOutNodes(nodeNumber);
}

/*
 * Returns the opinion of a node.
 */
template<typename GraphType>
inline bool LearningNetwork<GraphType>::getOpinion(const long int& scientistNumber) const {
	return listOfScientists.at(scientistNumber).getOpinion();
}

template<typename GraphType>
inline unsigned long int LearningNetwork<GraphType>::getNumberOfEdges() const {
	return graph.getNumberOfEdges();
}

/*
 * Calculates the averaged distribution of opinions for pairs of nodes.
 */
template<typename GraphType>
inline std::array<double, 4> LearningNetwork<GraphType>::calculatePairDistributions() const {

	// prepare empty array for counting (in order ++,+-,-+,--)
	std::array<double, 4> opinionCounts = { 0., 0., 0., 0. };

	// counter for pairs
	long int pairCounter = 0;

	// loop over all nodes
	for (Scientist scientist : listOfScientists) {

		// get the neighbors
		for (long int neighbor : getAdjacentOutNodes(scientist.getId())) {

			// compare the opinions
			if (scientist.getOpinion() == getOpinion(neighbor)) {
				if (scientist.getOpinion() == true) {
					opinionCounts.at(0)++; // pair is ++
					}
				else {
					opinionCounts.at(3)++; // pair is --
					}
				}
			else {
				if (scientist.getOpinion() == true) {
					opinionCounts.at(1)++; // pair is +-
					}
				else {
					opinionCounts.at(2)++; // pair is -+
					}
				}

				// increment counter
			pairCounter++;
		}
	}

	// divide by pairCounter
	for (unsigned long int index = 0; index < opinionCounts.size(); ++index) {
		opinionCounts.at(index) /= (pairCounter);
	}

	return opinionCounts;
}

/*
 * Performs a simple single flip Ising model simulation.
 * Convention: Spin up == True, Spin down == Down
 * Here, graphs are considered UNDIRECTED!
 */
template<typename GraphType>
inline void LearningNetwork<GraphType>::doIsingCycle(const double& J, const double& h, const double& beta) {

	// loop over nodes
	for (Scientist& scientist : listOfScientists) {

		//count same and different opinions in neighborhood
		double counterSame = 0.0;
		double counterDiff = 0.0;
		for (long int & neighbor : getAdjacentAllNodes(scientist.getId())) {
			if (listOfScientists.at(neighbor).getOpinion() == scientist.getOpinion()) {
				++counterSame;
			} else {
				++counterDiff;
			}
		}

		// for external field coupling convert  opinion to -1 or 1
		double spinOrientation = 0.0;
		if (scientist.getOpinion() == true) {
			spinOrientation = 1.0;
		} else {
			spinOrientation = -1.0;
		}

		// calculate the energy difference of the flip
		double energyDiff = 2.0 * J * (counterSame - counterDiff) + 2.0 * h * spinOrientation;

		if (false) {
			printf("Spin %ld is %f and has %f same neighbors and %f diff neighbors. dE =%f, Pflip=%f", scientist.getId(), spinOrientation, counterSame, counterDiff, energyDiff,
					fmin(1., exp(-beta * energyDiff)));
		}

		// decide whether to accept flip or not
		if (energyDiff < 0) {
			// accept flip instantly
			scientist.setOpinion(!scientist.getOpinion());

		} else if (RandomNumbers::getRandomReal() < exp(-beta * energyDiff)) { // do Metropolis
		//also accept flip
			scientist.setOpinion(!scientist.getOpinion());
		}
		if (false) {
			printf(", new orientation is %d \n", scientist.getOpinion());
		}
	}
}

/*
 * Performs a move for a successive umbrella  Ising model simulation.
 * Convention: Spin up == True, Spin down == Down
 * Here, graphs are considered UNDIRECTED!
 * lowerBound is the minimum number of up Spins allowed, and samples is the number of runs for the simulation.
 * Returns the logarithm of the ratio of N(upper)/N(lower)
 */
template<typename GraphType>
inline double LearningNetwork<GraphType>::doIsingCycleWindow(const double& J, const double& h, const double& beta, const long int& lowerBound, const long int& sampels) {
	// flag for print out
	const bool debug = false;

	// prepare counter for samples, where N_up=lowerBound+1
	double counterUp = 0;

	// === prepare network in right state

	// find distinct random indices....
	std::set<long int> indexSet;
	while (indexSet.size() < lowerBound) {
		// get a random index for a scientist
		long int randomIndex = RandomNumbers::getRandomLongInt(0,getSize());

		// store it if it is not yet inside the set
		if (indexSet.count(randomIndex) == 0) {
			indexSet.insert(randomIndex);
		}


	}

	// .... and set those with opinion up, the others with opinion down
	for (Scientist& scientist : listOfScientists) {
		// check if id is in index set
		if (indexSet.count(scientist.getId()) == 1) {
			scientist.setOpinion(true);
			if (debug) {
				printf("Set up  %ld \n",scientist.getId());
			}
		} else {
			scientist.setOpinion(false);
		}
	}

	if (debug) {
		printf("Number up after init %f \n",getFractionUp()*getSize());
	}

	// sample loop
	for (long int index = 0; index < sampels; ++index) {

		// loop over nodes
		for (Scientist& scientist : listOfScientists) {

			//count same and different opinions in neighborhood
			double counterSame = 0.0;
			double counterDiff = 0.0;
			for (long int & neighbor : getAdjacentAllNodes(scientist.getId())) {
				if (listOfScientists.at(neighbor).getOpinion() == scientist.getOpinion()) {
					++counterSame;
				} else {
					++counterDiff;
				}
			}

			// for external field coupling convert  opinion to -1 or 1
			double spinOrientation = 0.0;
			if (scientist.getOpinion() == true) {
				spinOrientation = 1.0;
			} else {
				spinOrientation = -1.0;
			}

			// calculate the energy difference of the flip
			double energyDiff = 2.0 * J * (counterSame - counterDiff) + 2.0 * h * spinOrientation;


			if (debug) {
				printf("Spin %ld is %f and has %f same neighbors and %f diff neighbors. dE =%f, P_flip=%f ", scientist.getId(), spinOrientation, counterSame, counterDiff,
						energyDiff, fmin(1., exp(-beta * energyDiff)));
			}

			// decide whether to accept flip or not
			if (energyDiff < 0) {

				// flip accept for energy...
				scientist.setOpinion(!scientist.getOpinion());

				// but check now umbrella sampling condition
				long int numberUp = getNumberUp();
				if (numberUp < lowerBound || numberUp > lowerBound + 1) {
					// reverse the change
					scientist.setOpinion(!scientist.getOpinion());
				}

			} else if (RandomNumbers::getRandomReal() < exp(-beta * energyDiff)) { // do Metropolis
			//also accept flip
				scientist.setOpinion(!scientist.getOpinion());

				// but check now umbrella sampling condition
				long int numberUp = getNumberUp();
				if (numberUp < lowerBound || numberUp > lowerBound + 1) {
					// reverse the change
					scientist.setOpinion(!scientist.getOpinion());
				}
			}
			if (debug) {
				printf(", new orientation is %d \n", scientist.getOpinion());
			}
		}

		// check counter
		if ( getNumberUp() == lowerBound+1 ) {
			++counterUp;
		}

		if (debug) {
			printf("counter %f samples %ld ratio %f \n",counterUp,sampels,log((counterUp)/(sampels-counterUp)));
		}
	}

	// prevent nans
	if (counterUp == 0.0) {
		counterUp = 1.0;
	}

	// return logarithm of ratios
	return log((counterUp)/(sampels-counterUp));
}

template<typename GraphType>
inline double LearningNetwork<GraphType>::getFractionUp() const {
	// counter for up spins
	double counterUp = 0.0;

	// loop over scientists
	for (const Scientist& scientist : listOfScientists) {
		// check opinion
		if (scientist.getOpinion() == true) {
			++counterUp;
		}
	}

	// normalize the result
	double result = counterUp / getSize();

	return result;
}

template<typename GraphType>
inline double LearningNetwork<GraphType>::getNumberUp() const {
	// counter for up spins
	double counterUp = 0.0;

	// loop over scientists
	for (const Scientist& scientist : listOfScientists) {
		// check opinion
		if (scientist.getOpinion() == true) {
			++counterUp;
		}
	}

	// return the result
	return counterUp ;

}

/*
 * Prints out the fraction of nodes up for two groups:
 * One above and equal to the degree threshold and one below the threshold.
 * Return value:
 * <fracAbove,fracBelow>
 * Complexity: O(n)
 */
template<typename GraphType>
inline std::array<double, 2> LearningNetwork<GraphType>::trackWithDegreeThreshold(const long int& degreeThreshold) const {
	// prepare counters
	double counterUpAbove = 0.;
	double counterUpBelow = 0.;

	double counterAbove = 0.0;
	double counterBelow = 0.0;

	for (const Scientist& scientist : listOfScientists) {
		// check the degree threshold
		if (getAdjacentOutNodes(scientist.getId()).size() >= degreeThreshold) {

			// is above threshold,...
			++counterAbove;

			//... check opinion
			if (scientist.getOpinion() == true) {
				++counterUpAbove;
			}

		} else {

			// is below threshold,...
			++counterBelow;

			// check opinion
			if (scientist.getOpinion() == true) {
				++counterUpBelow;
			}

		}
	}

	// prepare output
	std::array<double, 2> results;

	// calculate fractions
	results.at(0) = counterUpAbove / counterAbove;
	results.at(1) = counterUpBelow / counterBelow;

	return results;
}

/*
 * Creates an MetaGraph with correlation coefficients as edge weights.
 */
template<typename GraphType> inline MetaGraph LearningNetwork<GraphType>::createCorrelationGraph(const long int& timeSteps) {

	// prepare matrix for storing the state of the nodes
	arma::Mat<short> opinionMatrix(timeSteps, getSize());

	// prepare matrix for storing edge weights
	arma::SpMat<double> edgeWeights(getSize(), getSize());

	// === get the data ===

	// Loop for time evolution
	for (long int timeStep = 0; timeStep < timeSteps; ++timeStep) {

		// perform a cycle
		performCycle();

		// store the opinion of every node
		for (long int nodeNumber = 0; nodeNumber < getSize(); ++nodeNumber) {
			if (listOfScientists.at(nodeNumber).getOpinion() == true) {
				opinionMatrix(timeStep, nodeNumber) = 1;
			} else {
				opinionMatrix(timeStep, nodeNumber) = 0;
			}
		}
	}

	//  === do the analysis  ===

	// loop over edges
	for (long int nodeNumber = 0; nodeNumber < getSize(); ++nodeNumber) {
		for (long int neighbor : getAdjacentOutNodes(nodeNumber)) {

			// calculate the correlation
			edgeWeights(nodeNumber, neighbor) = getCorrelation(opinionMatrix, nodeNumber, neighbor);
		}
	}

	// === return result ===
	std::vector<double> ones(getSize(), 1.0);

	return MetaGraph(graph.getAdjacencyMatrix(), ones, ones, edgeWeights);
}

/*
 * Calculation of the correlation between rows i and j in matrix by library call.
 */
template<typename GraphType> inline double LearningNetwork<GraphType>::getCorrelation(const arma::Mat<short>& matrix, const long int& i, const long int& j) const {

	// prepare vectors
	arma::vec vecI = arma::conv_to<arma::vec>::from(matrix.col(i));
	arma::vec vecJ = arma::conv_to<arma::vec>::from(matrix.col(j));

	//call method
	arma::vec result = arma::cor(vecI, vecJ);

	//return result
	return result(0, 0);
}

#endif /* LEARNINGNETWORK_HPP_ */

