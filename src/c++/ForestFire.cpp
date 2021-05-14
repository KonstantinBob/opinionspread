/*
 * ForestFire.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#include "ForestFire.h"
#include "RandomNumbers.h"

#include <armadillo>

#include <iostream>
#include <vector>
#include <queue>
#include <set>

/*
 * Constructur. Main algorithm is in createEdges().
 * Complexity: superlinear, but subquadratic.
 */
ForestFire::ForestFire(const double& forwardBuringProbability, const double& backwardBurningProbability, const long int& numberOfNodes) {

	// check number of nodes
	if (numberOfNodes < 2) {
		std::cerr << " Number of nodes must be greater or equal two!" << std::endl;
		exit(1);
	}

	// copy values
	this->forwardBuringProbability = forwardBuringProbability;
	this->backwardBurningProbability = backwardBurningProbability;

	// set identifier
	typeIdentifier = "ForestFire";

	// construct the graph
	createEdges(numberOfNodes);

}

/*
 * Overriding the virtual method for resizing. Main algorithm is in createEdges().
 * Complexity: superlinear, but subquadratic.
 */
void ForestFire::resize(const long int& numberOfNodes) {
	// check number of nodes
	if (numberOfNodes < 2) {
		std::cerr << " Number of nodes must be greater or equal two!" << std::endl;
		exit(1);
	} else {
		// construct new graph
		createEdges(numberOfNodes);
	}
}

/*
 * Creates a forest fire graph with given parameters according to Leskovic et al.
 * Complexity: superlinear, but subquadratic.
 */
void ForestFire::createEdges(const long int& numberOfNodes) {

	// prepare the empty adjacency matrix
	adjacencyMatrix = arma::SpMat<short>(numberOfNodes, numberOfNodes);

	// declare and reserve space for listOfNeighborlists
	listOfInNeighborlists = std::vector<std::vector<long int>>(numberOfNodes, std::vector<long int>(0));
	listOfOutNeighborlists = std::vector<std::vector<long int>>(numberOfNodes, std::vector<long int>(0));

	// loop until desired  number is achieved. The loop starts at 1 since the first node (number 0) is already there.
	for (long int currentNodeNumber = 1; currentNodeNumber < numberOfNodes; ++currentNodeNumber) {

		// pick random node that has already been created
		long int ambassadorNodeNumber = RandomNumbers::getRandomLongInt(0, currentNodeNumber);

		// add new edge to the ambassador node.
		createEdge(currentNodeNumber, ambassadorNodeNumber);

		// update neighborhoods
		updateOutNeighborlists(currentNodeNumber, ambassadorNodeNumber);
		updateInNeighborlists(ambassadorNodeNumber, currentNodeNumber);

		// create the queue of reference nodes
		std::queue<long int> referenceNodeQueue;

		// insert the ambassador node
		referenceNodeQueue.push(ambassadorNodeNumber);

		// create the set of already visited nodes
		std::set<long int> visitedNodeSet;

		//do not visit the current node
		visitedNodeSet.insert(currentNodeNumber);

		do {
			// get the next reference node ...
			long int referenceNodeNumber = referenceNodeQueue.front();

			// ... and remove it from the queue.
			referenceNodeQueue.pop();

			// add reference node to visited node set
			visitedNodeSet.insert(referenceNodeNumber);

			// sample number of neighbors to create
			long int outNeighborsToCreate = RandomNumbers::getRandomGeometric(1.0 - forwardBuringProbability);
			long int inNeighborsToCreate = RandomNumbers::getRandomGeometric(1.0 - backwardBurningProbability);

			// set up list of neighbors to examine
			std::vector<long int> neighborhood;
			neighborhood.reserve(inNeighborsToCreate + outNeighborsToCreate);

			// pick the outgoing edges ...
			if (outNeighborsToCreate > 0) {
				std::vector<long int> outgoingNeighbors = getRandomNeighbors(getAdjacentOutNodes(referenceNodeNumber), outNeighborsToCreate, visitedNodeSet);
				neighborhood.insert(neighborhood.end(), outgoingNeighbors.begin(), outgoingNeighbors.end());
			}

			// ... and the the ingoing edges
			if (inNeighborsToCreate > 0) {
				std::vector<long int> ingoingNeighbors = getRandomNeighbors(getAdjacentInNodes(referenceNodeNumber), inNeighborsToCreate, visitedNodeSet);
				neighborhood.insert(neighborhood.end(), ingoingNeighbors.begin(), ingoingNeighbors.end());
			}

			// loop over neighborhood
			for (long int & nodeNumber : neighborhood) {

				// use node as a reference for new nodes
				referenceNodeQueue.push(nodeNumber);

				// create the edge
				createEdge(currentNodeNumber, nodeNumber);

				// update neighborhoods
				updateOutNeighborlists(currentNodeNumber, nodeNumber);
				updateInNeighborlists(nodeNumber, currentNodeNumber);

				// add node to visited nodes
				visitedNodeSet.insert(nodeNumber);

			}

		} while (referenceNodeQueue.size() > 0);

	}

	// store size
	size = numberOfNodes;
}
