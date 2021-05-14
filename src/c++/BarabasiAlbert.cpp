/*
 * BarabasiAlbert.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: konbob
 */

#include "BarabasiAlbert.h"
#include "RandomNumbers.h"

BarabasiAlbert::BarabasiAlbert(const long int & mZero, const long int & m, const long int& numberOfNodes) {

	if (0 < mZero && 0 < m && m <= mZero) {
		// store the values
		this->m = m;
		this->mZero = mZero;
	} else {
		std::cerr << "Illegal parameter m or mZero!" << std::endl;
		exit(1);
	}


	// check number of nodes
	if (numberOfNodes < 2 || numberOfNodes < mZero) {
		std::cerr << " Number of nodes must be greater or equal two and greater than mZero!" << std::endl;
		exit(1);
	} else {
		// construct new graph
		createEdges(numberOfNodes);
	}

	// set identifier
	typeIdentifier = "BarabasiAlbert";
}

void BarabasiAlbert::resize(const long int& numberOfNodes) {

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
 * Creates a directed graph according to the Barabasi-Albert model.
 */
void BarabasiAlbert::createEdges(const long int& numberOfNodes) {
	// set the final size
	size = numberOfNodes;



	// prepare the empty adjacency matrix
	adjacencyMatrix = arma::SpMat<short>(numberOfNodes, numberOfNodes);

	// === create mZero full connected nodes ====
	for (long int currentNodeNumber = 0; currentNodeNumber < mZero; ++currentNodeNumber) {
		for (long int currentNeighbor = 0; currentNeighbor < mZero; ++currentNeighbor) {
			createEdge(currentNodeNumber, currentNeighbor);
		}
		// remove  self loops
		removeEdge(currentNodeNumber, currentNodeNumber);

	}

	// prepare the neighbor lists
	createNeighborlists(false);

	// ==== create other nodes ====

	// list for storage of previous target nodes (idea by networkx implementation)
	std::vector<long int> targetNodes;

	// reserve for complete graph and edges to come
	targetNodes.reserve(mZero * mZero + (numberOfNodes - mZero) * m);

	// fill list for complete graph
	for (long int currentNodeNumber = 0; currentNodeNumber < mZero; ++currentNodeNumber) {
		for (long int counter = 0; counter < mZero - 1; ++counter) {
			targetNodes.push_back(currentNodeNumber);
		}
	}

	// loop until desired  number is achieved. The loop starts at mZero since the first nodes already there.
	for (long int currentNodeNumber = mZero; currentNodeNumber < numberOfNodes; ++currentNodeNumber) {

		// counter for valid trails
		long int counterValid = 0;

		// form m valid edges
		do {
			// pick random node that has already been created
			long int targetNodeNumber = RandomNumbers::getRandomEntry(targetNodes);

			// check if there is not already a connection
			if (adjacencyMatrix(currentNodeNumber, targetNodeNumber) == 0) {

				// increment counter
				++counterValid;

				// add new edge to the target node
				createEdge(currentNodeNumber, targetNodeNumber);

				// update neighborhoods
				updateOutNeighborlists(currentNodeNumber, targetNodeNumber);
				updateInNeighborlists(targetNodeNumber, currentNodeNumber);

				// push back current node
				targetNodes.push_back(currentNodeNumber);
			}


		} while (counterValid < m);

	}

}
