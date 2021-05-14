/*
 * WattsStrogatz.cpp
 *
 *  Created on: Jul 19, 2017
 *      Author: konbob
 */

#include "WattsStrogatz.h"
#include "RandomNumbers.h"

#include <armadillo>

/*
 * Constructor for WattsStrogatz network.
 */
WattsStrogatz::WattsStrogatz(const long int& k, const double& pRewire, const long int& numberOfNodes) {
	// check k input
	if (k > 0) {
		this->k = k;
	} else {
		std::cerr << "# k must be greater than zero! Exit." << std::endl;
		exit(1);
	}

	// check p input
	if (0. <= pRewire && pRewire <= 1.0) {
		this->pRewire = pRewire;
	} else {
		std::cerr << "# p must be between zero and one! Exit." << std::endl;
		exit(1);
	}

	// check nodes input
	if (numberOfNodes > k) {
		this->size = numberOfNodes;
	} else {
		std::cerr << "# numberOfNodes must be greater than k. Exit." << std::endl;
		exit(1);
	}

	// create the edges
	createEdges(numberOfNodes);

	// set identifier
	typeIdentifier = "WattsStrogatz";
}

/*
 * Resize method. Calls create Edges.
 */
void WattsStrogatz::resize(const long int& numberOfNodes) {
	// check nodes input
	if (numberOfNodes > k) {
		createEdges(numberOfNodes);
	} else {
		std::cerr << "# numberOfNodes must be greater than k! Exit." << std::endl;
		exit(1);
	}
}

/*
 * Creates the new edges to the given number of nodes.
 */
void WattsStrogatz::createEdges(const long int& numberOfNodes) {
	// check argument
	if (numberOfNodes < k) {
		std::cerr << "# numberOfNodes must be greater than k! Exit." << std::endl;
		exit(1);
	}

	// prepare the adjacency matrix
	adjacencyMatrix = arma::SpMat<short>(numberOfNodes, numberOfNodes);

	// ==== build up ring structure
	// connect every node to its k nearest neighbors by DIRECTED edges.
	for (long int nodeNumber = 0; nodeNumber < numberOfNodes; ++nodeNumber) {

		for (long int l = 1; l < k + 1; ++l) {
			createEdge(nodeNumber, (nodeNumber + l) % numberOfNodes);
		}
	}

	// === rewire the edges
	// loop over the distance
	for (long int l = 1; l < k + 1; ++l) {

		// loop over the nodes
		for (long int nodeNumber = 0; nodeNumber < numberOfNodes; ++nodeNumber) {
			// === find unconnected node

			// flag for finding unconnected node
			bool foundUnconnected = false;

			// storage for target node
			long int targetNode = -1;

			do {
				// pick random node
				targetNode = RandomNumbers::getRandomLongInt(0, numberOfNodes);

				// check if it is not already connected
				if (adjacencyMatrix(nodeNumber, targetNode) == 0) {
					foundUnconnected = true;
				}
			} while (!foundUnconnected);

			// === check rewiring

			if (RandomNumbers::getRandomReal() < pRewire) {
				// rewire the edge
				createEdge(nodeNumber, targetNode);

				// remove old edge
				removeEdge(nodeNumber, (nodeNumber + l) % numberOfNodes);
			}
		}
	}

	// create neighbor lists
	createNeighborlists(false);
}

