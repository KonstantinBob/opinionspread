/*
 * Erdos.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#include "Erdos.h"
#include "RandomNumbers.h"

#include <armadillo>

/*
 * Constructor for Erdos-Renyi graph.
 * Arguments:
 * - edge formation probability
 * - desired number of nodes
 * Complexity: O(numberOfNodes^2)
 */
Erdos::Erdos(const double& probability, const long int& numberOfNodes) {
	//check arguments
	if (probability < 0.0 || probability > 1.0 || numberOfNodes < 1.0) {
		std::cerr << "Invalid arguments for createErdosRenyiGraph! Exit." << std::endl;
		exit(1);
	} else {
		// store the parameter
		this->probability = probability;

		// set identifier
		typeIdentifier = "Erdos";

		// create the graph
		createEdges(numberOfNodes);

	}
}

/*
 * Resize the Erdos-Renyi graph.
 * Arguments:
 * - desired number of nodes
 * Complexity: O(numberOfNodes^2)
 */
void Erdos::resize(const long int& numberOfNodes) {
	// check arguments
	if (numberOfNodes < 1.0) {
		std::cerr << "Invalid arguments for createErdosRenyiGraph! Exit." << std::endl;
		exit(1);
	} else {
		createEdges(numberOfNodes);
	}

}

/*
 * Create the edges in the graph.
 * Arguments:
 * - desired number of nodes
 * Complexity: O(numberOfNodes^2)
 */
void Erdos::createEdges(const long int& numberOfNodes) {
	// prepare the empty adjacency matrix
	adjacencyMatrix = arma::SpMat<short>(numberOfNodes, numberOfNodes);
	// ========= loop over it and create the edges ===========
	for (long int columnIndex = 0; columnIndex < numberOfNodes; columnIndex++) {
		for (long int rowIndex = 0; rowIndex < numberOfNodes; rowIndex++) {
			// sample if the edge will be placed
			if (rowIndex != columnIndex && RandomNumbers::getRandomReal() < probability) {
				createEdge(columnIndex, rowIndex);

			}

		}
	}

	// store the size
	size = numberOfNodes;

	// store neighborhoods
	createNeighborlists(true);

}
