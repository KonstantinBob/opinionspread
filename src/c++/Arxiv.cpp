/*
 * Arxiv.cpp
 *
 *  Created on: Dec 6, 2016
 *      Author: konbob
 */

#include "Arxiv.h"

/*
 * Constructor.
 * Reads file in SNAP format.
 * Complexity: O(n)
 */
Arxiv::Arxiv(const std::string& filename) {
	// read it from file
	readGraphFromSnapFile(filename);

	// create neighborlists
	createNeighborlists(true);

	// set identifier
	typeIdentifier = "Arxiv";
}


/*
 * Resizes the network to given size.
 * -> Yields ERROR !
 * Complexity: O(1)
 */
void Arxiv::resize(const long int& numberOfNodes) {
	std::printf("Not callable. Exit");
	exit(1);
}
