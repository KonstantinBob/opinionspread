/*
 * MetaGraph.cpp
 *
 *  Created on: Jun 6, 2017
 *      Author: konbob
 */

#include "MetaGraph.h"
#include <armadillo>

/*
 * Default construtor.
 */

MetaGraph::MetaGraph() {

}

/*
 * Constructor by given adjacency matrix.
 */
MetaGraph::MetaGraph(const arma::SpMat<short>& adjacencyMatrix,std::vector<double> nodeWeights,std::vector<double> nodeLabels,const arma::SpMat<double>& edgeWeights) {
	this->adjacencyMatrix = adjacencyMatrix;
	this->nodeWeights = nodeWeights;
	this->nodeLabels = nodeLabels;
	this->edgeWeights = edgeWeights;
}

/*
 * Resizes the network to given size.
 * -> Yields ERROR !
 * Complexity: O(1)
 */
void MetaGraph::resize(const long int& numberOfNodes) {
	std::printf("Not callable. Exit");
	exit(1);
}

const std::vector<double>& MetaGraph::getNodeWeights() const {
	return nodeWeights;
}

void MetaGraph::setNodeWeights(const std::vector<double>& nodeWeights) {
	this->nodeWeights = nodeWeights;
}

/*
 * Writes a spreadsheet file for gephi.
 * It produces two files: One edge table and a node table.
 */
void MetaGraph::writeGephiFile(const std::string& filename) {
	std::ofstream file;

	// === write node table ===
	try {
		// Open file
		file.open("nodes_" + filename, std::ofstream::out);

		// Header
		file << "Id weight fraction" << std::endl;

		//write data to file
		for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {

			// write color value for node
			file << nodeNumber << " " << nodeWeights.at(nodeNumber)<< " "<< nodeLabels.at(nodeNumber) << std::endl;

		}

		// close file
		file.close();

	} catch (...) {
		std::cerr << "Could not open or write to file nodes_" << filename << " or index out of bounds." << std::endl;
	}

	// === write edge table ===
	try {
		// Open file
		file.open("edges_" + filename, std::ofstream::out);

		// Header
		file << "Source Target feature" << std::endl;

		//write data to file
		for (unsigned int nodeNumber = 0; nodeNumber < getNumberOfNodes(); ++nodeNumber) {
			for (unsigned int adjacentNodeNumber = 0; adjacentNodeNumber < getNumberOfNodes(); ++adjacentNodeNumber) {
				// if nodes are connected write edge to file. Self loops are ignored
				if (nodeNumber != adjacentNodeNumber && adjacencyMatrix(nodeNumber, adjacentNodeNumber) != 0) {
					file << nodeNumber << " " << adjacentNodeNumber << " "<< edgeWeights(nodeNumber,adjacentNodeNumber) << std::endl;
				}
			}
		}

		// close file
		file.close();

	} catch (...) {
		std::cerr << "Could not open or write to file edges_" << filename << " or index out of bounds." << std::endl;
	}

	std::cerr << "#Finished writing graph to both files for " << filename << std::endl;

}

const std::vector<double>& MetaGraph::getNodeLabels() const {
return nodeLabels;
}

void MetaGraph::setNodeLabels(const std::vector<double>& nodeLabels) {
this->nodeLabels = nodeLabels;
}
