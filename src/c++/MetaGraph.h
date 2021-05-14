/*
 * MetaGraph.h
 *
 *  Created on: Jun 6, 2017
 *      Author: konbob
 */

#ifndef METAGRAPH_H_
#define METAGRAPH_H_

#include "Graph.h"
#include <armadillo>


class MetaGraph: public Graph {
public:
	// methods
	MetaGraph();
	MetaGraph(const arma::SpMat<short>& adjacencyMatrix,std::vector<double> nodeWeights,std::vector<double> nodeLabels,const arma::SpMat<double>& edgeWeights);
	void resize(const long int& numberOfNodes);
	const std::vector<double>& getNodeWeights() const;
	void setNodeWeights(const std::vector<double>& nodeWeights);
	void writeGephiFile(const std::string& filename);
	const std::vector<double>& getNodeLabels() const;
	void setNodeLabels(const std::vector<double>& nodeLabels);

private:
	// members
	std::vector<double> nodeWeights;
	std::vector<double> nodeLabels;
	arma::SpMat<double> edgeWeights;

};

#endif /* METAGRAPH_H_ */
