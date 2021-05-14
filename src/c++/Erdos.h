/*
 * Erdos.h
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#ifndef ERDOS_H_
#define ERDOS_H_

#include "Graph.h"

class Erdos: public Graph {
public:
	Erdos(const double& probability,const long int& numberOfNodes);
	virtual void resize(const long int& numberOfNodes);

private:
	double probability;

	void createEdges(const long int& numberOfNodes);
};

#endif /* ERDOS_H_ */
