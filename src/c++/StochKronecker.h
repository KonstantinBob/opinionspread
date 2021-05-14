/*
 * StochKronecker.h
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#ifndef STOCHKRONECKER_H_
#define STOCHKRONECKER_H_

#include "Graph.h"


class StochKronecker: public Graph {
public:
	StochKronecker(const arma::mat& IteratorMatrix,const long int& numberOfNodes);
	virtual void resize(const long int& numberOfNodes);


private:
	void createDeterministicKroneckerGraph(const arma::mat& IteratorMatrix, const int& power);
	void createEdges(const long int& numberOfNodes);
	arma::mat IteratorMatrix;
};

#endif /* STOCHKRONECKER_H_ */
