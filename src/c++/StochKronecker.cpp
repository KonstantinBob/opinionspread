/*
 * StochKronecker.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: konbob
 */

#include "StochKronecker.h"
#include "RandomNumbers.h"

#include <armadillo>

#include <cmath>

StochKronecker::StochKronecker(const arma::mat& IteratorMatrix, const long int& numberOfNodes) {
	//check arguments
	if (numberOfNodes < 1.0) {
		std::cerr << "Invalid arguments for creating a stochastic Kronecker Graph! Exit." << std::endl;
		exit(1);
	} else {
		// store the parameters
		this->IteratorMatrix = IteratorMatrix;

		// set identifier
		typeIdentifier = "StochKron";

		// create the graph
		createEdges(numberOfNodes);

	}

}

void StochKronecker::resize(const long int& numberOfNodes) {
	//check arguments
	if (numberOfNodes < 1.0) {
		std::cerr << "Invalid arguments for creating a stochastic Kronecker Graph! Exit." << std::endl;
		exit(1);
	} else {

		// create the graph
		createEdges(numberOfNodes);

	}

}

void StochKronecker::createDeterministicKroneckerGraph(const arma::mat& IteratorMatrix, const int& power) {
	// check if iterator is a square matrix
	if (IteratorMatrix.n_rows == IteratorMatrix.n_cols && power > 1) {
		// store matrix size
		int matrixSize = IteratorMatrix.n_rows;
		// calculate final matrix size
		int finalSize = (int) pow(matrixSize, power);
		// create new matrix
		arma::SpMat<short> temporaryMatrix(finalSize, finalSize);

		// fill matrix
		for (int rowCount = 0; rowCount < finalSize; ++rowCount) {
			for (int columnCount = 0; columnCount < finalSize; ++columnCount) {
				// start from big matrix and go down recursively
				int subsize = finalSize; // current size of the submatrix in the current stage
				bool setEntry = true; // tells if entry will be nonzero
				int remainRow = rowCount; // stores the remainder of the row coordinate of the entry
				int remainColumn = columnCount; // stores the remainder of the column coordinate of the entry

				for (int stage = power; stage > 0; --stage) {
					// calculate size of submatrix
					subsize /= matrixSize;
					// get entry of iterator matrix
					int subrow = remainRow / subsize; // row coordinate in the submatrix of the current stage
					int subcolumn = remainColumn / subsize; // row coordinate in the submatrix of the current stage
					double subentry = IteratorMatrix(subrow, subcolumn); // entry of the submatrix in the current stage
					// calculate remaining coordinates
					remainColumn = remainColumn % subsize;
					remainRow = remainRow % subsize;

					if (subentry == 0) {
						setEntry = false;
						break;
					}
				}
				if (setEntry) {
					temporaryMatrix(rowCount, columnCount) = 1.0;

				} else {

				}

			}
		}
		adjacencyMatrix = temporaryMatrix;

		// prepare neighborlists
		createNeighborlists(true);


	} else {
		std::cout << "Iterator matrix must be a square matrix and power > 1 !" << std::endl;
		exit(1);
	}
}

void StochKronecker::createEdges(const long int& numberOfNodes) {
	//store the new size
	size = numberOfNodes;

	// get size of iterator matrix
	int iteratorMatrixSize = IteratorMatrix.n_rows;
	int iteratorMatrixEntryNumber = iteratorMatrixSize * iteratorMatrixSize;

	// calculate the power of the matrix
	// --> Creating bigger matrix by rounding up, the size will be reduced later!
	double power = ceil(log(numberOfNodes)/log(iteratorMatrixSize));

	// calculate final matrix size
	int finalSize = (int) pow(iteratorMatrixSize, power);

	// create new matrix
	arma::SpMat<short> temporaryMatrix(finalSize, finalSize);

	// get sum of all entries in Iterator
	double sumIteratorEntries = 0.0;

	// cumulative entries
	std::vector<double> cumulativeEntries(iteratorMatrixEntryNumber, 0.0);
	for (int i = 0; i < iteratorMatrixSize; ++i) {
		for (int j = 0; j < iteratorMatrixSize; ++j) {
			sumIteratorEntries += IteratorMatrix(i, j);
			try {
				cumulativeEntries.at(i * iteratorMatrixSize + j) = sumIteratorEntries;

			} catch (...) {
				std::cout << "Out of Range error at createStochasticKroneckerGraph() " << std::endl;
			}
		}
	}

	// normalize cumulative entries
	for (int i = 0; i < iteratorMatrixSize; ++i) {
		for (int j = 0; j < iteratorMatrixSize; ++j) {
			cumulativeEntries.at(i * iteratorMatrixSize + j) /= sumIteratorEntries;
		}
	}

	// Calculate number of expected Edges
	int ExpectedEdges = (int) pow(sumIteratorEntries, power);

	// counter of placed edges
	int placedEdgesCounter = 0;

	while (placedEdgesCounter < ExpectedEdges) {

		// reset start values
		int rowCoordinate = 0; // accumulates the row coordinate of the edge
		int columnCoordinate = 0; // accumulates the column coordinate of the edge
		double currentLengthScale = pow(iteratorMatrixSize, power);

		// go down recursively and decide where the edge will be placed
		for (int stage = power; stage > 0; --stage) {
			currentLengthScale /= iteratorMatrixSize;
			bool foundEntry = false;
			// sample on iterator matrix entries
			double randomNumber = RandomNumbers::getRandomReal();
			for (int i = 0; i < iteratorMatrixSize && !foundEntry; ++i) {
				for (int j = 0; j < iteratorMatrixSize && !foundEntry; ++j) {
					try {

						if (randomNumber <= cumulativeEntries.at(i * iteratorMatrixSize + j)) {
							// Right entry was found,
							foundEntry = true;
							// accumulate coordinates
							rowCoordinate += currentLengthScale * i;
							columnCoordinate += currentLengthScale * j;

						} // else the edge will not be placed and the loop can be left
					} catch (...) {
						std::cout << "Out of Range error at createStochasticKroneckerGraph() " << std::endl;
					}
				}
			}
		}

		if (temporaryMatrix(rowCoordinate, columnCoordinate) == 0) { // check if edge has not already been placed
			temporaryMatrix(rowCoordinate, columnCoordinate) = 1;
			placedEdgesCounter++;
		}

	}

	// copy only the desired number of nodes
	adjacencyMatrix = temporaryMatrix.submat( 0,0, size-1, size -1);

	createNeighborlists(true);
}
