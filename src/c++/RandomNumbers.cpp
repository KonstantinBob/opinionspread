/*
 * RandomNumbers.cpp
 *
 *  Created on: Sep 19, 2016
 *      Author: konbob
 */

#include "RandomNumbers.h"
#include <math.h>
#include <iostream>
#include <random>

std::mt19937 RandomNumbers::generator; // period of 2E19937

/*
 * returns random real in interval [0,1[
 */
double RandomNumbers::getRandomReal() {
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}

/*
 * returns normal distributed random real
 */
double RandomNumbers::getNormalReal(double mu, double sigma) {
	std::normal_distribution<double> distribution(mu, sigma);
	return distribution(generator);
}

/*
 * Returns a random integer in [a,b[
 */
int RandomNumbers::getRandomInt(int a, int b) {
	std::uniform_int_distribution<int> distribution(a, b-1);
	return distribution(generator);
}

/*
 * Returns a random long integer in [a,b[
 */
long int RandomNumbers::getRandomLongInt(long int a,long  int b) {
	std::uniform_int_distribution<long int> distribution(a, b-1);
	return distribution(generator);
}

/*
 * Returns random geometric number with
 * Prob(X = k) = (1-p)^{k}p
 */
long int RandomNumbers::getRandomGeometric(double p){
	std::geometric_distribution<int> distribution(p);
	return distribution(generator);
}

/*
 * Returns a random bool, based on a random int.
 */
bool RandomNumbers::getRandomBool() {
	// get random int
	int randomInt = getRandomInt(0,2);

	// now use 0 -> false, 1-> true;
	if (randomInt == 0) {
		return false;
	}
	return true;
}







