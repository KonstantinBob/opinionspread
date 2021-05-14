/*
 * RandomNumbers.h
 *
 *  Created on: Sep 19, 2016
 *      Author: konbob
 */

#ifndef RANDOMNUMBERS_H_
#define RANDOMNUMBERS_H_
#include <random>

class RandomNumbers {
private:
	static std::mt19937 generator;
public:
	static double getRandomReal();
	static double getNormalReal(double mu, double sigma);
	static int getRandomInt(int a, int b);
	static long int getRandomLongInt(long int a, long int b);
	static long int getRandomGeometric(double p);
	static bool getRandomBool();
	template<typename T> static T getRandomEntry(const std::vector<T>& list);

};

/*
 * Returns a random entry of the given list.
 */

template<typename T> inline T RandomNumbers::getRandomEntry(const std::vector<T>& list) {
	// get random index
	long int index = getRandomLongInt(0, list.size());

	// return the element
	return list.at(index);
}
#endif /* RANDOMNUMBERS_H_ */
