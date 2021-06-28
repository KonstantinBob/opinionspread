/*
 * Clustering.h
 *
 *  Created on: Jan 12, 2017
 *      Author: konbob
 */

#ifndef CLUSTERING_HPP_
#define CLUSTERING_HPP_

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <armadillo>

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <map>

#include "LearningNetwork.hpp"
#include "Histogram.h"
#include "MetaGraph.h"

typedef boost::bimap<boost::bimaps::set_of<long int>, boost::bimaps::multiset_of<long int>> multiBiMap_t;
typedef multiBiMap_t::value_type relation;

template<typename GraphType> class Clustering {
public:
	inline Clustering(const LearningNetwork<GraphType>& learningNetwork);

	// output
	inline void printClustering(const std::string& filename) const;
	inline void writeClusteringFileScikit(const std::string& filename) const;
	inline Histogram getSizeDistribution(const double& binWidth) const;

	// clustering methods
	inline void doFullClustering();
	inline void doEntropyOptimization(const long int& maxIterations);
	inline void doDissimilarityOptimization(const long int& maxIterations);
	inline void doEntropySplitting(const long int& clusterNumber);
	inline double calculateTotalModularity();
	inline double calculateTotalEntropy() const;
	inline void doModularityOptimization(const long int& maxIterations);
	inline double getOrderParameter(const long int& numberClustersStart) const;
	inline void updateClustering();

	// coarse graining methods
	inline MetaGraph createGraphOfClusterNodes(const long int& clusterNumber) const;
	inline MetaGraph createMetaCluster() const;

	// getter
	inline long int getNumberOfClusters() const;
	inline std::vector<long int> getNodes(const long int& clusterNumber) const;
	inline long int getNumberOfNodes(const long int& clusterNumber) const;
	inline std::vector<long int> getListOfNonEmptyClusters() const;

	// analysis
	inline double getUpFractionOfCluster(const long int& clusterNumber) const;
	inline void writeNodesAndClusters(const std::string& filename, const long int& timeStep) const;

private:
	// utility methods
	inline std::vector<long int> getNodesAndNeighbors(const long int& clusterNumber) const;
	inline long int getUnconnectedRandomCluster(const long int& sourceCluster, bool& flag);
	inline double getClusterEntropy(const long int& cluster) const;
	inline double getClusterDissimilarity(const long int& cluster) const;
	inline bool checkNodeExistence(const long int& node, const long int& clusterNumber) const;
	inline double getModularity(const long int& clusterNumber) const;
	inline std::vector<long int> getConnectedClusters(const long int& nodeNumber) const;
	inline void removeNode(const long int& nodeNumber, const long int& clusterNumber);
	inline void insertNode(const long int& nodeNumber, const long int& clusterNumber);
	inline long int getRandomNode(const long int& clusterNumber) const;
	inline long int getRandomCluster();


	// member
	multiBiMap_t clustering;
	multiBiMap_t clusteringModularity;
	const LearningNetwork<GraphType>& learningNetwork;
	arma::mat modularityMatrix;

};

/*
 * Constructor for Clustering.
 */
template<typename GraphType>
inline Clustering<GraphType>::Clustering(const LearningNetwork<GraphType>& learningNetwork) :
		learningNetwork(learningNetwork) {
	// create a cluster for every node
	for (long int var = 0; var < learningNetwork.getSize(); ++var) {
		insertNode(var, var);
	}

	// do the modularity optimization
	doModularityOptimization(3 * learningNetwork.getSize());

	// store the modulization based clustering
	clusteringModularity = clustering;
}

/*
 * Returns the number of clusters that hold a node.
 * Complexity: O(n)
 */
template<typename GraphType>
inline long int Clustering<GraphType>::getNumberOfClusters() const {
	// get the list
	std::vector<long int> listOfKeys = getListOfNonEmptyClusters();

	// return its size
	return listOfKeys.size();
}

/*
 * Returns the node numbers in a given cluster.
 * In case the cluster number does not point to a existing cluster an empty vector is returned.
 */
template<typename GraphType>
inline std::vector<long int> Clustering<GraphType>::getNodes(const long int& clusterNumber) const {
	auto it_low = clustering.right.lower_bound(clusterNumber);
	auto it_high = clustering.right.upper_bound(clusterNumber);

	// set up vector
	std::vector<long int> returnVector;
	returnVector.reserve(std::distance(it_low, it_high));

	// fill it
	for (; it_low != it_high; ++it_low) {
		returnVector.push_back(it_low->second);
	}

	return returnVector;
}

/*
 * Returns the number of node in the given cluster.
 */
template<typename GraphType>
inline long int Clustering<GraphType>::getNumberOfNodes(const long int& clusterNumber) const {
	auto it_low = clustering.right.lower_bound(clusterNumber);
	auto it_high = clustering.right.upper_bound(clusterNumber);

	return std::distance(it_low, it_high);
}

/*
 * Returns a random node.
 */
template<typename GraphType>
inline long int Clustering<GraphType>::getRandomNode(const long int& clusterNumber) const {
	// get the iterators
	auto it_low = clustering.right.lower_bound(clusterNumber);
	auto it_high = clustering.right.upper_bound(clusterNumber);

	// calculate distances
	long int distance = std::distance(it_low, it_high);

	// advance it
	long int randomDisplacement = RandomNumbers::getRandomLongInt(0, distance);
	std::advance(it_low, randomDisplacement);
	return it_low->second;
}

/*
 * Creates a graph with the clusters as nodes and the number of nodes therein as weights and the fraction of up nodes as labels.
 */
template<typename GraphType>
inline MetaGraph Clustering<GraphType>::createMetaCluster() const {
	// get list of clusters. It is already sorted.
	std::vector<long int> clusterList = getListOfNonEmptyClusters();


	// prepare the empty adjacency matrix
	arma::SpMat<short> adjacencyMatrix(clusterList.size(), clusterList.size());

	// prepare the weight list
	std::vector<double> weightList(clusterList.size(),0.0);

	// prepare the label list
	std::vector<double> labelList(clusterList.size(),0.0);

	// ==== For compressing the adjacency matrix the cluster numbers need to be mapped onto 0,...,Size-1 ====
	// create map object
	std::map<long int, long int> newClusterNumbers;

	// fill it
	for (long int index = 0; index < clusterList.size(); index++) {
		newClusterNumbers[clusterList.at(index)] = index;
	}

	// === loop over clusters ===
	for (long int& cluster : getListOfNonEmptyClusters()) {

		// list of connected clusters
		std::vector<long int> connectedClusters;

		// reserve some capacity
		connectedClusters.reserve(getNumberOfClusters());

		// loop over its nodes
		for (long int & node : getNodes(cluster)){
			// get the connected clusters
			std::vector<long int> connections = getConnectedClusters(node);

			// append the the data
			connectedClusters.insert(connectedClusters.end(),connections.begin(),connections.end());
		}

		// get unique cluster numbers
		std::sort(connectedClusters.begin(), connectedClusters.end());
		connectedClusters.erase(std::unique(connectedClusters.begin(), connectedClusters.end()), connectedClusters.end());

		// add the edges to the adjacency matrix
		for (long int&  neighborCluster : connectedClusters) {
			adjacencyMatrix(newClusterNumbers.find(cluster)->second,newClusterNumbers.find(neighborCluster)->second)=1;
		}

		// save the weight
		weightList.at(newClusterNumbers.find(cluster)->second) = getNumberOfNodes(cluster);

		// save the label
		labelList.at(newClusterNumbers.find(cluster)->second) = getUpFractionOfCluster(cluster);
	}

	// create the graph
	arma::SpMat<double> edgeWeights(adjacencyMatrix.size(),adjacencyMatrix.size());
	MetaGraph meta(adjacencyMatrix,weightList,labelList,edgeWeights);

	return meta;
}

/*
 * Returns the fraction of nodes with opinion "True" or "Up" in the given cluster.
 * Complexity: O(n)
 */
template<typename GraphType>
inline double Clustering<GraphType>::getUpFractionOfCluster(const long int& clusterNumber) const {
	// sanity check of input
	if (getNumberOfNodes(clusterNumber) <= 0) {
		printf(" # Empty cluster in getUpFractionOfCluster. Exit \n");
		exit(1);
	}

	// prepare counter
	double counterUp = 0.;

	// loop over cluster nodes
	for (long int & node : getNodes(clusterNumber)){
		if (learningNetwork.getOpinion(node) == true) {
			++counterUp;
		}
	}

	// return the fraction
	return counterUp/getNumberOfNodes(clusterNumber);

}

template<typename GraphType>
inline void Clustering<GraphType>::writeNodesAndClusters(const std::string& filename, const long int& timeStep) const {
	// file object
	std::ofstream file;

		try {
			// Open file
			file.open(filename + "_" + std::to_string(timeStep), std::ofstream::out);

			// write heade
			file << "# Node ClusterNumber" << std::endl;

			// loop over nodes
			for (auto it = clustering.left.begin(); it != clustering.left.end(); ++it) {
				// write node and cluster number
				file << it->first << " " << it->second << std::endl;
			}

			// close file
			file.close();

		} catch (...) {
			std::cerr << "Could not open or write to file " << filename << "." << std::endl;
		}
}

/*
 * Returns a random non empty cluster.
 * Complexity: O(n)
 */
template<typename GraphType>
inline long int Clustering<GraphType>::getRandomCluster() {
	// get the list
	std::vector<long int> listOfKeys = getListOfNonEmptyClusters();

	// get a random element from it
	return listOfKeys.at(RandomNumbers::getRandomLongInt(0, listOfKeys.size()));

}

/*
 * Removes a node from a cluster.
 */
template<typename GraphType>
inline void Clustering<GraphType>::removeNode(const long int& nodeNumber, const long int& clusterNumber) {
	clustering.erase(relation(nodeNumber, clusterNumber));
}

/*
 * Inserts a node into a cluster.
 */
template<typename GraphType>
inline void Clustering<GraphType>::insertNode(const long int& nodeNumber, const long int& clusterNumber) {
	clustering.insert(relation(nodeNumber, clusterNumber));
}

/*
 * Checks if a node is part of a cluster.
 */
template<typename GraphType>
inline bool Clustering<GraphType>::checkNodeExistence(const long int& node, const long int& clusterNumber) const {
	auto it = clustering.left.find(node);
	if (it->second == clusterNumber) {
		return true;
	} else {
		return false;
	}
}

/*
 * Writes a dot file with color information.
 * The maximum number of clusters that can be visualized is 729.
 */
template<typename GraphType>
inline void Clustering<GraphType>::printClustering(const std::string& filename) const {
	// store number of clusters
	unsigned long int clusteringSize = getNumberOfClusters();

	// http://stackoverflow.com/questions/309149/generate-distinctly-different-rgb-colors-in-graphs
	const std::vector<std::string> colorList = { "#9BC4E5", "#04640D", "#FEFB0A", "#FB5514", "#E115C0", "#00587F", "#0BC582", "#FEB8C8", "#9E8317", "#01190F", "#847D81", "#58018B",
			"#B70639", "#703B01", "#F7F1DF", "#118B8A", "#4AFEFA", "#FCB164", "#796EE6", "#000D2C", "#53495F", "#F95475", "#61FC03", "#5D9608", "#DE98FD", "#98A088", "#4F584E",
			"#248AD0", "#5C5300", "#9F6551", "#BCFEC6", "#932C70", "#2B1B04", "#B5AFC4", "#D4C67A", "#AE7AA1", "#C2A393", "#0232FD", "#6A3A35", "#BA6801", "#168E5C", "#16C0D0",
			"#C62100", "#014347", "#233809", "#42083B", "#82785D", "#023087", "#B7DAD2", "#196956", "#8C41BB", "#ECEDFE", "#2B2D32", "#94C661", "#F8907D", "#895E6B", "#788E95",
			"#FB6AB8", "#576094", "#DB1474", "#8489AE", "#860E04", "#FBC206", "#6EAB9B", "#F2CDFE", "#645341", "#760035", "#647A41", "#496E76", "#E3F894", "#F9D7CD", "#876128",
			"#A1A711", "#01FB92", "#FD0F31", "#BE8485", "#C660FB", "#120104", "#D48958", "#05AEE8", "#C3C1BE", "#9F98F8", "#1167D9", "#D19012", "#B7D802", "#826392", "#5E7A6A",
			"#B29869", "#1D0051", "#8BE7FC", "#76E0C1", "#BACFA7", "#11BA09", "#462C36", "#65407D", "#491803", "#F5D2A8", "#03422C", "#72A46E", "#128EAC", "#47545E", "#B95C69",
			"#A14D12", "#C4C8FA", "#372A55", "#3F3610", "#D3A2C6", "#719FFA", "#0D841A", "#4C5B32", "#9DB3B7", "#B14F8F", "#747103", "#9F816D", "#D26A5B", "#8B934B", "#F98500",
			"#002935", "#D7F3FE", "#FCB899", "#1C0720", "#6B5F61", "#F98A9D", "#9B72C2", "#A6919D", "#2C3729", "#D7C70B", "#9F9992", "#EFFBD0", "#FDE2F1", "#923A52", "#5140A7",
			"#BC14FD", "#6D706C", "#0007C4", "#C6A62F", "#000C14", "#904431", "#600013", "#1C1B08", "#693955", "#5E7C99", "#6C6E82", "#D0AFB3", "#493B36", "#AC93CE", "#C4BA9C",
			"#09C4B8", "#69A5B8", "#374869", "#F868ED", "#E70850", "#C04841", "#C36333", "#700366", "#8A7A93", "#52351D", "#B503A2", "#D17190", "#A0F086", "#7B41FC", "#0EA64F",
			"#017499", "#08A882", "#7300CD", "#A9B074", "#4E6301", "#AB7E41", "#547FF4", "#134DAC", "#FDEC87", "#056164", "#FE12A0", "#C264BA", "#939DAD", "#0BCDFA", "#277442",
			"#1BDE4A", "#826958", "#977678", "#BAFCE8", "#7D8475", "#8CCF95", "#726638", "#FEA8EB", "#EAFEF0", "#6B9279", "#C2FE4B", "#304041", "#1EA6A7", "#022403", "#062A47",
			"#054B17", "#F4C673", "#02FEC7", "#9DBAA8", "#775551", "#835536", "#565BCC", "#80D7D2", "#7AD607", "#696F54", "#87089A", "#664B19", "#242235", "#7DB00D", "#BFC7D6",
			"#D5A97E", "#433F31", "#311A18", "#FDB2AB", "#D586C9", "#7A5FB1", "#32544A", "#EFE3AF", "#859D96", "#2B8570", "#8B282D", "#E16A07", "#4B0125", "#021083", "#114558",
			"#F707F9", "#C78571", "#7FB9BC", "#FC7F4B", "#8D4A92", "#6B3119", "#884F74", "#994E4F", "#9DA9D3", "#867B40", "#CED5C4", "#1CA2FE", "#D9C5B4", "#FEAA00", "#507B01",
			"#A7D0DB", "#53858D", "#588F4A", "#FBEEEC", "#FC93C1", "#D7CCD4", "#3E4A02", "#C8B1E2", "#7A8B62", "#9A5AE2", "#896C04", "#B1121C", "#402D7D", "#858701", "#D498A6",
			"#B484EF", "#5C474C", "#067881", "#C0F9FC", "#726075", "#8D3101", "#6C93B2", "#A26B3F", "#AA6582", "#4F4C4F", "#5A563D", "#E83005", "#32492D", "#FC7272", "#B9C457",
			"#552A5B", "#B50464", "#616E79", "#DCE2E4", "#CF8028", "#0AE2F0", "#4F1E24", "#FD5E46", "#4B694E", "#C5DEFC", "#5DC262", "#022D26", "#7776B8", "#FD9F66", "#B049B8",
			"#988F73", "#BE385A", "#2B2126", "#54805A", "#141B55", "#67C09B", "#456989", "#DDC1D9", "#166175", "#C1E29C", "#A397B5", "#2E2922", "#ABDBBE", "#B4A6A8", "#A06B07",
			"#A99949", "#0A0618", "#B14E2E", "#60557D", "#D4A556", "#82A752", "#4A005B", "#3C404F", "#6E6657", "#7E8BD5", "#1275B8", "#D79E92", "#230735", "#661849", "#7A8391",
			"#FE0F7B", "#B0B6A9", "#629591", "#D05591", "#97B68A", "#97939A", "#035E38", "#53E19E", "#DFD7F9", "#02436C", "#525A72", "#059A0E", "#3E736C", "#AC8E87", "#D10C92",
			"#B9906E", "#66BDFD", "#C0ABFD", "#0734BC", "#341224", "#8AAAC1", "#0E0B03", "#414522", "#6A2F3E", "#2D9A8A", "#4568FD", "#FDE6D2", "#FEE007", "#9A003C", "#AC8190",
			"#DCDD58", "#B7903D", "#1F2927", "#9B02E6", "#827A71", "#878B8A", "#8F724F", "#AC4B70", "#37233B", "#385559", "#F347C7", "#9DB4FE", "#D57179", "#DE505A", "#37F7DD",
			"#503500", "#1C2401", "#DD0323", "#00A4BA", "#955602", "#FA5B94", "#AA766C", "#B8E067", "#6A807E", "#4D2E27", "#73BED7", "#D7BC8A", "#614539", "#526861", "#716D96",
			"#829A17", "#210109", "#436C2D", "#784955", "#987BAB", "#8F0152", "#0452FA", "#B67757", "#A1659F", "#D4F8D8", "#48416F", "#DEBAAF", "#A5A9AA", "#8C6B83", "#403740",
			"#70872B", "#D9744D", "#151E2C", "#5C5E5E", "#B47C02", "#F4CBD0", "#E49D7D", "#DD9954", "#B0A18B", "#2B5308", "#EDFD64", "#9D72FC", "#2A3351", "#68496C", "#C94801",
			"#EED05E", "#826F6D", "#E0D6BB", "#5B6DB4", "#662F98", "#0C97CA", "#C1CA89", "#755A03", "#DFA619", "#CD70A8", "#BBC9C7", "#F6BCE3", "#A16462", "#01D0AA", "#87C6B3",
			"#E7B2FA", "#D85379", "#643AD5", "#D18AAE", "#13FD5E", "#B3E3FD", "#C977DB", "#C1A7BB", "#9286CB", "#A19B6A", "#8FFED7", "#6B1F17", "#DF503A", "#10DDD7", "#9A8457",
			"#60672F", "#7D327D", "#DD8782", "#59AC42", "#82FDB8", "#FC8AE7", "#909F6F", "#B691AE", "#B811CD", "#BCB24E", "#CB4BD9", "#2B2304", "#AA9501", "#5D5096", "#403221",
			"#F9FAB4", "#3990FC", "#70DE7F", "#95857F", "#84A385", "#50996F", "#797B53", "#7B6142", "#81D5FE", "#9CC428", "#0B0438", "#3E2005", "#4B7C91", "#523854", "#005EA9",
			"#F0C7AD", "#ACB799", "#FAC08E", "#502239", "#BFAB6A", "#2B3C48", "#0EB5D8", "#8A5647", "#49AF74", "#067AE9", "#F19509", "#554628", "#4426A4", "#7352C9", "#3F4287",
			"#8B655E", "#B480BF", "#9BA74C", "#5F514C", "#CC9BDC", "#BA7942", "#1C4138", "#3C3C3A", "#29B09C", "#02923F", "#701D2B", "#36577C", "#3F00EA", "#3D959E", "#440601",
			"#8AEFF3", "#6D442A", "#BEB1A8", "#A11C02", "#8383FE", "#A73839", "#DBDE8A", "#0283B3", "#888597", "#32592E", "#F5FDFA", "#01191B", "#AC707A", "#B6BD03", "#027B59",
			"#7B4F08", "#957737", "#83727D", "#035543", "#6F7E64", "#C39999", "#52847A", "#925AAC", "#77CEDA", "#516369", "#E0D7D0", "#FCDD97", "#555424", "#96E6B6", "#85BB74",
			"#5E2074", "#BD5E48", "#9BEE53", "#1A351E", "#3148CD", "#71575F", "#69A6D0", "#391A62", "#E79EA0", "#1C0F03", "#1B1636", "#D20C39", "#765396", "#7402FE", "#447F3E",
			"#CFD0A8", "#3A2600", "#685AFC", "#A4B3C6", "#534302", "#9AA097", "#FD5154", "#9B0085", "#403956", "#80A1A7", "#6E7A9A", "#605E6A", "#86F0E2", "#5A2B01", "#7E3D43",
			"#ED823B", "#32331B", "#424837", "#40755E", "#524F48", "#B75807", "#B40080", "#5B8CA1", "#FDCFE5", "#CCFEAC", "#755847", "#CAB296", "#C0D6E3", "#2D7100", "#D5E4DE",
			"#362823", "#69C63C", "#AC3801", "#163132", "#4750A6", "#61B8B2", "#FCC4B5", "#DEBA2E", "#FE0449", "#737930", "#8470AB", "#687D87", "#D7B760", "#6AAB86", "#8398B8",
			"#B7B6BF", "#92C4A1", "#B6084F", "#853B5E", "#D0BCBA", "#92826D", "#C6DDC6", "#BE5F5A", "#280021", "#435743", "#874514", "#63675A", "#E97963", "#8F9C9E", "#985262",
			"#909081", "#023508", "#DDADBF", "#D78493", "#363900", "#5B0120", "#603C47", "#C3955D", "#AC61CB", "#FD7BA7", "#716C74", "#8D895B", "#071001", "#82B4F2", "#B6BBD8",
			"#71887A", "#8B9FE3", "#997158", "#65A6AB", "#2E3067", "#321301", "#FEECCB", "#3B5E72", "#C8FE85", "#A1DCDF", "#CB49A6", "#B1C5E4", "#3E5EB0", "#88AEA7", "#04504C",
			"#975232", "#6786B9", "#068797", "#9A98C4", "#A1C3C2", "#1C3967", "#DBEA07", "#789658", "#E7E7C6", "#A6C886", "#957F89", "#752E62", "#171518", "#A75648", "#01D26F",
			"#0F535D", "#047E76", "#C54754", "#5D6E88", "#AB9483", "#803B99", "#FA9C48", "#4A8A22", "#654A5C", "#965F86", "#9D0CBB", "#A0E8A0", "#D3DBFA", "#FD908F", "#AEAB85",
			"#A13B89", "#F1B350", "#066898", "#948A42", "#C8BEDE", "#19252C", "#7046AA", "#E1EEFC", "#3E6557", "#CD3F26", "#2B1925", "#DDAD94", "#C0B109", "#37DFFE", "#039676",
			"#907468", "#9E86A5", "#3A1B49", "#BEE5B7", "#C29501", "#9E3645", "#DC580A", "#645631", "#444B4B", "#FD1A63", "#DDE5AE", "#887800", "#36006F", "#3A6260", "#784637",
			"#FEA0B7", "#A3E0D2", "#6D6316", "#5F7172", "#B99EC7", "#777A7E", "#E0FEFD", "#E16DC5", "#01344B", "#F8F8FC", "#9F9FB5", "#182617", "#FE3D21", "#7D0017", "#822F21",
			"#EFD9DC", "#6E68C4", "#35473E", "#007523", "#767667", "#A6825D", "#83DC5F", "#227285", "#A95E34", "#526172", "#979730", "#756F6D", "#716259", "#E8B2B5", "#B6C9BB",
			"#9078DA", "#4F326E", "#B2387B", "#888C6F", "#314B5F", "#E5B678", "#38A3C6", "#586148", "#5C515B", "#CDCCE1", "#C8977F" };

	if (colorList.size() < clusteringSize) {
		printf("Not enough color values. Exit.\n");
		exit(1);
	}

	// ==== write to file ====
	std::ofstream file;
	try {
		// Open file
		file.open(filename, std::ofstream::out);

		file << "Digraph {" << std::endl;
		file << "node [style=filled]\n"
				//<<	"node [label=\"\"]\n"
				 << "node [shape=\"circle\"]\n" << "overlap=false\n" << "edge [color=\"#8598CC\"]\n"
				<< "edge [arrowsize=0.6]\n"
				<< "splines=true;\n"
				<< "sep=\"+5,5\";\n"
				<< "node  [fontcolor=\"#ffffff\"] \n"
				<< std::endl;

		//===== write data to file
		// counter for cluster index
		long int count = 0;

		// loop over clusters
		for (long int clusterNumber : getListOfNonEmptyClusters()) {

			// get the nodes
			std::vector<long int> nodes = getNodes(clusterNumber);

			// write cluster to file
			//file << "subgraph cluster_" + std::to_string(count) + "{" << std::endl;
			//file << "style=invisible;"<< std::endl;
			for (long int node : nodes) {
			// write label
			std::string label;
			if (learningNetwork.getOpinion(node)==true) {
				label="1";
			} else {
				label="0";
			}
			//file <<std::to_string(node) << "[label=\""<< std::to_string(node)<<" \n "<<label<<"\", color=\"" + colorList.at(clusterNumber) + "\"]" << std::endl;
			file <<std::to_string(node) << "[label=\""<<label<<"\"]" << std::endl;
			file <<std::to_string(node) << "[ weight=\"" + std::to_string(clusterNumber) + "\"]" << std::endl;

			}
			//file << "}" << std::endl;

			// write edges to file
			for (long int node : nodes) {
				// write edges
				for (long int outNode : learningNetwork.getAdjacentOutNodes(node)) {
					file << std::to_string(node) << "->" << std::to_string(outNode) << std::endl;
					// check if the node is in the same cluster
					if (clustering.left.find(outNode)->second == clusterNumber) {
						// add same color for edge
						//file <<" [color=\""<< colorList.at(clusterNumber)<<" \"];\n"<< std::endl;
						file <<" [weight=5 ];\n"<< std::endl;
					} else {
						// require minimum length
						//file <<" [minlen=25];\n"<< std::endl;
						file <<" [weight=1 ];\n"<< std::endl;

					}

				}

			}

			++count;

		}

		// write closing
		file << "}" << std::endl;

		// close file
		file.close();
		std::cerr << "# Finished writing graph to " << filename << std::endl;

	} catch (...) {
		std::cerr << "Could not open or write to file " << filename << " or index out of bounds." << std::endl;
	}

}

/*
 * Performs the clustering.
 * Complexity O(n^2).
 */
template<typename GraphType>
inline void Clustering<GraphType>::doFullClustering() {

	// do a modularity optimization step
	doModularityOptimization(3 * learningNetwork.getSize());

	// do the entropy optimization
	doDissimilarityOptimization(3*getNumberOfClusters());


	/*printf("# mod %ld entro %ld \n",numberModularityClusters,numberEntropyClusters);

		return 1.0-((double) numberModularityClusters-(double) numberEntropyClusters)/( (double) numberModularityClusters);*/

}

/*
 * Calculation of cluster information entropy.
 * Empty clusters and clusters with only one element have zero entropy.
 * Complexity: O(n)
 */
template<typename GraphType>
inline double Clustering<GraphType>::getClusterEntropy(const long int& cluster) const {
	// start with zero entropy
	double entropy = 0.0;

	if (getNumberOfNodes(cluster) > 0) {
		// get the node numbers of the cluster
		std::vector<long int> nodes = getNodes(cluster);

		// loop over pairs
		for (unsigned long int i = 0; i < nodes.size(); ++i) {
			for (unsigned long int l = i + 1; l < nodes.size(); ++l) {

				double sim = learningNetwork.getSimilarity(nodes.at(i), nodes.at(l));

				/*if (sim != 0.0 && sim != 1.0) {
					entropy -= sim * log(sim) + (1 - sim) * log(1 - sim);
				} */// else zero would be subtracted

				entropy+=1-sim;
			}
		}

		entropy/=getNumberOfNodes(cluster);
	}

	return entropy;
}

/*
 * Method from original paper. Moves nodes to unconnected cluster -> unwanted!
 * Complexity: O(maxIterations)
 */
template<typename GraphType>
inline void Clustering<GraphType>::doEntropyOptimization(const long int& maxIterations) {

	// storage of entropies to avoid unneccessary calculations
	std::map<long int,double> entropies;

	// set up entropies
	for (long int clusterId : getListOfNonEmptyClusters()){
		entropies[clusterId] = getClusterEntropy(clusterId);
	}

	// counter for acceptance rate
	double acceptedMoves = 0;

	for (long int iterations = 0; iterations < maxIterations; ++iterations) {

		// get the first cluster at random
		long int sourceCluster = getRandomCluster();

		// try to find another unconnected cluster
		bool successfulSearch = false;
		long int destinationCluster = getUnconnectedRandomCluster(sourceCluster, successfulSearch);

		if (successfulSearch == true) {

			// === try to decrease the entropy ===

			// save the node
			long int movingNode = getRandomNode(sourceCluster);

			// add nodes to new cluster
			removeNode(movingNode, sourceCluster);
			insertNode(movingNode, destinationCluster);

			// calculate entropy of new partition
			double entropyA = getClusterEntropy(sourceCluster);
			double entropyB = getClusterEntropy(destinationCluster);

			// check if move decreased the entropy
			if (entropyA + entropyB < entropies.at(sourceCluster) + entropies.at(destinationCluster)) {

				// accept move and store the new entropies
				entropies.at(sourceCluster) = entropyA;
				entropies.at(destinationCluster) = entropyB;

				// counter up
				++acceptedMoves;

			} else {

				// discard the move
				removeNode(movingNode, destinationCluster);
				insertNode(movingNode, sourceCluster);

			}

		} // -- if statement success

	} // -- iterations loop

}

/*
 * Calculation of modularity for whole clustering.
 * Complexity O(n^2)
 */
template<typename GraphType>
inline double Clustering<GraphType>::calculateTotalModularity() {

	// set up modularity matrix
	modularityMatrix = arma::mat(clustering.left.size(), clustering.left.size(), arma::fill::zeros);

	for (auto it = clustering.left.begin(); it != clustering.left.end(); ++it) {
		// get the cluster number of current node
		long int clusterNumberNode = it->second;

		for (long int neighbour : learningNetwork.getAdjacentAllNodes(it->first)) {
			// get the cluster number of the neighbour
			long int clusterNumberNeighbour = clustering.left.find(neighbour)->second;

			// increment the count in the modularity matrix
			modularityMatrix(clusterNumberNode, clusterNumberNeighbour) += 1;
		}
	}

	// normalize it (the factor two arises because the edges are counted twice)
	modularityMatrix = modularityMatrix / (2.0 * learningNetwork.getNumberOfEdges()); // XXX factor 2 ?

	// return the modularity
	return arma::trace(modularityMatrix) - arma::accu(modularityMatrix * modularityMatrix);
}

/*
 * Performs a modularity optimization.
 * Complexity O(maxIterations)
 */
template<typename GraphType>
inline void Clustering<GraphType>::doModularityOptimization(const long int& maxIterations) {

	for (long int iteration = 0; iteration < maxIterations; ++iteration) {
		// pick a random cluster
		long int sourceCluster = getRandomCluster();

		// store the nodes
		std::vector<long int> nodes = getNodes(sourceCluster);

		for (long int node : nodes) {
			// get all connected clusters
			std::vector<long int> connectedClusters = getConnectedClusters(node);

			// remove duplicate entries
			std::sort(connectedClusters.begin(), connectedClusters.end());
			connectedClusters.erase(std::unique(connectedClusters.begin(), connectedClusters.end()), connectedClusters.end());

			// variable for maximum  modularity gain and corresponding cluster
			double maxModularityGain = 0;
			long int bestClusterID = -1;

			// loop over destination clusters
			for (long int destinationCluster : connectedClusters) {

				if (destinationCluster != sourceCluster) {

					// calculate the modularity
					double modularity = getModularity(sourceCluster) + getModularity(destinationCluster);

					// move the node
					removeNode(node, sourceCluster);
					insertNode(node, destinationCluster);

					// calculate the modularity again
					double newModularity = 0;
					if (getNumberOfNodes(sourceCluster) == 0) {
						// implies mod = 0 for empty cluster
						newModularity = getModularity(destinationCluster);
					} else {
						newModularity = getModularity(sourceCluster) + getModularity(destinationCluster);
					}

					if (newModularity - modularity > maxModularityGain) {
						// store the results
						maxModularityGain = newModularity - modularity;
						bestClusterID = destinationCluster;
					}

					// undo the move
					removeNode(node, destinationCluster);
					insertNode(node, sourceCluster);
				}

			}

			// check if there is a improving move
			if (bestClusterID != -1) {
				// make the move with the best result
				removeNode(node, sourceCluster);
				insertNode(node, bestClusterID);
			}
		}

	}
}

/*
 * Returns a list with unique numbers of nodes and neighbors of the given cluster.
 */
template<typename GraphType>
inline std::vector<long int> Clustering<GraphType>::getNodesAndNeighbors(const long int& clusterNumber) const {
	// declare empty vector
	std::vector<long int> list;

	if (getNumberOfNodes(clusterNumber) > 0) {
		// add all the nodes
		std::vector<long int> nodes = getNodes(clusterNumber);
		list.insert(list.end(), std::make_move_iterator(nodes.begin()), std::make_move_iterator(nodes.end()));

		// list of neighbors
		std::vector<long int> neighbors;

		// iterate over nodes and get their neighbors
		for (long int node : nodes) {
			neighbors = learningNetwork.getAdjacentAllNodes(node);
			list.insert(list.end(), std::make_move_iterator(neighbors.begin()), std::make_move_iterator(neighbors.end()));
		}

		// remove multiply occurrences of nodes
		std::sort(list.begin(), list.end());
		list.erase(std::unique(list.begin(), list.end()), list.end());
	}

	return list;
}

/*
 * Returns a unconnected cluster.
 * It is possible that there is no unconnected cluster. Therefore the flag variable for success has to be passed as a reference
 * and CHECKED BY THE USER.
 * In case of failure the method will return then id of the source Cluster, which is valid data, but might lead to LOCIGAL ERRORS.
 */
template<typename GraphType>
inline long int Clustering<GraphType>::getUnconnectedRandomCluster(const long int& sourceCluster, bool& flag) {

	// flag for finding it
	bool foundUnconnected = false;

	// flag to prevent infinite loop
	bool triedManyPossibilities = false;

	// store the already visited clusters
	std::set<long int> visitedClusters;
	visitedClusters.insert(sourceCluster);

	do {

		// get another cluster
		long int destinationCluster = getRandomCluster();

		// check if has not already been visited
		if (visitedClusters.count(destinationCluster) == 0) {
			// add it to the visited cluster set
			visitedClusters.insert(destinationCluster);

			// set the flag here to save instructions
			foundUnconnected = true;

			if (getNumberOfNodes(destinationCluster) > 0) {
				// get the nodes and neighbors
				std::vector<long int> nodesAndNeighbors = getNodesAndNeighbors(destinationCluster);

				// loop through nodes and neighbors
				for (long int node : nodesAndNeighbors) {
					if (checkNodeExistence(node, sourceCluster) == true) {
						foundUnconnected = false;
						break;
					}
				}

			} else {
				// it is empty and thus unconnected
				flag = true;
				return destinationCluster;

			}

		}

		if (visitedClusters.size() >= ceil(0.5 * getNumberOfClusters())) {
			triedManyPossibilities = true;
			flag = false;
			return sourceCluster;

		}

		if (foundUnconnected == true) {
			flag = true;
			return destinationCluster;

		}

	} while (!foundUnconnected && !triedManyPossibilities);
}

/*
 * Returns the modularity of the given cluster.
 * Complexity: O(n)
 */
template<typename GraphType>
inline double Clustering<GraphType>::getModularity(const long int& clusterNumber) const {
	// get list of non empty clusters
	std::vector<long int> clusterList = getListOfNonEmptyClusters();

	// it is already sorted so you can use it as a map:
	// the index of the array will be used as the new index for the cluster

	// set up vector for row in modularity matrix
	arma::vec row(clusterList.size(), arma::fill::zeros);

	// loop over nodes
	for (long int node : getNodes(clusterNumber)) {
		// get the clusters
		std::vector<long int> connectedClusters = getConnectedClusters(node);

		// loop over clusters
		for (long int cluster : connectedClusters) {

			// find the new index of the cluster
			auto it = std::find(clusterList.begin(),clusterList.end(),cluster);
			long int newIndex = std::distance(clusterList.begin(),it);

			// increment counter for connection
			row(newIndex) += 1;
		}
	}

	// normalize it
	row = row / (2.0 * learningNetwork.getNumberOfEdges());

	// find new index of clusterNumber
	long int newClusterNumber = std::distance(clusterList.begin(),std::find(clusterList.begin(),clusterList.end(),clusterNumber));

	// calculate and return the modularity
	return row(newClusterNumber) - pow(arma::accu(row), 2.0);

}

/*
 * Returns a list of the cluster numbers of the clusters the given node is connected to.
 * The same clusters may appear several times!
 * Complexity:
 */
template<typename GraphType>
inline std::vector<long int> Clustering<GraphType>::getConnectedClusters(const long int& nodeNumber) const {
	// get the neighbors
	std::vector<long int> neighbours = learningNetwork.getAdjacentAllNodes(nodeNumber);

	// empty list of cluster numbers
	std::vector<long int> clusterNumbers;

	// get the cluster numbers of the nodes
	for (long int node : neighbours) {
		auto it = clustering.left.find(node);
		clusterNumbers.push_back(it->second);
	}

	// return it
	return clusterNumbers;
}

/*
 * Performs a splitting of a cluster into two clusters by optimization of the entropy.
 */
template<typename GraphType>
inline void Clustering<GraphType>::doEntropySplitting(const long int& clusterNumber) {
	// =====create a new cluster with an unused id
	long int newCluster = -1;
	bool foundId = false;
	do {
		newCluster = RandomNumbers::getRandomLongInt(0, learningNetwork.getSize());
		if (clustering.right.count(newCluster) == 0) {
			foundId = true;
		}
	} while (!foundId);

	for (long int movingNode : getNodes(clusterNumber)) {
		// calculate old entropies
		double entropySourceOld = getClusterEntropy(clusterNumber);
		double entropyDestinationOld = getClusterEntropy(newCluster);

		// move node to new cluster
		removeNode(movingNode, clusterNumber);
		insertNode(movingNode, newCluster);

		// calculate entropy of new partition
		double entropySourceNew = getClusterEntropy(clusterNumber);
		double entropyDestinationNew = getClusterEntropy(newCluster);



		// check if move decreased the entropy
		if (entropySourceNew + entropyDestinationNew > entropySourceOld +entropyDestinationOld ) {
			// discard the move
			removeNode(movingNode, newCluster);
			insertNode(movingNode, clusterNumber);
		} //else  accept move, do nothing
	}

	// back again
	for (long int movingNode : getNodes(newCluster)) {
		// calculate old entropies
		double entropySourceOld = getClusterEntropy(newCluster);
		double entropyDestinationOld = getClusterEntropy(clusterNumber);

		// move node to new cluster
		removeNode(movingNode, newCluster);
		insertNode(movingNode, clusterNumber);

		// calculate entropy of new partition
		double entropySourceNew = getClusterEntropy(newCluster);
		double entropyDestinationNew = getClusterEntropy(clusterNumber);

		// check if move decreased the entropy
		if (entropySourceNew + entropyDestinationNew > entropySourceOld + entropyDestinationOld) {
			// discard the move
			removeNode(movingNode, clusterNumber);
			insertNode(movingNode, newCluster);
		} //else  accept move, do nothing
	}

}

/*
 * Writes a file for use with scikit learn.
 */
template<typename GraphType>
inline void Clustering<GraphType>::writeClusteringFileScikit(const std::string& filename) const {
	std::ofstream file;
		try {
			// Open file
			file.open(filename, std::ofstream::out);

			// fill it
			for (auto it = clustering.left.begin(); it != clustering.left.end(); ++it) {
					file<<it->second<<std::endl;
			}

			// close file
			file.close();
			std::cerr << "Finished writing clustering to " << filename << std::endl;

		} catch (...) {
			std::cerr << "Could not open or write to file " << filename << " or index out of bounds." << std::endl;
		}

}

/*
 * Returns the total entropy of the clustering.
 * Complexity: O(n)
 */
template<typename GraphType>
inline double Clustering<GraphType>::calculateTotalEntropy() const {

	double entro = 0;
	for (long int cluster : getListOfNonEmptyClusters()) {
		entro += getClusterEntropy(cluster);
	}

	return entro;

}

template<typename GraphType>
inline void Clustering<GraphType>::doDissimilarityOptimization(const long int& maxIterations) {
	// storage of dissimilarities to avoid unneccessary calculations
	std::map<long int, double> dissimilarities;

	// set up dissimilarities
	for (long int clusterId : getListOfNonEmptyClusters()) {
		dissimilarities[clusterId] = getClusterEntropy(clusterId);
	}

	for (long int iterations = 0; iterations < maxIterations; ++iterations) {
		// get a random cluster
		long int sourceCluster = getRandomCluster();

			for (long int movingNode : getNodes(sourceCluster)) {

				// list of connected cluster for node
				std::vector<long int> clusters  =  getConnectedClusters(movingNode);

				// get size of clusters list
				unsigned long int size = clusters.size();

				if (size > 0) {
					// get random cluster (more connected have higher chance)
					long int destinationCluster = clusters.at(RandomNumbers::getRandomLongInt(0,size));

					// add nodes to new cluster
					removeNode(movingNode, sourceCluster);
					insertNode(movingNode, destinationCluster);

					// calculate entropy of new partition
					double dissimilarityA = getClusterDissimilarity(sourceCluster);
					double dissimilarityB = getClusterDissimilarity(destinationCluster);

					// check if move decreased the entropy
					if (dissimilarityA + dissimilarityB <= dissimilarities.at(sourceCluster) + dissimilarities.at(destinationCluster)) {

						// accept move and store the new dissimilarities
						dissimilarities.at(sourceCluster) = dissimilarityA;
						dissimilarities.at(destinationCluster) = dissimilarityB;

					} else {

						// discard the move
						removeNode(movingNode, destinationCluster);
						insertNode(movingNode, sourceCluster);

					}
				}// if size of clusters

			}// for loop nodes



	} // -- iterations loop


}

/*
 * Returns the order parameter.
 * Argument: number of clusters by modularity optimization at the start.
 * Complexity: O(n)
 */
template<typename GraphType>
inline double Clustering<GraphType>::getOrderParameter(const long int& numberClustersStart) const {
	return 1.0-((double) numberClustersStart-(double) getNumberOfClusters()-1.0)/((double) numberClustersStart);
}

template<typename GraphType>
inline double Clustering<GraphType>::getClusterDissimilarity(const long int& cluster) const {
	// start with zero entropy
	double dissimilarity = 0.0;

	if (getNumberOfNodes(cluster) > 0) {
		// get the node numbers of the cluster
		std::vector<long int> nodes = getNodes(cluster);

		// loop over pairs
		for (unsigned long int i = 0; i < nodes.size(); ++i) {
			for (unsigned long int l = i + 1; l < nodes.size(); ++l) {

				double sim = learningNetwork.getSimilarity(nodes.at(i), nodes.at(l));

				dissimilarity += 1 - sim;
			}
		}

		dissimilarity /= getNumberOfNodes(cluster);
	}

	return dissimilarity;
}

/*
 * Updates the clustering by doing the dissimilarity reduction only.
 * Complexity: O(n^2)
 */
template<typename GraphType>
inline void Clustering<GraphType>::updateClustering() {

	// use modularity based clustering
	clustering=clusteringModularity;

	// do the attribute based optimization
	doDissimilarityOptimization(3*getNumberOfClusters());
}

/*
 * Returns a histogram of the size distribution of the clusters.
 */
template<typename GraphType>
inline Histogram Clustering<GraphType>::getSizeDistribution(const double& binWidth) const {
	// set up storage for cluster sizes
	std::vector<double> rawSizes;
	rawSizes.reserve(getNumberOfClusters());

	for (long int cluster : getListOfNonEmptyClusters()){
		rawSizes.push_back(getNumberOfNodes(cluster));
	}

	// create Histogram
	Histogram hist(rawSizes,binWidth);

	return hist;
}

/*
 * Returns a list of cluster numbers that correspond to non empty clusters.
 * BONUS: Already sorted.
 * Complexity: O(n)
 */
template<typename GraphType>
inline std::vector<long int> Clustering<GraphType>::getListOfNonEmptyClusters() const {
	// set up a list
	std::vector<long int> listOfKeys;
	listOfKeys.reserve(clustering.right.size());

	// fill it
	for (auto it = clustering.right.begin(); it != clustering.right.end(); ++it) {
		listOfKeys.push_back(it->first);
	}

	// get list of unique numbers
	std::sort(listOfKeys.begin(), listOfKeys.end());
	listOfKeys.erase(std::unique(listOfKeys.begin(), listOfKeys.end()), listOfKeys.end());

	// return the list
	return listOfKeys;
}

/*
 * Creates an graph object of all nodes in the given cluster and edges within the cluster. Edges to nodes in other clusters are ignored.
 */
template<typename GraphType>
inline MetaGraph Clustering<GraphType>::createGraphOfClusterNodes(const long int& clusterNumber) const {
	// get all the nodes in the cluster
	std::vector<long int> clusterNodes = getNodes(clusterNumber);

	// sort for later speed up of search
	std::sort(clusterNodes.begin(), clusterNodes.end());

	// prepare the empty adjacency matrix
	arma::SpMat<short> adjacencyMatrix(clusterNodes.size(), clusterNodes.size());

	// ==== For compressing the adjacency matrix the node numbers need to be mapped onto 0,...,Size-1 ====
	// create map object
	std::map<long int,long int> newNodeNumbers;

	// fill it
	for (long int index  = 0; index < clusterNodes.size(); index++){
		newNodeNumbers[clusterNodes.at(index)] = index;
	}

	// ==== Iterate over the nodes to get the edges =========
	for (long int& node : clusterNodes){

		// get the neighbors
		std::vector<long int> neighbors = learningNetwork.getAdjacentInNodes(node);

		// check if they are inside the cluster
		for (long int& neighbor : neighbors){
			if(std::binary_search(clusterNodes.begin(), clusterNodes.end(), neighbor) == true) {

			    // Then the neighbor is inside the cluster and the edge can be added in the adjacency matrix
				adjacencyMatrix(newNodeNumbers.find(node)->second,newNodeNumbers.find(neighbor)->second) = 1;

			}

		}


	}

	// ==== create the graph object ====
	std::vector<double> dummyWeights(0,clusterNodes.size());
	MetaGraph clusterGraph(adjacencyMatrix,dummyWeights);

	return clusterGraph;
}

#endif /* CLUSTERING_HPP_ */
