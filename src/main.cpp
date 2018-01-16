#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>
#include "vec.hpp"

struct GraphNode {
	size_t id;
	vec coords;
	std::vector<size_t> connected{};
	std::string toString() {
		std::stringstream ss;
		ss << "[" << id << "] {" << coords.x() << "," << coords.y() << "," << coords.z()
		   << "}";
		for (auto n : connected) ss << " " << n;
		ss << std::endl;
		return ss.str();
	}
};

// just a quick&dirty hash function so that unordered_pair is happy to store pairs of uint
struct pair_hash {
	inline std::size_t operator()(const std::pair<size_t, size_t>& v) const {
		return v.first * 1000000000u + v.second;
	}
};

size_t getNumberOfEdges(const std::vector<GraphNode>& graph) {
	std::unordered_set<std::pair<size_t, size_t>, pair_hash> visitedEdges;
	for (const auto& n : graph)
		for (size_t c : n.connected) visitedEdges.insert(std::minmax(n.id, c));
	return visitedEdges.size();
}

template <typename F>
std::vector<size_t> eulerOpt(const std::vector<GraphNode>& graph, const F& heuristic) {
	// a variation of hierholzer's algorithm with an added heuristic function

	size_t startNode = 0;

	const size_t totalEdgeNumber = getNumberOfEdges(graph);
	// first we check that the graph is (semi)eulerian
	std::vector<size_t> oddNodes;
	for (const auto& n : graph)
		if (n.connected.size() % 2 == 1) oddNodes.push_back(n.id);
	if (oddNodes.size() > 2) {
		std::cout << "This graph has " << oddNodes.size() << " odds. Aborting." << std::endl;
		for (auto on : oddNodes) {
			std::cout << " " << on;
		}
		return {{}};
	} else if (oddNodes.size() > 0)
		startNode = oddNodes[0];

	// then we start looking for subcicle
	std::vector<size_t> eulerPath;
	std::unordered_set<std::pair<size_t, size_t>, pair_hash> visitedEdges;
	eulerPath.push_back(startNode);
	size_t startNodeIndex = 0;









	while (visitedEdges.size() < totalEdgeNumber) {
		std::vector<size_t> subCircle;
		subCircle.push_back(startNode);
		do {
			// we find the next node that minimises the heuristic
			size_t nextNode;
			bool foundNext = false;
			double minHeuristic = std::numeric_limits<double>::max();
			size_t currentNode = subCircle.back();
			for (auto& nId : graph[currentNode].connected) {
				if (!visitedEdges.count(std::minmax(currentNode, nId))) {
					double eurVal = heuristic(graph, subCircle, nId);
					if (eurVal < minHeuristic) {
						nextNode = nId;
						minHeuristic = eurVal;
						foundNext = true;
					}
				}
			}
			std::cerr << "Chose " << nextNode << std::endl;
			assert(foundNext);
			visitedEdges.insert(std::minmax(currentNode, nextNode));
			subCircle.push_back(nextNode);
		} while (subCircle.back() != startNode);
		// subcircle formed
		// we first replace startNode in eulerGraph by this new subcircle
		eulerPath.insert(eulerPath.begin() + (long)startNodeIndex, subCircle.begin(),
		                 subCircle.end());
		// we then need to find the next unvisited edge from which to start
		// TODO: use the heuristic to find the node with the best unvisited edge... for now we
		// just take the first with any unvisited edge
		for (startNodeIndex = 0; startNodeIndex < eulerPath.size(); ++startNodeIndex) {
			bool hasUnvisitedEdge = false;
			startNode = eulerPath[startNodeIndex];
			for (auto& nId : graph[startNode].connected) {
				if (!visitedEdges.count(std::minmax(startNode, nId))) {
					hasUnvisitedEdge = true;
					break;
				}
			}
			if (hasUnvisitedEdge)
				break;
		}
		// now startNodeIndex points to the next startNode
	}
	return eulerPath;
}

std::vector<GraphNode> readGraph(std::string filePath) {
	std::vector<GraphNode> graph;
	std::ifstream infile(filePath);
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string rawInput;
		std::vector<std::string> numbers;
		while (std::getline(iss, rawInput, ' ')) {
			numbers.push_back(rawInput);
		}
		assert(numbers.size() > 3);
		GraphNode node;
		node.id = (size_t)std::stoi(numbers[0]);
		node.coords =
		    vec(std::stod(numbers[1]), std::stod(numbers[2]), std::stod(numbers[3]));
		for (size_t i = 4; i < numbers.size(); ++i) {
			if (!std::all_of(numbers[i].begin(), numbers[i].end(), isspace)) {
				node.connected.push_back((size_t)std::stoi(numbers[i]));
			}
		}
		graph.push_back(node);
	}
	return graph;
}

int main(int argc, char** argv) {
	if (argc <= 1) {
		std::cerr << "Usage: " << argv[0] << " graphFilePath" << std::endl;
		return 1;
	}
	auto graph = readGraph(argv[1]);
	auto dumbHeuristic = [](const auto&, const auto&, size_t) { return 0.0; };
	auto angleHeuristic = [](const auto& g, const auto& p, size_t n) {
		if (p.size() < 2) return 0.0;
		auto prevEdge = g[p.back()].coords - g[p[p.size() - 2]].coords;
		auto nextEdge = g[n].coords - g[p.back()].coords;
		double d = prevEdge.normalized().dot(nextEdge.normalized());
		return -d;
	};
	auto eulerPath = eulerOpt(graph, angleHeuristic);
	for (auto i : eulerPath) {
		std::cout << i << std::endl;
	}
	return 0;
}
