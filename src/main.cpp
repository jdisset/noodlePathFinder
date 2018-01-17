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

// finds the best unvisited edge. Returns -1 if all visited
template <class G, class V, class F>
std::tuple<size_t, double, bool> findBestNext(size_t prevNode, size_t currentNode,
                                              const G& gr, const V& visited,
                                              const F& heuristic) {
	int nextNode = -1;
	double minHeuristic = std::numeric_limits<double>::max();
	for (auto& nId : gr[currentNode].connected) {
		if (!visited.count(std::minmax(currentNode, nId))) {
			double eurVal = heuristic(gr, prevNode, currentNode, nId);
			if (eurVal < minHeuristic) {
				nextNode = nId;
				minHeuristic = eurVal;
			}
		}
	}
	bool hasUnvisitedEdges = (nextNode >= 0);
	return std::make_tuple(nextNode, minHeuristic, hasUnvisitedEdges);
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

	// then we start looking for subcircles
	std::vector<size_t> eulerPath;
	std::unordered_set<std::pair<size_t, size_t>, pair_hash> visitedEdges;
	eulerPath.push_back(startNode);
	size_t startNodeIndex = 0;

	while (visitedEdges.size() < totalEdgeNumber) {
		std::vector<size_t> subCircle;
		subCircle.push_back(startNode);
		do {
			// we need the previous node for the heuristic to work
			size_t prevNode = subCircle.front();
			if (subCircle.size() < 2) {
				if (startNodeIndex > 0) prevNode = eulerPath[startNodeIndex - 1];
			} else {
				prevNode = subCircle[subCircle.size() - 2];
			}
			size_t currentNode = subCircle.back();
			// we find the next node that minimises the heuristic
			auto bestNext = findBestNext(prevNode, currentNode, graph, visitedEdges, heuristic);
			size_t nextNode = std::get<0>(bestNext);
			assert(std::get<2>(bestNext));
			visitedEdges.insert(std::minmax(currentNode, nextNode));
			subCircle.push_back(nextNode);
		} while (subCircle.back() != startNode);

		// subcircle formed
		// we first replace startNode in eulerGraph by this new subcircle
		eulerPath.insert(eulerPath.begin() + (long)startNodeIndex, subCircle.begin(),
		                 subCircle.end());
		// we then need to find the next unvisited edge with the best possible next edge from
		// which to start
		size_t bestStartNodeIndex = 0;
		double bestHeuristic = std::numeric_limits<double>::max();
		for (size_t i = 0; i < eulerPath.size(); ++i) {
			size_t prevNode = eulerPath[0];
			if (i > 1) prevNode = eulerPath[i - 1];
			size_t currentNode = eulerPath[i];
			auto bestNext = findBestNext(prevNode, currentNode, graph, visitedEdges, heuristic);
			if (std::get<2>(bestNext)) {  // if it has unconnected edges
				if (std::get<1>(bestNext) <= bestHeuristic) {
					bestHeuristic = std::get<1>(bestNext);
					bestStartNodeIndex = i;
				}
			}
		}
		startNodeIndex = bestStartNodeIndex;
		startNode = eulerPath[startNodeIndex];
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

std::vector<size_t> readPath(std::string filePath) {
	std::fstream infile(filePath);
	size_t i;
	std::vector<size_t> p;
	while (infile >> i) p.push_back(i);
	return p;
}

int main(int argc, char** argv) {
	if (argc <= 1) {
		std::cerr << "Usage: " << argv[0] << " graphFilePath" << std::endl;
		return 1;
	}

	auto dumbHeuristic = [](const auto& g, size_t p, size_t c, size_t n) { return 0; };
	auto angleHeuristic = [](const auto& g, size_t p, size_t c, size_t n) {
		auto prevEdge = g[c].coords - g[p].coords;
		auto nextEdge = g[n].coords - g[c].coords;
		double d = prevEdge.normalized().dot(nextEdge.normalized());
		return -d + 1;
	};
	auto graph = readGraph(argv[1]);
	if (argc == 3) {
		std::cerr << "Computing score (lower = better)..." << std::endl;
		double sum = 0.0;
		double dist = 0;
		auto p = readPath(argv[2]);
		for (size_t i = 1; i < p.size() - 1; ++i) {
			sum += angleHeuristic(graph, p[i - 1], p[i], p[i + 1]);
			dist += (graph[p[i]].coords - graph[p[i - 1]].coords).length();
		}
		dist += (graph[p.back()].coords - graph[p[p.size() - 2]].coords).length();
		std::cerr << " total = " << sum << ". Distance = " << dist
		          << ", normalized cost = " << sum / dist << std::endl;
		return 0;
	}

	auto eulerPath = eulerOpt(graph, dumbHeuristic);
	for (auto i : eulerPath) {
		std::cout << i << std::endl;
	}
	return 0;
}
