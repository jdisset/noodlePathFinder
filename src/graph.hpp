#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <unordered_set>
#include <utility>
#include "vec.hpp"

struct pair_hash {
	// just a quick&dirty hash function so that unordered_pair is happy to store pairs of
	// uint
	inline std::size_t operator()(const std::pair<size_t, size_t>& v) const {
		return v.first * 1000000000u + v.second;
	}
};

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

class Graph {
	std::vector<GraphNode> graph;

 public:
	static Graph fromFile(std::string filePath) {
		Graph g;
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
			g.graph.push_back(node);
		}
		return g;
	}

	size_t size() const { return graph.size(); }
	void addNode(const GraphNode& n) { graph.push_back(n); }

	template <typename H>
	double computePathScore(const std::vector<size_t>& p, const H& heuristic) const {
		double sum = 0.0;
		double dist = 0;
		for (size_t i = 1; i < p.size() - 1; ++i) {
			sum += heuristic(graph, p[i - 1], p[i], p[i + 1]);
			dist += (graph[p[i]].coords - graph[p[i - 1]].coords).length();
		}
		dist += (graph[p.back()].coords - graph[p[p.size() - 2]].coords).length();
		return sum / dist;
	}
	double computePathDistance(const std::vector<size_t>& p) const {
		double dist = 0;
		for (size_t i = 1; i < p.size(); ++i) {
			dist += (graph[p[i]].coords - graph[p[i - 1]].coords).length();
		}
		return dist;
	}

	size_t nbOdds() {
		size_t nbOdds = 0;
		for (const auto& n : graph)
			if (n.connected.size() % 2 == 1) ++nbOdds;
		return nbOdds;
	}
	std::set<size_t> listOdds() {
		std::set<size_t> res;
		for (size_t i = 0; i < graph.size(); ++i) {
			if (graph[i].connected.size() % 2 == 1) res.insert(graph[i].id);
		}
		return res;
	}

	double computeEulerDistance() const {
		// just the sum of all edges distance
		std::unordered_set<std::pair<size_t, size_t>, pair_hash> visitedEdges;
		for (const auto& n : graph)
			for (size_t c : n.connected) visitedEdges.insert(std::minmax(n.id, c));
		double dist = 0;
		for (auto& e : visitedEdges)
			dist += (graph[e.first].coords - graph[e.second].coords).length();
		return dist;
	}

	size_t getNumberOfEdges() const {
		std::unordered_set<std::pair<size_t, size_t>, pair_hash> visitedEdges;
		for (const auto& n : graph)
			for (size_t c : n.connected) visitedEdges.insert(std::minmax(n.id, c));
		return visitedEdges.size();
	}

	void removeEdge(size_t a, size_t b) {
		auto& va = graph[a].connected;
		auto& vb = graph[b].connected;
		va.erase(std::remove(va.begin(), va.end(), b), va.end());
		vb.erase(std::remove(vb.begin(), vb.end(), a), vb.end());
	}

	void addEdge(size_t a, size_t b) {
		auto& aCo = graph[a].connected;
		if (std::find(aCo.begin(), aCo.end(), b) == aCo.end()) {
			aCo.push_back(b);
			graph[b].connected.push_back(a);
		}
	}

	void getReachable(size_t n, std::unordered_set<size_t>& visited) const {
		visited.insert(n);
		for (const auto& adj : graph[n].connected) {
			if (!visited.count(adj)) {
				getReachable(adj, visited);
			}
		}
	}

	std::unordered_set<size_t> getReachable(size_t n) const {
		std::unordered_set<size_t> reachable;
		getReachable(n, reachable);
		return reachable;
	}

	bool isBridge(size_t a, size_t b, size_t N) {
		// to see if the edge a -- b is a bridge, we delete it and
		// check if the graph is still fully connected
		removeEdge(a, b);
		auto sa = getReachable(a).size();
		auto sb = getReachable(b).size();
		addEdge(a, b);
		return (sa != N || sb != N);
	}

	bool isBridge(size_t a, size_t b) { return isBridge(a, b, graph.size()); }

	GraphNode& operator[](std::size_t idx) { return graph[idx]; }

	const GraphNode& operator[](std::size_t idx) const { return graph[idx]; }
	/********************************************************/
	/*                                                      */
	/*                 Hierholzer's version                 */
	/*                                                      */
	/********************************************************/

	// finds the best unvisited edge. Returns -1 if all visited
	template <class V, class F>
	std::tuple<size_t, double, bool> findBestNext(size_t prevNode, size_t currentNode,
	                                              const V& visited, const F& heuristic) {
		long long int nextNode = -1;
		double minHeuristic = std::numeric_limits<double>::max();
		for (auto& nId : graph[currentNode].connected) {
			if (!visited.count(std::minmax(currentNode, nId))) {
				double eurVal = heuristic(graph, prevNode, currentNode, nId);
				if (eurVal < minHeuristic) {
					nextNode = static_cast<long long int>(nId);
					minHeuristic = eurVal;
				}
			}
		}
		bool hasUnvisitedEdges = (nextNode >= 0);
		return std::make_tuple(nextNode, minHeuristic, hasUnvisitedEdges);
	}
	template <typename F>
	std::vector<size_t> eulerOpt(const F& heuristic, size_t startNode = 0) {
		// a variation of hierholzer's algorithm with an added heuristic function

		const size_t totalEdgeNumber = getNumberOfEdges();
		// first we check that the graph is (semi)eulerian
		std::vector<size_t> oddNodes;
		for (const auto& n : graph)
			if (n.connected.size() % 2 == 1) oddNodes.push_back(n.id);
		if (oddNodes.size() > 2) {
			std::cout << "This graph has " << oddNodes.size() << " odds. Aborting."
			          << std::endl;
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
				auto bestNext = findBestNext(prevNode, currentNode, visitedEdges, heuristic);
				size_t nextNode = std::get<0>(bestNext);
				assert(std::get<2>(bestNext));
				visitedEdges.insert(std::minmax(currentNode, nextNode));
				subCircle.push_back(nextNode);
			} while (subCircle.back() != startNode);

			// subcircle formed
			// we first replace startNode in eulerGraph by this new subcircle
			eulerPath.insert(eulerPath.begin() + (long)startNodeIndex, subCircle.begin(),
			                 subCircle.end());
			// we then need to find the next unvisited edge with the best possible next edge
			// from which to start
			size_t bestStartNodeIndex = 0;
			double bestHeuristic = std::numeric_limits<double>::max();
			for (size_t i = 0; i < eulerPath.size(); ++i) {
				size_t prevNode = eulerPath[0];
				if (i > 1) prevNode = eulerPath[i - 1];
				size_t currentNode = eulerPath[i];
				auto bestNext = findBestNext(prevNode, currentNode, visitedEdges, heuristic);
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

	/********************************************************/
	/*                                                      */
	/*                   Fleury's version                   */
	/*                                                      */
	/********************************************************/

	template <typename F>
	std::vector<size_t> eulerFleury(const F& heuristic, size_t startNode = 0) const {
		Graph gcopy = *this;
		size_t N = graph.size();

		std::vector<size_t> eulerPath;
		size_t prevNode = startNode;
		size_t currentNode = startNode;
		eulerPath.push_back(currentNode);

		do {
			auto& adj = gcopy[currentNode].connected;
			double minHeuristic = std::numeric_limits<double>::max();
			size_t bestNext = adj[0];
			for (size_t candidateNext : adj) {
				if (!gcopy.isBridge(currentNode, candidateNext, N)) {
					double h = heuristic(gcopy.graph, prevNode, currentNode, candidateNext);
					if (h <= minHeuristic) {
						minHeuristic = h;
						bestNext = candidateNext;
					}
				}
			}
			gcopy.removeEdge(currentNode, bestNext);
			if (gcopy[currentNode].connected.size() == 0) --N;
			if (gcopy[bestNext].connected.size() == 0) --N;
			prevNode = currentNode;
			currentNode = bestNext;
			eulerPath.push_back(currentNode);
		} while (gcopy[currentNode].connected.size() > 0);

		return eulerPath;
	}

	// std::vector<std::vector<size_t>> getAllCycles(Graph g, size_t currentNode = 0,
	// size_t startNode = 0) {
	// if(currentNode == startNode && acc.size()>0) acc.
	// for (auto& n : g[startNode].connected) {
	// g.removeEdge(currentNode, n);
	// getAllCycles(g, n, startNode);
	//}
	/*}*/
};

#endif
