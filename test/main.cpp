#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../src/graph.hpp"
#include "catch.hpp"

TEST_CASE("Factorials are computed", "[factorial]") {
	// https://www.geeksforgeeks.org/wp-content/uploads/Bridge1.png
	Graph graph;
	graph.addNode(GraphNode{0, {1, 1, 0}});
	graph.addNode(GraphNode{1, {0, 1, 0}});
	graph.addNode(GraphNode{2, {0, 0, 0}});
	graph.addNode(GraphNode{3, {2, 1, 0}});
	graph.addNode(GraphNode{4, {2, 0, 0}});

	graph[0].connected = {1, 2, 3};
	graph[1].connected = {0, 2};
	graph[2].connected = {1, 0};
	graph[3].connected = {0, 4};
	graph[4].connected = {3};

	REQUIRE(graph.getNumberOfEdges() == 5);

	REQUIRE(graph.getReachable(0).size() == 5);
	REQUIRE(graph.getReachable(1).size() == 5);
	REQUIRE(graph.getReachable(2).size() == 5);
	REQUIRE(graph.getReachable(3).size() == 5);
	REQUIRE(graph.getReachable(4).size() == 5);

	graph.removeEdge(0, 3);
	REQUIRE(graph[0].connected == std::vector<size_t>{1, 2});
	REQUIRE(graph[3].connected == std::vector<size_t>{4});
	REQUIRE(graph.getNumberOfEdges() == 4);
	REQUIRE(graph.getReachable(3).size() == 2);
	REQUIRE(graph.getReachable(4).size() == 2);
	REQUIRE(graph.getReachable(0).size() == 3);
	REQUIRE(graph.getReachable(1).size() == 3);
	REQUIRE(graph.getReachable(2).size() == 3);

	graph.addEdge(0, 3);

	REQUIRE(graph.isBridge(0, 3));
	REQUIRE(graph.isBridge(3, 0));
	REQUIRE(graph.isBridge(4, 3));
	REQUIRE(graph.isBridge(3, 4));

	REQUIRE(!graph.isBridge(1, 0));
	REQUIRE(!graph.isBridge(2, 0));
	REQUIRE(!graph.isBridge(0, 1));
	REQUIRE(!graph.isBridge(0, 2));
	REQUIRE(!graph.isBridge(2, 1));
	REQUIRE(!graph.isBridge(1, 2));

	auto angleHeuristic = [](const auto& g, size_t p, size_t c, size_t n) {
		auto prevEdge = g[c].coords - g[p].coords;
		auto nextEdge = g[n].coords - g[c].coords;
		double d = prevEdge.normalized().dot(nextEdge.normalized());
		return -d + 1;
	};

	graph.addEdge(0, 4);

	auto p = graph.eulerFleury(angleHeuristic, 1);
	for (auto n : p) std::cerr << n << " ";
	std::cerr << std::endl;
	REQUIRE(graph.computePathDistance(p) == graph.computeEulerDistance());
}
