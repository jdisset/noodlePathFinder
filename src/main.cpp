#include <functional>
#include <random>
#include <vector>
#include "cxxopts.hpp"
#include "graph.hpp"
#include "vec.hpp"

std::vector<size_t> readPath(std::string filePath) {
	std::fstream infile(filePath);
	size_t i;
	std::vector<size_t> p;
	while (infile >> i) p.push_back(i);
	return p;
}

int main(int argc, const char** argv) {
	std::string pathFile, graphFile;
	std::string algorithm = "hierholzer";
	std::string heuristic = "angle";
	int startNode = -1;

	cxxopts::Options options("NoodlePath",
	                         "When beetles fight these battles in a bottle with their "
	                         "paddles and the bottle's on a poodle and the poodle's eating "
	                         "noodles, they call this a muddle puddle tweetle poodle "
	                         "beetle noodle bottle paddle battle.");
	options.add_options()("a,algorithm", "Algorithm to use: fleury or |hierholzer|",
	                      cxxopts::value(algorithm))(
	    "h,heuristic", "Heuristic to use: |angle|, angleDistance or random",
	    cxxopts::value(heuristic))("g,graphFile", "Graph file path",
	                               cxxopts::value(graphFile))(
	    "p,pathFile", "Path file: will print stats about the given path",
	    cxxopts::value(pathFile))(
	    "s,startNode", "Choose start node. Default is 0 or first odd degree node",
	    cxxopts::value(startNode))("help", "Print help");

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	// heuristic takes graph, prevNode, currentNode, nextNode
	static std::random_device
	    rd;  // Will be used to obtain a seed for the random number engine
	static std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
	static std::uniform_real_distribution<double> dis(0.0, 1.0);
	auto randomHeuristic = [](const auto&, size_t, size_t, size_t) { return dis(gen); };
	auto angleDistHeuristic = [](const auto& g, size_t p, size_t c, size_t n) {
		auto prevEdge = g[c].coords - g[p].coords;
		auto nextEdge = g[n].coords - g[c].coords;
		double pl = prevEdge.length();
		double nl = nextEdge.length();
		double d = prevEdge.normalized().dot(nextEdge.normalized());
		double l = pl + nl;
		if (l == 0) l = 1;
		return (1 - d) / l;
	};
	auto angleHeuristic = [](const auto& g, size_t p, size_t c, size_t n) {
		auto prevEdge = g[c].coords - g[p].coords;
		auto nextEdge = g[n].coords - g[c].coords;
		double d = prevEdge.normalized().dot(nextEdge.normalized());
		return (1 - d);
	};

	if (graphFile.empty()) {
		std::cerr << "Please give a valid graph file path" << std::endl;
		return 1;
	}
	auto graph = Graph::fromFile(graphFile);

	std::function<double(const std::vector<GraphNode>&, size_t, size_t, size_t)> H;
	if (heuristic == "angle")
		H = angleHeuristic;
	else if (heuristic == "angleDistance")
		H = angleDistHeuristic;
	else if (heuristic == "random")
		H = randomHeuristic;
	else {
		std::cerr << "Invalid heuristic" << std::endl;
		return 1;
	}

	auto odds = graph.listOdds();
	size_t nbOdds = odds.size();

	if (nbOdds > 2) {
		std::cerr << "Too many odd degree nodes (" << nbOdds
		          << "), no Euler path possible. Aborting." << std::endl;
		return 1;
	}

	if (!pathFile.empty()) {
		std::cerr << "Number of odds vertex = " << nbOdds
		          << " ; Euler path total distance = " << graph.computeEulerDistance()
		          << std::endl;
		std::cerr << "Score with " << heuristic
		          << " heuristic = " << graph.computePathScore(readPath(pathFile), H)
		          << std::endl;
		return 0;
	}

	if (startNode == -1) {
		if (nbOdds > 0)
			startNode = *(odds.begin());
		else
			startNode = 0;
	} else {
		if (nbOdds > 0) {
			if (!odds.count(startNode)) {
				std::cerr << "Specified start node is not odd degree but graph has odd degrees "
				             "node. Odd degree nodes are :";
				for (auto& n : odds) std::cerr << n << " ";
				std::cerr << std::endl;
				return 1;
			}
		}
	}

	std::vector<size_t> eulerPath;
	if (algorithm == "fleury") {
		eulerPath = graph.eulerFleury(H, startNode);
	} else if (algorithm == "hierholzer") {
		eulerPath = graph.eulerHierholzer(H, startNode);
	} else if (algorithm == "exhaustive") {
		eulerPath = graph.eulerExhaustive(H, startNode);
	} else {
		std::cerr << "invalid algorithm" << std::endl;
		return 1;
	}
	for (auto& e : eulerPath) std::cout << e << std::endl;
	return 0;
}
