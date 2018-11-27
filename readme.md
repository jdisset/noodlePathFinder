# Build
Requirements: a C++14 compiler (g++ 5.1+, clang 3.4+, Visual Studio 2015+)

```
mkdir build
cd build && cmake .. && make
```

# Usage

`./noodlepath graphFile  > pathFile` to generate a pathFile from a graphFile

`./noodlepath graphFile pathFile` to print out an evaluation of the pathFile (turning costs + total distance) 

Complete Usage:
  NoodlePath [OPTION...]

```
  -a, --algorithm arg  Algorithm to use: fleury or |hierholzer|
  -h, --heuristic arg  Heuristic to use: |angle|, angleDistance or random
  -g, --graphFile arg  Graph file path
  -p, --pathFile arg   Path file: will print stats about the given path
  -s, --startNode arg  Choose start node. Default is 0 or first odd degree
                       node
      --help           Print help
```

# Syntax of a graph file

``` 
NODE_ID X Y Z [list of connected NODE_ID]
```

## example
```
0 0 0 0 1 2 4
1 1 0 0 0 4 2
2 1 1 0 3 4 0 1
3 0.5 2 0 4 2
4 0 1 0 0 1 2 3
```

