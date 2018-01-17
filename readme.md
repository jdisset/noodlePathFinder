# Build
Requirements: a C++14 compiler (g++ 5.1+, clang 3.4+, Visual Studio 2015+)

```
mkdir build
cd build && cmake .. && make
```

# Usage

`./noodlepath graphFile  > pathFile` to generate a pathFile from a graphFile

`./noodlepath graphFile pathFile` to print out an evaluation of the pathFile (turning costs + total distance) 
