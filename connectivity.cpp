#include <iostream>
#include <string>

#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include "connectivity.h"

// Driver
int main(int argc, char* argv[]) {
  using graph = parlay::sequence<parlay::sequence<int>>;
  auto usage = "Usage: connectivity <beta> <n> <permute> || connectivity <beta> <filename> <permute>";
  if (argc != 4) std::cout << usage << std::endl;
  else {
    float beta = 0.5; // defualt beta
    long n = 0;
    int permute = 0;
    try {permute = std::stoi(argv[3]);} catch (...) {}
    try {beta = std::atof(argv[1]);} catch(...) {}
    graph G;
    try { n = std::stol(argv[2]); }
    catch (...) {}
    if (n == 0) {
      G = graph_utils<int>::read_symmetric_graph_from_file(argv[2]);
      n = G.size();
    } else {
      G = graph_utils<int>::rmat_symmetric_graph(n, 20*n);
    }

    graph_utils<int>::print_graph_stats(G);
    parlay::sequence<int> result;
    parlay::internal::timer t("Time");
    std::cout << std::endl;

    
    for (int i=0; i < 5; i++) {
      result = connectivity(G, beta, 0, permute);
      t.next("");
    }
    
    auto diff_labels = parlay::remove_duplicates(result);
    std::cout << "number of connected componets: " << diff_labels.size() << std::endl;
  }
}
