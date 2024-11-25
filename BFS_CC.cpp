#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/internal/get_time.h>

#include "BFS_CC.h"
#include "examples/helper/graph_utils.h"

// Driver
int main(int argc, char* argv[]) {
  using graph = parlay::sequence<parlay::sequence<int>>;
  using utils = graph_utils<int>;
  auto usage = "Usage: BFS_CC <n> || BFS_CC <filename>";
  if (argc != 2) std::cout << usage << std::endl;
  else {
    long n = 0;
    graph G;
    try { n = std::stol(argv[1]); }
    catch (...) {}
    if (n == 0) {
      G = utils::read_symmetric_graph_from_file(argv[1]);
      n = G.size();
    } else {
      G = graph_utils<int>::rmat_symmetric_graph(n, 20*n);
    }
    utils::print_graph_stats(G);

    std::cout << std::endl;

    auto result = parlay::sequence<int>();
    parlay::internal::timer t("Time");
    for (int i=0; i < 5; i++) {
      result = CC(G);
      t.next("");
    }
  
    auto diff_labels = parlay::remove_duplicates(result);
    std::cout << "number of connected componets: " << diff_labels.size() << std::endl;
  }
}
