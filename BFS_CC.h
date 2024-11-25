#include <atomic>

#include "parlay/primitives.h"
#include "parlay/sequence.h"

// Convential Top-Down BFS
auto BFS(int start, const parlay::sequence<parlay::sequence<int>>& G, parlay::sequence<std::atomic<int>>& labels) {
  parlay::sequence<int> frontier(1,start);
  while (frontier.size() > 0) {
    auto out = flatten(map(frontier, [&] (int u) {return G[u];}));
    frontier = filter(out, [&] (int v) {
      int expected = -1;
      return (labels[v] == -1) && labels[v].compare_exchange_strong(expected, start);});
  }
}

// Connected Components Labeling - Iterate over the vertices of the graph,
// if a vertex is unlabeled, it is labeled with it's own id and a parallel
// BFS is performed starting at the vertex and propgate the id to it's
// connected component
auto CC(const parlay::sequence<parlay::sequence<int>>& G) {
  auto labels = parlay::tabulate<std::atomic<int>>(G.size(), [&] (long i) {return -1;});

  for (int i = 0; i < G.size(); i++) {
    if (labels[i].load() == -1) {
      labels[i] = int(i);
      BFS(i, G, labels);
    }
  }

  auto label_values = parlay::tabulate(G.size(), [&] (long i) {return labels[i].load();});
  return label_values;
}
