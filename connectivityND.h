#include <atomic>
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "examples/helper/graph_utils.h"

// Same algorithms as connectivity.h minus the use of delayed sequences

// Relabel Clusters
auto relabel(parlay::sequence<int>& clusters) {
	long n = clusters.size();

	auto map = parlay::tabulate(n+1, [&] (int i) {return 0;});
	// generate an array of size n+1 with 0's in each entry and set each index that 
	// corresponds to a cluster id to 1
	parlay::parallel_for(0, n, [&] (int i) {
		if (!map[clusters[i]]) { 
			map[clusters[i]] = 1;
		}
	});

	// generates new labels between 0 and number of unique clusters minus one
	parlay::scan_inplace(parlay::make_slice(map));

	int num_unique = map[n];
	parlay::parallel_for(0, n, [&] (int i) {
		// Update cluster labels
		clusters[i] = map[clusters[i]];
	});

	map.clear();
	return num_unique;
}

// Low Diameter Decomopostion
auto cluster(const parlay::sequence<parlay::sequence<int>>& G, float beta, bool permute) {
	long n = G.size();

	parlay::sequence<int> vertex_permutation;
	if (permute) {
		vertex_permutation = parlay::random_permutation<int>(n);
	}

	// Cluster labels
	auto cluster = parlay::tabulate<std::atomic<int>>(n, [&] (int i) {return -1;});

	parlay::random_generator gen(0);
	std::exponential_distribution<float> exp(beta);
	auto exp_samples = parlay::tabulate(n, [&] (long i) {
		auto r = gen[i];
		return (int) std::floor(exp(r));
	});

	int max_sample = parlay::reduce(exp_samples, parlay::maximum<int>());

	// Bucket Vertices based off exponential shifts
	auto buckets = parlay::group_by_index(parlay::tabulate(n, [&] (int i) {
		int id = (permute) ? vertex_permutation[i] : i;
		return std::pair(max_sample - exp_samples[id], id);
	}), max_sample + 1);

	auto cond = [&] (int i) {return cluster[i].load() == -1;};
	auto update = [&] (int s, int d) -> bool {
		int expected = -1;
		return cluster[d].compare_exchange_strong(expected, cluster[s]);
	};

	parlay::sequence<int> frontier;

	for (int i = 0; i < buckets.size(); i++) {
		auto new_clusters = parlay::filter(buckets[i], [&] (int j) {
			if (cluster[j].load() != -1) return false;
			cluster[j] = j;
			return true;
		});

		frontier = parlay::append(frontier, new_clusters);
		// Collect neighbors of vertices on the frontier 
		auto pairs = parlay::flatten(parlay::map(frontier, [&] (int v) {
			return parlay::map(G[v], [=] (int u) {
				return std::pair(v,u);
			});
		}));

		// Keep vertices that are unvisited and assign each to the cluster that succeded in setting 
		// the vertices id
		frontier = parlay::map_maybe(pairs, [&] (auto p) {
			auto [v,u] = p;
			return ((cond(u) && update(v,u)) ? std::make_optional(u) : std::nullopt);
		});

	}
	return parlay::map(cluster, [] (auto& c) {return c.load();});
}

//***********************************************************************************************************//
// Given a sequence of cluster id's between 0 and the number of unique clusters minus one, this              //
// function contracts each cluster into a single node and removes all duplicate edges and self               //
// loop. This function returns the contracted graph and two additional sequences. The first                  //
// is a sequence of flags that can be used to determine if a contracted node is singular (no adjacency).     //
// Specifically, flags[i] == flags[i+1] iff cluster i contracts into a singular node. The second is the      //
// the inverse map for flags. Specifically, vertex i (assuming non-singular) in the contracted graph         //
// represents cluster map[i]. Additionally, we maintain another mapping which is an extension of map such    //
// that inv_map[i] denotes the new label of cluster i if cluster i is non-singular. Otherwise inv_map[i]     //
// is undefined. The techniques for maintaining labels during contraction steps are orginally from the GBBS  //
// implementation of this LLD based connectivity algorithm. See https://github.com/ParAlg/gbbs/ for more     //
// details.                                                                                                  //
//***********************************************************************************************************//
auto contract(const parlay::sequence<parlay::sequence<int>>& G, const parlay::sequence<int>& clusters, int num_clusters, bool* stop) {
	// fetch inter-cluster edges
	auto potential_edges = parlay::flatten(parlay::tabulate(G.size(), [&](int u) {
        return parlay::map_maybe(G[u], [&](int v) {
            if (clusters[u] < clusters[v]) { // inter-cluster edge
                return std::optional<std::pair<int, int>>(std::pair(clusters[u], clusters[v]));
            }
            return std::optional<std::pair<int, int>>(std::nullopt);
        });
    }));
	
	auto edges = parlay::remove_duplicates(potential_edges);
	auto flags = parlay::tabulate(num_clusters + 1, [&] (int i) {return 0;});

	parlay::parallel_for(0, edges.size(), [&] (int i) {
		auto [u,v] = edges[i];
		if (!flags[u]) flags[u] = 1;
		if (!flags[v]) flags[v] = 1; 
	});
	parlay::scan_inplace(parlay::make_slice(flags));

	int num_ns_clusters = flags[num_clusters];

	if (num_ns_clusters == 0) { // all contracted nodes are singular
		parlay::sequence<parlay::sequence<int>> Gc = parlay::tabulate(num_clusters, [&] (int i) {
			return parlay::sequence<int>();
		});
		*stop = true;
		return std::make_tuple(std::move(Gc), std::move(flags), std::move(parlay::sequence<int>()));
	} else {
		auto map = parlay::sequence<int>::uninitialized(num_ns_clusters);

		parlay::parallel_for(0, num_clusters, [&] (int i) {
			if (flags[i] != flags[i+1]) {
				map[flags[i]] = i;
			}
		});

		auto inv_map = parlay::tabulate(num_clusters, [&] (int i) {return -1;});
		parlay::parallel_for(0, num_ns_clusters, [&] (int i) {
			inv_map[map[i]] = i;
		});

		// update edge labels
		parlay::parallel_for(0, edges.size(), [&] (int i) {
			auto [u,v] = edges[i];
			edges[i] = std::pair(inv_map[u], inv_map[v]);
		});

		// init the contracted graph
		auto Gc = graph_utils<int>::symmetrize(parlay::group_by_index(edges, num_ns_clusters));
		edges.clear();
		inv_map.clear(); 
		return std::make_tuple(std::move(Gc), std::move(flags), std::move(map));
	}
}

//********************************************************************************************//
// This function is an implementation of the main algorithm given in the paper:				  //
// "A simple and practical linear-work parallel algorithm for connectivity" by                //
// Julian Shun, Laxman Dhulipala, and Guy E. Blelloch. Additioanlly, the main structure       //
// of the function takes heavily from the GBBS implemntation of the algorithm.                //
//********************************************************************************************//
parlay::sequence<int> connectivity(const parlay::sequence<parlay::sequence<int>>& G, float beta, int level, int permute) {
	long n = G.size();
	
	permute = (level > 0 && permute == 1) ? 1 : 0; 
	auto clusters = cluster(G, beta, permute);

	int num_clusters = relabel(clusters);

	bool stop = false;
	auto out = contract(G, clusters, num_clusters, &stop);
	if (stop) {
		return clusters;
	}
	auto Gc = std::move(std::get<0>(out));
	auto& flags = std::get<1>(out);
	auto& map = std::get<2>(out);

	auto new_labels = connectivity(Gc, beta, level + 1, permute);
	// update labels
	parlay::parallel_for(0, n, [&] (int i) {
		int cluster = clusters[i];
		int gc_cluster = flags[cluster];
		// if the contracted cluster was non-singular we need to update its label
		if (gc_cluster != flags[cluster + 1]) {
			clusters[i] = map[new_labels[gc_cluster]];
		}
	});

	new_labels.clear();
	flags.clear();
	map.clear();
	return clusters;
}
