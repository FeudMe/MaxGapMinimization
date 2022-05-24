public abstract class TreeSolver {
	/**
	 * Computes an l-satisfying selection for a tree of minimal size
	 * 
	 * @param T          the tree
	 * @param S          the selection array (size needs to be the number of edges)
	 * @param d          distances array (size needs to be number of vertices)
	 * @param l          the length to satisfy
	 * @param root_index the key of the subtree's root
	 * @return the longest unblocked distance from the subtree's root to a recursive
	 *         child
	 */
	public static double select_tree_DFS(Tree T, boolean[] S, double[] d, double l, int root_index) {
		Node root = T.get(root_index);
		if (root.degree() == 1) {
			return 0.0;
		}
		int e_min = -1;
		int e_max = -1;
		double d_min = Double.MAX_VALUE;
		double d_max = 0.0;
		for (Edge e : root.incident_edges.values()) {
			// if edge is edge to parent ...
			if (root.key != Math.min(e.first, e.second)) {
				continue;
			}
			// else process child ...
			Node c = T.get(Math.max(e.second, e.first));
			double d_c = select_tree_DFS(T, S, d, l, c.key) + e.length;
			// store in array ...
			d[c.key] = d_c;
			if (d_c > l) {
				S[e.index] = true;
				d[c.key] = 0.0;
			}
			if (d_c > l / 2) {
				if (d_c < d_min) {
					if (e_min != -1) {
						S[e_min] = true;
					}
					e_min = e.index;
					d_min = d_c;
				} else {
					S[e.index] = true;
				}
			} else if (d_c > d_max) {
				e_max = e.index;
				d_max = d_c;
			}
		}
		if (e_max != e_min && d_max + d_min > l) {
			// System.out.println(d_max + " is dmax and " + d_min + " is dmin");
			if (e_min != -1) {
				S[e_min] = true;
			}
			return d_max;
		} else {
			return d_min;
		}
	}
	
	public static void parametric_search(Tree T, int k) {
		double[] possible_distances = new double[T.nodes.size() * T.nodes.size()];
	}
	
	public static void possible_gap_lengths(Tree T, double possible_distances) {
		
	}

	public static void select_greedy(Tree T, boolean[] S, int k) {
		int first_sel = (int) Math.random() * T.edges.size();
		S[first_sel] = true;
		for (int i = 1; i < k; i++) {
			Tupel<Integer, Double> t = farthest_point_edge(T, S);
			S[t.first] = true;
		}
	}
	
	public static double eval_LSP(Tree T, boolean[] S) {
		double lsp = 0.0;

		for (Node n : T.nodes) {
			double d_max = max_DFS(n.key, T, S);
			lsp = Math.max(d_max, lsp);
		}
		
		return lsp;
	}
	
	public static double max_DFS(int root_index, Tree T, boolean[] S) {
		Node root = T.get(root_index);
		if (root.degree() == 1) {
			return 0.0;
		}
		double d_max = 0.0;
		for (Edge e : root.incident_edges.values()) {
			if (root.key != Math.min(e.first, e.second) || S[e.index]) {
				continue;
			} 
			
			Node c = T.get(Math.max(e.second, e.first));
			double d = max_DFS(c.key, T, S) + e.length;
			d_max = Math.max(d_max, d);
		}
		
		return d_max;
	}

	public static Tupel<Integer, Double> farthest_point_edge(Tree T, boolean[] S) {
		double[] d = new double[T.nodes.size()];
		run_dijkstra(T, S, d);
		Edge e_max = null;
		double d_max = -0.0;
		for (Edge e : T.edges) {
			if (!S[e.index]) {
				if (d[e.first] + e.length + d[e.second] > d_max) {
					e_max = e;
					d_max = d[e.first] + e.length + d[e.second];
				}
			}
		}
		return new Tupel<Integer, Double>(e_max.index, d_max);
	}

	static void run_dijkstra(Tree G, boolean[] S, double[] distance) {
		MHPriorityQueue Q = new MHPriorityQueue(G.nodes.size());
		for (int i = 0; i < distance.length; i++) {
			distance[i] = Double.MAX_VALUE;
		}
		for (Edge e : G.edges) {
			if (S[e.index]) {
				distance[e.first] = 0.0;
				distance[e.second] = 0.0;
			}
		}
		for (int i = 0; i < distance.length; i++) {
			if (distance[i] == 0.0) {
				Q.insert(G.get(i), 0.0);
			}
		}

		while (Q.size() > 0) {
			Tupel<Double, Node> t = Q.extract_min();
			Node n = t.second;
			int key = n.key;

			if (distance[key] < t.first) {
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.nodes.get(e.first).key;
				int index_second = G.nodes.get(e.second).key;

				double alt = distance[key] + e.length;

				if (e.first == key) {
					if (distance[index_second] > alt) {
						Q.insert(G.nodes.get(e.second), distance[index_second] = alt);
					}
				} else if (e.second == key) {
					if (distance[index_first] > alt) {
						Q.insert(G.nodes.get(e.first), distance[index_first] = alt);
					}
				}
			}
		}
	}

}
