import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;

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
	public static double select_tree_DFS(Tree T, boolean[] S, double[] d, double l, int root_index, int abs_root,
			boolean[] visited) {
		Node root = T.get(root_index);
		if (root.degree() == 1 && root_index != abs_root) {
			// System.out.println(root_index + " returned early");
			return 0.0d;
		}
		visited[root.key] = true;

		int e_above = -1;
		int e_below = -1;
		double d_above = Double.MAX_VALUE;
		double d_below = 0.0;
		for (Edge e : root.incident_edges.values()) {
			// if edge is edge to parent ...
			if (visited[e.first] && visited[e.second]) {
				continue;
			}
			// else process child ...
			Node c = (visited[e.first]) ? T.nodes.get(e.second) : T.nodes.get(e.first);
			double d_c = select_tree_DFS(T, S, d, l, c.key, abs_root, visited) + e.length;
			// store in array ...
			d[c.key] = d_c;
			if (d_c > l) {
				S[e.index] = true;
				d[c.key] = 0.0;
				d_c = 0.0;
			}
			if (d_c > l / 2) {
				if (d_c < d_above) {
					if (e_above != -1) {
						S[e_above] = true;
					}
					e_above = e.index;
					d_above = d_c;
				} else {
					S[e.index] = true;
					d[c.key] = 0.0;
				}
			} else if (d_c >= d_below) {
				e_below = e.index;
				d_below = d_c;
			}
		}
		if (e_below != e_above && d_below + d_above > l) {
			// System.out.println(d_max + " is dmax and " + d_min + " is dmin");
			if (e_above != -1) {
				S[e_above] = true;
			}
			return d_below;
		} else {
			return d_above;
		}
	}

	public static Tupel<Double, Double> eval_tree_DFS(Tree T, boolean[] S, int root_index, int abs_root,
			boolean[] visited) {
		Node root = T.get(root_index);
		if (root.degree() == 1 && root_index != abs_root) { // 
			// System.out.println(root_index + " returned early");
			return new Tupel<Double, Double>(0.0, 0.0);
		}
		visited[root.key] = true;

		double d_max = 0.0;
		double d_max_2 = 0.0;
		double rec_max = 0.0;
		for (Edge e : root.incident_edges.values()) {
			// if edge is edge to parent ...
			if ((visited[e.first] && visited[e.second])) {
				continue;
			}
			// else process child ...
			Node c = (visited[e.first]) ? T.nodes.get(e.second) : T.nodes.get(e.first);
			Tupel<Double, Double> d_c = eval_tree_DFS(T, S, c.key, abs_root, visited);
			d_c.first = (S[e.index]) ? 0.0 : d_c.first + e.length;

			rec_max = Math.max(rec_max, d_c.second);
			if (d_c.first > d_max) {
				d_max_2 = d_max;
				d_max = d_c.first;
			} else if (d_c.first > d_max_2) {
				d_max_2 = d_c.first;
			}

		}
		return new Tupel<Double, Double>(d_max, Math.max(d_max + d_max_2, rec_max));

	}

	public static boolean[] parametric_search(Tree T, int k, int root) {
		double[] possible_distances = new double[T.nodes.size() * T.nodes.size()];
		possible_gap_lengths(T, possible_distances);
		Arrays.sort(possible_distances);
		LinkedList<Double> dists = new LinkedList<Double>();
		for (double d : possible_distances) {
			dists.add(d);
		}
		Utility.removeDuplicates(dists);

		int left = 0;
		int right = dists.size() - 1;
		double d = 0.0;
		boolean[] S_opt = new boolean[T.nodes.size()];
		while (left < right) {
			int middle = (left + right) / 2;
			boolean[] S = new boolean[T.edges.size()];
			d = dists.get(middle);
			select_tree_DFS(T, S, new double[T.nodes.size()], d, root, root, new boolean[T.nodes.size()]);
			int n = Utility.num_selected(S);
			// System.out.println(n + ", dist: " + d);
			if (n <= k && d > 0.0) {
				S_opt = Arrays.copyOf(S, S.length);
				right = middle - 1;
			} else {
				left = middle + 1;
			}
		}
		// System.out.println(left + " " + right + " " + d + " " + dists.get(left));
		return S_opt;

	}

	public static boolean[] k_largest(Tree T, int k) {
		boolean[] S = new boolean[T.edges.size()];
		T.edges.sort((Edge e_1, Edge e_2) -> {
			return (e_1.length > e_2.length) ? 1 : -1;
		});
		for (int i = 0; i < k; i++) {
			S[T.edges.get(i).index] = true;
		}
		return S;
	}

	public static void possible_gap_lengths(Tree T, double[] possible_distances) {
		for (Node n : T.nodes) {
			boolean[] visited = new boolean[T.nodes.size()];
			distancesDFS(T, possible_distances, n.key, visited, n.key, 0);
		}
	}

	public static void distancesDFS(Tree T, double[] distances, int root_index, boolean[] visited, int source,
			double d) {
		visited[root_index] = true;
		distances[T.nodes.size() * source + root_index] = d;
		distances[T.nodes.size() * root_index + source] = d;
		if (T.get(root_index).degree() == 1) {
			return;
		}
		for (Edge e : T.get(root_index).incident_edges.values()) {
			int c_key = (e.first == root_index) ? e.second : e.first;
			if (visited[c_key]) {
				continue;
			}
			distancesDFS(T, distances, c_key, visited, source, d + e.length);
		}
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
			boolean[] visited = new boolean[T.nodes.size()];
			double d_max = max_DFS(n.key, n.key, T, S, visited);
			lsp = Math.max(d_max, lsp);
		}

		return lsp;
	}

	public static double max_DFS(int root_index, int root_abs, Tree T, boolean[] S, boolean[] visited) {
		visited[root_index] = true;
		Node root = T.get(root_index);
		if (root.degree() == 1 && root_index != root_abs) {
			return 0.0;
		}
		double d_max = 0.0;
		for (Edge e : root.incident_edges.values()) {
			int c_key = (e.first == root_index) ? e.second : e.first;
			if (visited[c_key] || S[e.index]) {
				continue;
			}

			Node c = T.get(c_key);
			double d = max_DFS(c.key, root_abs,T, S, visited) + e.length;
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
