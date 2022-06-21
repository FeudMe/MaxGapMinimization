import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;

public abstract class Graph_Eval {

	public static double evaluate_LSP(Graph G, boolean S[]) {
		double dist[][] = new double[G.num_nodes()][G.num_nodes()];
		int pre[][] = new int[G.num_nodes()][G.num_nodes()];
		short masked[][] = new short[G.num_nodes()][G.num_nodes()];

		for (int i = 0; i < G.num_nodes(); i++) {
			for (int j = 0; j < G.num_nodes(); j++) {
				if (j == i) {
					dist[i][j] = 0;
				} else {
					dist[i][j] = Double.MAX_VALUE;
				}
				pre[i][j] = -1;
				masked[i][j] = 0;
			}
		}

		for (Node v : G.get_nodes().values()) {
			int key = v.key;
			// define comparator

			// initialize priority queue
			MHPriorityQueue Q = new MHPriorityQueue(G.num_nodes());
			Q.insert(v, dist[key][key]);
			// until queue is empty ...
			while (Q.size() > 0) {
				// pop first element from queue
				Tupel<Double, Node> t = Q.extract_min();
				Node n = t.second;
				// if shortest path to node is not optimal -> skip
				if (dist[key][n.reduced_key] < t.first) {
					continue;
				}
				// forall incident edges ...
				for (Edge e : n.incident_edges.values()) {
					// get indices
					int index_first = G.get_nodes().get(e.first).reduced_key;
					int index_second = G.get_nodes().get(e.second).reduced_key;
					// get distance over edge
					double alt = dist[key][n.reduced_key] + e.length;
					// check if distance of n unvisited endpoint is smaller or larger
					if (e.first == n.key) {
						if (dist[key][index_second] > alt || (dist[key][index_second] >= alt
								&& masked[key][index_second] > masked[key][n.reduced_key])) {
							// update distance array
							dist[key][index_second] = alt;
							// update predecessor array
							pre[key][index_second] = n.reduced_key;
							// update masked array
							if (masked[key][n.reduced_key] != 0 || S[e.index]) {
								masked[key][index_second] = 1;
							} else {
								masked[key][index_second] = 0;
							}
							// push updated end vertex to priority queue
							if (!Q.contains(e.second)) {
								Q.insert(G.get_nodes().get(e.second), dist[key][index_second]);
							} else {
								Q.decrease_key(G.get_nodes().get(e.second), dist[key][index_second]);
							}
						}
					} else if (e.second == n.key) {
						if (dist[key][index_first] > alt || (dist[key][index_second] >= alt
								&& masked[key][index_second] > masked[key][n.reduced_key])) {
							// update distance array
							dist[key][index_first] = alt;
							// update predecessor array
							pre[key][index_first] = n.reduced_key;
							// update masked array
							if (masked[key][n.reduced_key] != 0 || S[e.index]) {
								masked[key][index_first] = 1;
							} else {
								masked[key][index_first] = 0;
							}
							// push updated end vertex to priority queue
							if (!Q.contains(e.first)) {
								Q.insert(G.get_nodes().get(e.first), dist[key][index_first]);
							} else {
								Q.decrease_key(G.get_nodes().get(e.first), dist[key][index_first]);
							}
						}
					}
				}
			}
		}
		double max_val = 0.0;
		for (int i = 0; i < G.num_nodes(); i++) {
			for (int j = 0; j < G.num_nodes(); j++) {
				if (masked[i][j] == 0 && dist[i][j] > max_val) {
					max_val = dist[i][j];
				}
			}
		}
		return max_val;
	}

	public static double evaluate_LSP_memory_optimised(Graph G, boolean S[]) {
		double dist[] = new double[G.num_nodes()];
		int pre[] = new int[G.num_nodes()];
		short masked[] = new short[G.num_nodes()];

		double max_val = 0.0;

		for (Node v : G.get_nodes().values()) {
			int key = v.key;
			MHPriorityQueue Q = new MHPriorityQueue(G.num_nodes());

			for (int i = 0; i < G.num_nodes(); i++) {
				if (key == i) {
					dist[i] = 0;
				} else {
					dist[i] = Double.MAX_VALUE;
				}
				pre[i] = -1;
				masked[i] = 0;
				// Q.insert(G.get_nodes().get(i), dist[i]);
			}

			// initialize priority queue
			Q.insert(v, dist[key]);

			int num_masked_in_Q = 0;

			// until queue is empty ...
			while (Q.size() > num_masked_in_Q) {
				// pop first element from queue
				Tupel<Double, Node> t = Q.extract_min();
				Node n = t.second;
				if (masked[n.reduced_key] == 1) {
					num_masked_in_Q--;
				}

				// forall incident edges ...
				for (Edge e : n.incident_edges.values()) {
					// get indices
					int index_first = G.get_nodes().get(e.first).reduced_key;
					int index_second = G.get_nodes().get(e.second).reduced_key;
					// get distance over edge
					double alt = dist[n.reduced_key] + e.length;
					// check if distance of n unvisited endpoint is smaller or larger
					if (e.first == n.key) {
						if (dist[index_second] > alt
								|| (dist[index_second] >= alt && masked[index_second] > masked[n.reduced_key])) {
							// update distance array
							dist[index_second] = alt;
							// update predecessor array
							pre[index_second] = n.reduced_key;
							boolean was_masked = masked[index_second] == 1;
							// update masked array
							if (masked[n.reduced_key] != 0 || S[e.index]) {
								masked[index_second] = 1;
							} else {
								masked[index_second] = 0;
							}
							// push updated end vertex to priority queue
							if (!Q.contains(e.second)) {
								Q.insert(G.get_nodes().get(e.second), dist[index_second]);
								if (masked[index_second] == 1) {
									num_masked_in_Q++;
								}
							} else {
								Q.decrease_key(G.get_nodes().get(e.second), dist[index_second]);
								if ((!was_masked) && masked[index_second] == 1) {
									num_masked_in_Q++;
								} else if (was_masked && masked[index_second] == 0) {
									num_masked_in_Q--;
								}
							}
						}
					} else if (e.second == n.key) {
						if (dist[index_first] > alt
								|| (dist[index_second] >= alt && masked[index_second] > masked[n.reduced_key])) {
							// update distance array
							dist[index_first] = alt;
							// update predecessor array
							boolean was_masked = masked[index_first] == 1;
							pre[index_first] = n.reduced_key;
							// update masked array
							if (masked[n.reduced_key] != 0 || S[e.index]) {
								masked[index_first] = 1;
							} else {
								masked[index_first] = 0;
							}
							// push updated end vertex to priority queue
							if (!Q.contains(e.first)) {
								Q.insert(G.get_nodes().get(e.first), dist[index_first]);
								if (masked[index_first] == 1) {
									num_masked_in_Q++;
								}
							} else {
								Q.decrease_key(G.get_nodes().get(e.first), dist[index_first]);
								if ((!was_masked) && masked[index_first] == 1) {
									num_masked_in_Q++;
								} else if (was_masked && masked[index_first] == 0) {
									num_masked_in_Q--;
								}
							}
						}
					}
				}
			}
			for (int i = 0; i < G.num_nodes(); i++) {
				if (masked[i] == 0 && dist[i] > max_val && dist[i] < Double.MAX_VALUE) {
					max_val = dist[i];
				}
			}
		}

		return max_val;
	}

	/**
	 * Computes an all pairs shortest path array given a graph with positive edge
	 * weights
	 * 
	 * @param G    graph
	 * @param dist 2-dimensional distance array
	 * @param pre  2-dimensional predecessor arra
	 */
	public static void all_pairs_shortest_paths(Graph G, double dist[][], int pre[][]) {
		// for each node ...
		for (Node v : G.get_nodes().values()) {
			int key = v.reduced_key;
			// define comparator
			Comparator<Tupel<Double, Node>> comp = new Comparator<Tupel<Double, Node>>() {
				public int compare(Tupel<Double, Node> n_1, Tupel<Double, Node> n_2) {
					if (dist[key][n_2.second.reduced_key] - dist[key][n_1.second.reduced_key] > 0) {
						return -1;
					} else if (dist[key][n_2.second.reduced_key] - dist[key][n_1.second.reduced_key] < 0) {
						return 1;
					} else
						return 0;
				}
			};
			// initialize priority queue
			PriorityQueue<Tupel<Double, Node>> Q = new PriorityQueue<Tupel<Double, Node>>(G.num_nodes(), comp);
			Q.add(new Tupel<Double, Node>(dist[key][key], v));
			// until queue is empty ...
			while (Q.size() > 0) {
				// pop first element from queue
				Tupel<Double, Node> t = Q.remove();
				Node n = t.second;
				// if shortest path to node is not optimal -> skip
				if (dist[key][n.reduced_key] < t.first) {
					continue;
				}
				// forall incident edges ...
				for (Edge e : n.incident_edges.values()) {
					// get indices
					int index_first = G.get_nodes().get(e.first).reduced_key;
					int index_second = G.get_nodes().get(e.second).reduced_key;
					// get distance over edge
					double alt = dist[key][n.reduced_key] + e.length;
					// check if distance of n unvisited endpoint is smaller or larger
					if (e.first == n.key) {
						if (dist[key][index_second] > alt) {
							// update distance array
							dist[key][index_second] = alt;
							// update predecessor array
							pre[key][index_second] = n.reduced_key;
							// push updated end vertex to priority queue
							Q.add(new Tupel<Double, Node>(dist[key][index_second], G.get_nodes().get(e.second)));
						}
					} else if (e.second == n.key) {
						if (dist[key][index_first] > alt) {
							// update distance array
							dist[key][index_first] = alt;
							// update predecessor array
							pre[key][index_first] = n.reduced_key;
							// push updated end vertex to priority queue
							Q.add(new Tupel<Double, Node>(dist[key][index_first], G.get_nodes().get(e.first)));
						}
					}
				}
			}
		}
		return;
	}

	static void run_dijkstra(Graph G, double[] distance) {

		MHPriorityQueue Q = new MHPriorityQueue(G.num_nodes());
		for (Node n : G.get_nodes().values()) {
			if (distance[n.reduced_key] == 0.0) {
				Q.insert(n, distance[n.reduced_key]);
			}
		}

		while (Q.size() > 0) {
			Tupel<Double, Node> t = Q.extract_min();
			Node n = t.second;
			int key = n.reduced_key;

			if (distance[key] < t.first) {
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.get_nodes().get(e.first).reduced_key;
				int index_second = G.get_nodes().get(e.second).reduced_key;

				double alt = distance[key] + e.length;

				if (e.first == key) {
					if (distance[index_second] > alt) {
						Q.insert(G.get_nodes().get(e.second), distance[index_second] = alt);
					}
				} else if (e.second == key) {
					if (distance[index_first] > alt) {
						Q.insert(G.get_nodes().get(e.first), distance[index_first] = alt);
					}
				}
			}
		}
	}

	static void run_dijkstra_with_cutoff(Graph G, double[] distance, double max_dist) {

		MHPriorityQueue Q = new MHPriorityQueue(G.num_nodes());
		for (Node n : G.get_nodes().values()) {
			if (distance[n.reduced_key] == 0.0) {
				Q.insert(n, distance[n.reduced_key]);
			}
		}

		while (Q.size() > 0) {
			Tupel<Double, Node> t = Q.extract_min();
			Node n = t.second;
			int key = n.reduced_key;
			if (distance[key] > max_dist) {
				return;
			}
			if (distance[key] < t.first) {
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.get_nodes().get(e.first).reduced_key;
				int index_second = G.get_nodes().get(e.second).reduced_key;
				double alt = distance[key] + e.length;

				if (e.first == key) {
					if (distance[index_second] > alt) {
						Q.insert(G.get_nodes().get(e.second), distance[index_second] = alt);
					}
				} else if (e.second == key) {
					if (distance[index_first] > alt) {
						Q.insert(G.get_nodes().get(e.first), distance[index_first] = alt);
					}
				}
			}
		}
	}

	/*
	 * static void run_dijkstra(Graph G, double[] distance) {
	 * Comparator<Tupel<Double, Node>> comp = new Comparator<Tupel<Double, Node>>()
	 * { public int compare(Tupel<Double, Node> n_1, Tupel<Double, Node> n_2) { if
	 * (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] > 0) {
	 * return -1; } else if (distance[n_2.second.reduced_key] -
	 * distance[n_1.second.reduced_key] < 0) { return 1; } else return 0; } };
	 * 
	 * PriorityQueue<Tupel<Double, Node>> Q = new PriorityQueue<Tupel<Double,
	 * Node>>(G.num_nodes(), comp); for (Node n : G.get_nodes().values()) {
	 * Q.add(new Tupel<Double, Node>(distance[n.reduced_key], n)); }
	 * 
	 * while (Q.size() > 0) { Tupel<Double, Node> t = Q.remove(); Node n = t.second;
	 * if (distance[n.reduced_key] < t.first) { continue; } for (Edge e :
	 * n.incident_edges.values()) { int index_first =
	 * G.get_nodes().get(e.first).reduced_key; int index_second =
	 * G.get_nodes().get(e.second).reduced_key; double alt = distance[n.reduced_key]
	 * + e.length;
	 * 
	 * if (e.first == n.key) { if (distance[index_second] > alt) {
	 * distance[index_second] = alt; Q.add(new Tupel<Double,
	 * Node>(distance[index_second], G.get_nodes().get(e.second))); } } else if
	 * (e.second == n.key) { if (distance[index_first] > alt) {
	 * distance[index_first] = alt; Q.add(new Tupel<Double,
	 * Node>(distance[index_first], G.get_nodes().get(e.first))); } } else {
	 * System.out.println("alert"); } } } }
	 * 
	 * static void run_dijkstra(TwoDGraph G, double[] distance) {
	 * Comparator<Tupel<Double, TwoDNode>> comp = new Comparator<Tupel<Double,
	 * TwoDNode>>() { public int compare(Tupel<Double, TwoDNode> n_1, Tupel<Double,
	 * TwoDNode> n_2) { if (distance[n_2.second.reduced_key] -
	 * distance[n_1.second.reduced_key] > 0) { return -1; } else if
	 * (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] < 0) {
	 * return 1; } else return 0; } };
	 * 
	 * PriorityQueue<Tupel<Double, TwoDNode>> Q = new PriorityQueue<Tupel<Double,
	 * TwoDNode>>(G.num_nodes(), comp); for (TwoDNode n : G.get_nodes().values()) {
	 * Q.add(new Tupel<Double, TwoDNode>(distance[n.reduced_key], n)); }
	 * 
	 * while (Q.size() > 0) { Tupel<Double, TwoDNode> t = Q.remove(); TwoDNode n =
	 * t.second; if (distance[n.reduced_key] < t.first) { // here lies the
	 * problem!!! or does it? continue; } for (Edge e : n.incident_edges.values()) {
	 * int index_first = G.get_nodes().get(e.first).reduced_key; int index_second =
	 * G.get_nodes().get(e.second).reduced_key; double alt = distance[n.reduced_key]
	 * + e.length;
	 * 
	 * if (e.first == n.key) { if (distance[index_second] > alt) {
	 * distance[index_second] = alt; Q.add(new Tupel<Double,
	 * TwoDNode>(distance[index_second], G.get_nodes().get(e.second))); } } else if
	 * (e.second == n.key) { if (distance[index_first] > alt) {
	 * distance[index_first] = alt; Q.add(new Tupel<Double,
	 * TwoDNode>(distance[index_first], G.get_nodes().get(e.first))); } } else {
	 * System.out.println("alert"); } } } }
	 */
	static void run_dijkstra(TwoDGraph G, double[] distance) {
		MHPriorityQueue Q = new MHPriorityQueue(G.num_nodes());
		for (TwoDNode n : G.get_nodes().values()) {
			if (distance[n.reduced_key] == 0.0) {
				Q.insert(n, distance[n.reduced_key]);
			}
		}

		while (Q.size() > 0) {
			Tupel<Double, Node> t = Q.extract_min();
			Node n = t.second;
			int key = n.reduced_key;
			if (distance[key] < t.first) { // here lies the problem!!! or does it?
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.get_nodes().get(e.first).reduced_key;
				int index_second = G.get_nodes().get(e.second).reduced_key;
				double alt = distance[key] + e.length;

				if (e.first == key) {
					if (distance[index_second] > alt) {
						Q.insert(G.get_nodes().get(e.second), distance[index_second] = alt);
					}
				} else if (e.second == key) {
					if (distance[index_first] > alt) {
						Q.insert(G.get_nodes().get(e.first), distance[index_first] = alt);
					}
				}
			}
		}
	}

	/**
	 * Computes a distance array from a single source to all other vertices in a
	 * graph
	 * 
	 * @param S source
	 * 
	 * @param G graph G
	 * 
	 * @return distance array
	 */
	public static double[] dijkstra_single_source(Node S, Graph G) {
		double[] distance = new double[G.num_nodes()];
		for (Node n : G.get_nodes().values()) {
			distance[n.key] = Double.MAX_VALUE;
		}
		distance[S.reduced_key] = 0;
		run_dijkstra(G, distance);
		return distance;
	}

	/**
	 * Computes a distance array from a single source to all other vertices in a
	 * graph
	 * 
	 * @param S source
	 * @param G graph G
	 * @return distance array
	 */
	public static double[] dijkstra_single_source_update(Node S, double[] distance, Graph G) {
		// initialisiere
		distance[S.reduced_key] = 0;
		run_dijkstra(G, distance);
		return distance;
	}

	static double[] dijkstra_gap_eval_rm_upd_subgraph(boolean[] selection, Graph G, double[] distance,
			Edge removed_edge, double nuke_radius) {
		TwoDNode r_1 = (TwoDNode) G.get_nodes().get(removed_edge.first);
		TwoDNode r_2 = (TwoDNode) G.get_nodes().get(removed_edge.second);
		// create aubgraph
		TwoDGraph S = new TwoDGraph();
		double[] correct = Graph_Eval.dijkstra_gap_eval(selection, G);
		// add nodes to subgraph
		for (Node n : G.get_nodes().values()) {
			TwoDNode n_2d = (TwoDNode) n;
			double d_1 = 0.0d;
			double d_2 = 0.0d;
			if (G.TWO_D_TYPE == 1) {
				d_1 = Utility.dist_lat_lon(r_1.y, r_1.x, n_2d.y, n_2d.x);
				d_2 = Utility.dist_lat_lon(r_2.y, r_2.x, n_2d.y, n_2d.x);
			} else if (G.TWO_D_TYPE == 2) {
				d_1 = Math.sqrt((r_1.x - n_2d.x) * (r_1.x - n_2d.x) + (r_1.y - n_2d.y) * (r_1.y - n_2d.y));
				d_2 = Math.sqrt((r_2.x - n_2d.x) * (r_2.x - n_2d.x) + (r_2.y - n_2d.y) * (r_2.y - n_2d.y));
			} else {
				return null;
			}
			if (Math.min(d_1, d_2) <= nuke_radius) {
				S.add_node(new TwoDNode(n_2d.reduced_key, n_2d.x, n_2d.y));
				if (n_2d.degree() > 1) {
					distance[n.key] = Double.MAX_VALUE;
				}
			}
		}
		// add edges to subgraph
		ArrayList<Edge> add_edges = new ArrayList<Edge>();
		ArrayList<TwoDNode> add_nodes = new ArrayList<TwoDNode>();
		boolean indicator = false;
		for (Edge e : G.get_edges()) {
			// fully contained edges
			if (S.get_nodes().containsKey(e.first) && S.get_nodes().containsKey(e.second)) {
				add_edges.add(e);
			} else if (S.get_nodes().containsKey(e.first)) {
				// edges with second node outside
				TwoDNode second = (TwoDNode) G.get_nodes().get(e.second);
				add_nodes.add(new TwoDNode(second.reduced_key, second.x, second.y));
				add_edges.add(e);
			} else if (S.get_nodes().containsKey(e.second)) {
				// edges with first node outside
				TwoDNode first = (TwoDNode) G.get_nodes().get(e.first);
				add_nodes.add(new TwoDNode(first.reduced_key, first.x, first.y));
				add_edges.add(e);
			}
			// set selected edge endpoints as sources
			if (selection[e.index]) {
				if (e.index == 67 || e.index == 68) {
					System.out.println("here");
					indicator = true;
				}
				distance[e.first] = 0.0;
				distance[e.second] = 0.0;
			}
		}
		for (TwoDNode n : add_nodes) {
			S.add_node(n);
		}
		for (Edge e : add_edges) {
			S.add_edge(e);
		}

		// System.out.println(S.num_nodes() + " " + S.num_edges());
		run_dijkstra(S, distance);
		boolean exit = false;
		for (int i = 0; i < distance.length; i++) {
			if (distance[i] != correct[i]) {
				System.err.println(
						"Error: distance[" + i + "] contains value " + distance[i] + " != correct " + correct[i]);
				System.err.println("G.|V| = " + G.num_nodes() + ", S.|V| = " + S.num_nodes());
				TwoDNode n = (TwoDNode) G.get_nodes().get(i);
				System.err.println("Affected node:\n" + n.toString());
				System.err.println("Node contained in S: " + S.get_nodes().containsKey(i));
				double d_1 = 0.0d;
				double d_2 = 0.0d;
				if (G.TWO_D_TYPE == 1) {
					d_1 = Utility.dist_lat_lon(r_1.y, r_1.x, n.y, n.x);
					d_2 = Utility.dist_lat_lon(r_2.y, r_2.x, n.y, n.x);
				} else if (G.TWO_D_TYPE == 2) {
					d_1 = Math.sqrt((r_1.x - n.x) * (r_1.x - n.x) + (r_1.y - n.y) * (r_1.y - n.y));
					d_2 = Math.sqrt((r_2.x - n.x) * (r_2.x - n.x) + (r_2.y - n.y) * (r_2.y - n.y));
				} else {
					return null;
				}
				System.err.println("Distance from removed edge: " + Math.min(d_1, d_2) + " vs radius "
						+ (nuke_radius + removed_edge.length));
				System.err.println("Incident edges:");
				for (Edge e : n.incident_edges.values()) {
					System.err.println(e.index + " -> " + e.toString());
				}
				if (indicator) {
					System.err.println("edge was marked!");
				}
				System.err.println();
				exit = true;
			}

		}
		if (exit) {
			for (Edge e : G.get_edges()) {
				if (selection[e.index]) {
					System.err.print(e.index + ", ");
				}
			}
			System.err.println();
			System.err.println("Removed edge was: " + removed_edge.toString());
			System.exit(1);
		}
		return distance;
	}

	// TODO fix correctness
	@SuppressWarnings("unused")
	private static double[] dijkstra_gap_eval_rm_upd(boolean[] selection, Graph G, double[] distance, Edge removed_edge,
			double nuke_radius) {
		TwoDNode r_1 = (TwoDNode) G.get_nodes().get(removed_edge.first);
		TwoDNode r_2 = (TwoDNode) G.get_nodes().get(removed_edge.second);
		for (Node n : G.get_nodes().values()) {
			TwoDNode n_2d = (TwoDNode) n;
			double d_1 = 0.0d;
			double d_2 = 0.0d;
			if (G.TWO_D_TYPE == 1) {
				d_1 = Utility.dist_lat_lon(r_1.y, r_1.x, n_2d.y, n_2d.x);
				d_2 = Utility.dist_lat_lon(r_2.y, r_2.x, n_2d.y, n_2d.x);
			} else if (G.TWO_D_TYPE == 2) {
				d_1 = Math.sqrt((r_1.x - n_2d.x) * (r_1.x - n_2d.x) + (r_1.y - n_2d.y) * (r_1.y - n_2d.y));
				d_2 = Math.sqrt((r_2.x - n_2d.x) * (r_2.x - n_2d.x) + (r_2.y - n_2d.y) * (r_2.y - n_2d.y));
			}
			if (Math.min(d_1, d_2) <= nuke_radius) {
				distance[n.reduced_key] = Double.MAX_VALUE;
			}
		}
		for (Edge e : G.get_edges()) {
			if (selection[e.index]) {
				distance[e.first] = 0.0;
				distance[e.second] = 0.0;
			}
		}
		run_dijkstra(G, distance);
		return distance;
	}

	public static double[] dijkstra_two_source_removal_upd(boolean[] selection, Graph G, double[] distance,
			Edge removed_edge, double nuke_radius) {
		if (!(G.get_nodes().entrySet().iterator().next().getValue() instanceof TwoDNode)) {
			return Graph_Eval.dijkstra_gap_eval(selection, G);
		} else {
			return dijkstra_gap_eval_rm_upd_subgraph(selection, G, distance, removed_edge, nuke_radius);
		}
	}

	public static double[] dijkstra_default(boolean[] selection, Graph G) {
		double[] distance = new double[G.num_nodes()];
		for (Node n : G.get_nodes().values()) {
			distance[n.reduced_key] = Double.MAX_VALUE;
		}
		for (int i = 0; i < selection.length; i++) {
			if (selection[i]) {
				Edge e = G.get_edges().get(i);
				distance[G.get_nodes().get(e.first).reduced_key] = 0;
				distance[G.get_nodes().get(e.second).reduced_key] = 0;
			}
		}
		run_dijkstra(G, distance);
		return distance;
	}

	/**
	 * Computes a distance array from a given selection of edges in a graph G
	 * 
	 * @param selection selection of edges
	 * @param G         graph G
	 * @return distance array
	 */
	public static double[] dijkstra_gap_eval(boolean[] selection, Graph G) {
		double[] distance = new double[G.num_nodes()];
		for (Node n : G.get_nodes().values()) {
			distance[n.reduced_key] = Double.MAX_VALUE;
		}
		for (int i = 0; i < selection.length; i++) {
			if (selection[i]) {
				Edge e = G.get_edges().get(i);
				distance[G.get_nodes().get(e.first).reduced_key] = 0;
				distance[G.get_nodes().get(e.second).reduced_key] = 0;
			}
		}
		// System.out.println(Arrays.toString(distance));
		run_dijkstra(G, distance);
		return distance;

	}

	/**
	 * Computes a distance array for a graph G based on a pre-calculated distance
	 * array and two sources.
	 * 
	 * @param S_1      source 1
	 * @param S_2      source 2
	 * @param distance distance array
	 * @param G        graph G
	 * @return distance array computed using Dijkstra's algorithm
	 */
	public static double[] dijkstra_two_source_update(Node S_1, Node S_2, double[] distance, Graph G) {
		distance[S_1.reduced_key] = 0;
		distance[S_2.reduced_key] = 0;

		run_dijkstra(G, distance);
		return distance;
	}

	public static Tupel<Edge, Double> evaluate_selection_farthest_midpoint(boolean[] sel, Graph G, double[] dist) {
		double max_gap = 0.0;
		Edge max_edge = null;
		for (Edge e : G.get_edges()) {
			if (!sel[e.index]) {
				double d_e = Math.min(dist[e.first], dist[e.second]) + 0.5 * e.length;
				if (max_gap < d_e || max_gap == 0.0) {
					max_edge = e;
					max_gap = d_e;
				}
			}
		}
		return new Tupel<Edge, Double>(max_edge, max_gap);
	}

	public static Tupel<Edge, Double> evaluate_selection_farthest_midpoint(boolean[] sel, Graph G) {
		double[] dist = new double[G.num_nodes()];
		for (int i = 0; i < dist.length; i++) {
			dist[i] = Double.MAX_VALUE;
		}
		for (int i = 0; i < sel.length; i++) {
			if (sel[i]) {
				dist[G.get_edges().get(i).first] = 0.0;
				dist[G.get_edges().get(i).second] = 0.0;
			}
		}
		run_dijkstra(G, dist);
		double max_gap = 0.0;
		Edge max_edge = null;
		for (Edge e : G.get_edges()) {
			if (!sel[e.index]) {
				double d_e = Math.min(dist[e.first], dist[e.second]) + 0.5 * e.length;
				if (max_gap < d_e || max_gap == 0.0) {
					max_edge = e;
					max_gap = d_e;
				}
			}
		}
		return new Tupel<Edge, Double>(max_edge, max_gap);
	}

	/**
	 * Function to evaluate a selection boolean array based on the gap defined by
	 * d(u) + |(u,v)| + d(v) on a given distance array
	 * 
	 * @param sel  selection
	 * @param G    graph G
	 * @param dist distance array
	 * @return Tupel of edge with maximum gap and the length of the max gap
	 */
	public static Tupel<Edge, Double> evaluate_selection(boolean[] sel, Graph G, double[] dist) {
		double max_gap = 0.0;
		Edge max_edge = null;
		for (Edge e : G.get_edges()) {
			if (!sel[e.index]) {
				double d_e = dist[e.first] + dist[e.second] + e.length;
				if (max_gap < d_e || max_gap == 0.0) {
					max_edge = e;
					max_gap = d_e;
				}
			}
		}
		return new Tupel<Edge, Double>(max_edge, max_gap);
	}

	/**
	 * Function to evaluate a selection boolean array based on the gap defined by
	 * d(u) + |(u,v)| + d(v)
	 * 
	 * @param sel selection to evaluate
	 * @param G   graph G
	 * @return Tupel of edge with maximum gap and the length of the max gap
	 */
	public static Tupel<Edge, Double> evaluate_selection(boolean[] sel, Graph G) {
		double[] dist = new double[G.num_nodes()];
		dist = dijkstra_gap_eval(sel, G);
		double max_gap = 0.0;
		Edge max_edge = null;
		for (Edge e : G.get_edges()) {
			if (!sel[e.index]) {
				double d_e = dist[e.first] + dist[e.second] + e.length;
				if (max_gap < d_e) {
					max_edge = e;
					max_gap = d_e;
				}
			}
		}
		return new Tupel<Edge, Double>(max_edge, max_gap);
	}

}
