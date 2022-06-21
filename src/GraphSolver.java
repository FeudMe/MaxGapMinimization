import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import me.tongfei.progressbar.ProgressBar;

/**
 * Algorithms for solving Maximum Gap Minimization on graphs
 * 
 * @author Oliver Leenders
 *
 */
public abstract class GraphSolver {
	/**
	 * Private constructor
	 */
	private GraphSolver() {
	}

	/**
	 * Distance based sequential addition algorithm
	 * 
	 * @param G      graph G
	 * @param k      number of edges to be selected
	 * @param notify whether to print a status message every 100 edges
	 * @return selection
	 */
	public static boolean[] sequential_addition_dist(Graph G, int k, boolean notify) {
		boolean[] sel = new boolean[G.num_edges()];

		Random rd = new Random(System.currentTimeMillis());
		Edge max_len = G.get_edges().get(rd.nextInt(G.num_edges()));
		sel[max_len.index] = true;

		double[] dist = Graph_Eval.dijkstra_default(sel, G);

		int prev_sel_1 = max_len.first;
		int prev_sel_2 = max_len.second;

		if (notify) {
			// display progress bar (make sure, eclipse console interprets ASCII control
			// characters -- especially '\r' -- correctly) and that console encoding is set
			// to utf-8
			try {
				ProgressBar pb = new ProgressBar("Distance Minimizer", k);
				for (int i = 1; i < k; i++) {
					dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
							G.get_nodes().get(prev_sel_2), dist, G);
					int max_index = Utility.index_of_max(dist);

					max_len = new Edge(-1, -1, -1, -1);

					for (Edge e_1 : G.get_nodes().get(max_index).incident_edges.values()) {
						if (e_1.length > max_len.length) {
							max_len = e_1;
						}
					}

					prev_sel_1 = max_len.first;
					prev_sel_2 = max_len.second;
					sel[max_len.index] = true;
					pb.step();
				}
				pb.close();

			} catch (Exception e) {
				e.printStackTrace();
			}
		} else {

			for (int i = 1; i < k; i++) {
				dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
						G.get_nodes().get(prev_sel_2), dist, G);
				int max_index = Utility.index_of_max(dist);

				max_len = new Edge(-1, -1, -1, -1);

				for (Edge e_1 : G.get_nodes().get(max_index).incident_edges.values()) {
					if (e_1.length > max_len.length) {
						max_len = e_1;
					}
				}

				prev_sel_1 = max_len.first;
				prev_sel_2 = max_len.second;
				sel[max_len.index] = true;
			}
		}
		return sel;
	}

	public static boolean[] k_seq_add_dist(boolean[] sel, Graph G, int k) {
		int i = 1;

		double[] dist = Graph_Eval.dijkstra_default(sel, G);
		int max_index = Utility.index_of_max(dist);

		Edge max_len = new Edge(-1, -1, -1, -1);

		for (Edge e_1 : G.get_nodes().get(max_index).incident_edges.values()) {
			if (e_1.length > max_len.length) {
				max_len = e_1;
			}
		}

		int prev_sel_1 = max_len.first;
		int prev_sel_2 = max_len.second;
		sel[max_len.index] = true;

		while (i < k) {
			dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1), G.get_nodes().get(prev_sel_2),
					dist, G);
			max_index = Utility.index_of_max(dist);

			max_len = new Edge(-1, -1, -1, -1);

			for (Edge e_1 : G.get_nodes().get(max_index).incident_edges.values()) {
				if (e_1.length > max_len.length) {
					max_len = e_1;
				}
			}

			prev_sel_1 = max_len.first;
			prev_sel_2 = max_len.second;
			sel[max_len.index] = true;

			i++;
		}
		return sel;
	}
	
	

	/**
	 * Gap-based sequential addition algorithm
	 * 
	 * @param G      graph G
	 * @param k      number of edges to be selected
	 * @param notify whether to print a status message every 100 edges
	 * @return selection
	 */
	public static boolean[] sequential_addition_gap(Graph G, int k, boolean notify) {
		boolean[] sel = new boolean[G.num_edges()];
		Random rd = new Random();
		
		Edge start_edge = G.get_edges().get(rd.nextInt(G.num_edges()));
		sel[start_edge.index] = true;
		double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
		int prev_sel_1 = start_edge.first;
		int prev_sel_2 = start_edge.second;
		if (notify) {
			try {
				ProgressBar pb = new ProgressBar("Maximum Gap", k);
				for (int i = 1; i < k; i++) {
					dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
							G.get_nodes().get(prev_sel_2), dist, G);
					Edge max_len = Graph_Eval.evaluate_selection(sel, G, dist).first;
					prev_sel_1 = max_len.first;
					prev_sel_2 = max_len.second;
					sel[max_len.index] = true;
					pb.step();
				}
				pb.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else {
			for (int i = 1; i < k; i++) {
				dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
						G.get_nodes().get(prev_sel_2), dist, G);
				Edge max_len = Graph_Eval.evaluate_selection(sel, G, dist).first;
				prev_sel_1 = max_len.first;
				prev_sel_2 = max_len.second;
				sel[max_len.index] = true;
			}
		}
		return sel;
	}

	public static boolean[] sequential_addition_farthest_midpoint(Graph G, int k, boolean notify) {
		boolean[] sel = new boolean[G.num_edges()];
		Random rd = new Random();
		
		Edge start_edge = G.get_edges().get(rd.nextInt(G.num_edges()));
		sel[start_edge.index] = true;
		double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
		int prev_sel_1 = start_edge.first;
		int prev_sel_2 = start_edge.second;
		if (notify) {
			try {
				ProgressBar pb = new ProgressBar("Maximum Gap", k);
				for (int i = 1; i < k; i++) {
					dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
							G.get_nodes().get(prev_sel_2), dist, G);
					Edge max_len = Graph_Eval.evaluate_selection_farthest_midpoint(sel, G, dist).first;
					prev_sel_1 = max_len.first;
					prev_sel_2 = max_len.second;
					sel[max_len.index] = true;
					pb.step();
				}
				pb.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else {
			for (int i = 1; i < k; i++) {
				dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1),
						G.get_nodes().get(prev_sel_2), dist, G);
				Edge max_len = Graph_Eval.evaluate_selection_farthest_midpoint(sel, G, dist).first;
				prev_sel_1 = max_len.first;
				prev_sel_2 = max_len.second;
				sel[max_len.index] = true;
			}
		}
		return sel;
	}

	
	/**
	 * Computes a heuristic for MGMP. Performs a parallel variant of the two-source
	 * shortest path problem.
	 * 
	 * @param G      Graph G
	 * @param k      number of vertices to select
	 * @param notify whether to print status reports
	 * @return selection
	 */
	public static boolean[] sequential_addition_gap_p(Graph G, int k, boolean notify) {
		int i = 1;
		boolean[] sel = new boolean[G.num_edges()];
		Random rd = new Random(System.currentTimeMillis());
		Edge start_edge = G.get_edges().get(rd.nextInt(G.num_edges()));
		sel[start_edge.index] = true;
		double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
		int prev_sel_1 = start_edge.first;
		int prev_sel_2 = start_edge.second;
		while (i < k) {
			if (i % 100 == 0 && notify) {
				System.out.println(i + "/" + k + " edges selected -- " + (k - i) + " to go");
			}
			double[] dist_1 = dist.clone();
			double[] dist_2 = dist.clone();
			final int index = prev_sel_1;

			Thread t = new Thread(() -> {
				Graph_Eval.dijkstra_single_source_update(G.get_nodes().get(index), dist_1, G);
			});
			t.start();
			dist_2 = Graph_Eval.dijkstra_single_source_update(G.get_nodes().get(prev_sel_2), dist_2, G);
			try {
				t.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			dist = Utility.min(dist_1, dist_2);
			Edge max_len = Graph_Eval.evaluate_selection(sel, G, dist).first;
			prev_sel_1 = max_len.first;
			prev_sel_2 = max_len.second;
			sel[max_len.index] = true;

			i++;
		}
		return sel;
	}

	public static boolean[] one_opt(Graph G, boolean[] sel) {
		boolean has_changed = true;
		// System.out.println(prev_max_gap);
		double curr_max_gap_length = Graph_Eval.evaluate_selection(sel, G).second;
		while (has_changed) {
			has_changed = false;
			double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
			for (int i = 0; i < sel.length; i++) {
				if (sel[i]) {
					sel[i] = false;
					Edge e = G.get_edges().get(i);
					dist = Graph_Eval.dijkstra_two_source_removal_upd(sel, G, dist, e, curr_max_gap_length);
					Tupel<Edge, Double> t_0 = Graph_Eval.evaluate_selection(sel, G, dist);
					int j = t_0.first.index;
					if (j != i) {
						sel[j] = true;
						dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(G.get_edges().get(j).first),
								G.get_nodes().get(G.get_edges().get(j).second), dist, G);
						double new_max_gap_length = Graph_Eval.evaluate_selection(sel, G, dist).second;
						if (new_max_gap_length < curr_max_gap_length) {
							has_changed = true;
							curr_max_gap_length = new_max_gap_length;
						} else {
							sel[j] = false;
							sel[i] = true;
						}
					} else {
						sel[i] = true;
					}
				}
			}
		}
		return sel;
	}
	

	public static boolean[] one_opt_naive(Graph G, boolean[] sel) {
		boolean has_changed = true;
		// System.out.println(prev_max_gap);
		double curr_max_gap_length = Graph_Eval.evaluate_selection(sel, G).second;
		while (has_changed) {
			has_changed = false;
			double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
			for (int i = 0; i < sel.length; i++) {
				if (sel[i]) {
					sel[i] = false;
					dist = Graph_Eval.dijkstra_gap_eval(sel, G);
					Tupel<Edge, Double> t_0 = Graph_Eval.evaluate_selection(sel, G, dist);
					int j = t_0.first.index;
					if (j != i) {
						sel[j] = true;
						dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(G.get_edges().get(j).first),
								G.get_nodes().get(G.get_edges().get(j).second), dist, G);
						double new_max_gap_length = Graph_Eval.evaluate_selection(sel, G, dist).second;
						if (new_max_gap_length < curr_max_gap_length) {
							has_changed = true;
							curr_max_gap_length = new_max_gap_length;
							// System.out.println(new_max_gap_length);
						} else {
							sel[j] = false;
							sel[i] = true;
						}
					} else {
						sel[i] = true;
					}
				}
			}
		}
		return sel;
	}

	public static boolean[] k_seq_add_gap(boolean[] sel, Graph G, int k) {
		int i = 1;
		double[] dist = Graph_Eval.dijkstra_gap_eval(sel, G);
		Edge start_edge = Graph_Eval.evaluate_selection(sel, G, dist).first;
		sel[start_edge.index] = true;
		int prev_sel_1 = start_edge.first;
		int prev_sel_2 = start_edge.second;

		while (i < k) {
			dist = Graph_Eval.dijkstra_two_source_update(G.get_nodes().get(prev_sel_1), G.get_nodes().get(prev_sel_2),
					dist, G);
			Edge max_len = Graph_Eval.evaluate_selection(sel, G, dist).first;
			prev_sel_1 = max_len.first;
			prev_sel_2 = max_len.second;
			sel[max_len.index] = true;

			i++;
		}
		return sel;
	}

	/**
	 * Greedy algorithm selecting the k longest edges
	 * 
	 * @param G graph G
	 * @param k number of edges to be selected
	 * @return selection
	 */
	public static boolean[] greedy_heuristic(Graph G, int k) {
		@SuppressWarnings("unchecked")
		ArrayList<Edge> sel_edges = (ArrayList<Edge>) G.get_edges().clone();
		Comparator<Edge> comp = new Comparator<Edge>() {
			public int compare(Edge e1, Edge e2) {
				if (e2.length > e1.length) {
					return 1;
				} else if (e1.length > e2.length) {
					return -1;
				}
				return 0;
			}
		};
		Collections.sort(sel_edges, comp);
		boolean[] sel = new boolean[G.num_edges()];
		for (int i = 0; i < k; i++) {
			Edge e = sel_edges.get(i);
			sel[e.index] = true;
		}
		return sel;
	}

}
