import java.io.File;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

import me.tongfei.progressbar.ProgressBar;

/**
 * Driver Code for Maximum Gap Minimization testing
 * 
 * @author Oliver Leenders
 *
 */
public class Main {
	public static void main(String[] args) {
		// lsp_eval_test();
		// eval_runtime_test();
		test_tree();
		// runtime_test_polyline();
		// quality_test_polyline();
		// quality_test_graph();
		// graph_test();
		// test_heap();
	}

	private static void test_tree() {
//		Tree T = Graph_Util.read_tree(new File("src/Benchmark_Graphs/test.tree"));
//		System.out.println(T.toString());
//		boolean[] S = new boolean[T.edges.size()];
//		double[] d = new double[T.nodes.size()];
//		S = TreeSolver.parametric_search(T, 8);
//		System.out.println(Utility.num_selected(S));
//		for (int i = 0; i < S.length; i++) {
//			if (S[i]) {
//				System.out.println(T.edges.get(i).toString());
//			}
//		}
//		System.out.println("LSP-length: " + TreeSolver.eval_LSP(T, S));
//		boolean[] S_greedy = new boolean[T.edges.size()];
//		TreeSolver.select_greedy(T, S_greedy, 8);
//		System.out.println();
//		System.out.println(Utility.num_selected(S_greedy));
//		for (int i = 0; i < S_greedy.length; i++) {
//			if (S_greedy[i]) {
//				System.out.println(T.edges.get(i).toString());
//			}
//		}
//		System.out.println("LSP-length: " + TreeSolver.eval_LSP(T, S_greedy));
//		
		Graph G = Graph_Util.read_osm_graph(new File("src/Benchmark_Graphs/larger_map.osm"));

		G = Graph_Util.largest_connected_component(G);
		G.reduce_graph();

		Tree T = new Tree(G);
		// Tree T_2 = Graph_Util.read_tree(new File("src/Benchmark_Graphs/test.tree"));

		String csv = "k, random, k-largest, greedy, exact\n";

		int stepSize = T.edges.size() / 50;
		for (int k = stepSize; k < T.edges.size(); k += stepSize) {
			System.out.print(k / stepSize + ", ");
			double cutoff = Double.MAX_VALUE;

			csv += k + ", ";

			double q_rd = 0.0;
			for (int i = 0; i < 10; i++) {
				boolean[] random = Utility.random_selection(T.edges.size(), k);
				double random_d = TreeSolver.eval_tree_DFS(T, random, 0, 0, new boolean[T.nodes.size()]).second;
				cutoff = Math.min(random_d, cutoff);
				q_rd += (random_d) / 10;
			}
			csv += q_rd + ", ";

			boolean[] k_largest = TreeSolver.k_largest(T, k);
			double k_largest_d = TreeSolver.eval_tree_DFS(T, k_largest, 0, 0, new boolean[T.nodes.size()]).second;
			cutoff = Math.min(k_largest_d, cutoff);
			csv += k_largest_d + ", ";

			boolean[] S_greedy = new boolean[T.edges.size()];
			TreeSolver.select_greedy(T, S_greedy, k);
			double greedy_d = TreeSolver.eval_tree_DFS(T, S_greedy, 0, 0, new boolean[T.nodes.size()]).second;
			cutoff = Math.min(cutoff, greedy_d);
			csv += greedy_d + ", ";

			boolean[] S_exact = TreeSolver.parametric_search(T, k, (int) (Math.random() * T.nodes.size()));

			csv += TreeSolver.eval_tree_DFS(T, S_exact, 0, 0, new boolean[T.nodes.size()]).second + "\n";
		}
		System.out.println();
		System.out.println();
		System.out.println(csv);

	}

	private static void eval_runtime_test() {
		String csv = "";

//		for (double i = 0.01; i < 0.2; i += 0.01) {
//			Graph G = Graph_Util.subgraph_100_k(i);
//			System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//			csv += G.num_nodes() + ", ";
//			int k = (int) Math.sqrt(G.num_edges());
//			boolean[] S = Utility.random_selection(G.num_edges(), k);
//			csv += test_ten_times(G, S);
//		}

//		Graph G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.NY.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		int k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));
//		
//		

//		Graph G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.FLA.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		int k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));

//		Graph G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.LKS.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		int k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));
//		

//		Graph G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.CAL.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		int k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));
//		
//		G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.NE.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));
//		
//		G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.NW.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));

		Graph G = Graph_Util.read_usa_graph(new File("src/Benchmark_Graphs/TX.tmp"));
		G = Graph_Util.largest_connected_component(G);
		G.reduce_graph();
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());

		int k = (int) Math.sqrt(G.num_edges());
		csv += G.num_nodes() + ", ";
		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));

//		Graph G = Graph_Util.read_dimacs_graph(new File("src/Benchmark_Graphs/USA-road-d.W.gr"));
//		G = Graph_Util.largest_connected_component(G);
//		G.reduce_graph();
//		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
//		
//		k = (int) Math.sqrt(G.num_edges());
//		csv += G.num_nodes() + ", ";
//		csv += test_ten_times(G, Utility.random_selection(G.num_edges(), k));
//		

		System.out.println(csv);
	}

	private static String test_ten_times(Graph G, boolean[] S) {
		long start;
		long duration = 0;
		String lcsv = "";
		// evaluate dist
		for (int j = 0; j < 10; j++) {
			start = System.nanoTime();
			Utility.max(Graph_Eval.dijkstra_default(S, G));
			duration += (System.nanoTime() - start) / 10;
		}
		lcsv += duration + ", ";
		// evaluate CG
		duration = 0;
		for (int j = 0; j < 10; j++) {
			start = System.nanoTime();
			Graph_Eval.evaluate_selection(S, G);
			duration += (System.nanoTime() - start) / 10;
		}
		lcsv += duration + ", ";
		// evaluate farthest midpoint
		duration = 0;
		for (int j = 0; j < 10; j++) {
			start = System.nanoTime();
			// Graph_Eval.evaluate_selection_farthest_midpoint(S, G);
			duration += (System.nanoTime() - start) / 10;
		}
		lcsv += duration + "\n";
		return lcsv;
	}

	private static void test_heap() {
		Node n_0 = new Node(0);
		Node n_1 = new Node(1);
		Node n_2 = new Node(2);
		Node n_3 = new Node(3);
		Node n_4 = new Node(4);
		Node n_5 = new Node(5);

		MHPriorityQueue Q = new MHPriorityQueue(6);

		Q.insert(n_0, 4.2);
		Q.insert(n_1, 2.2);
		Q.insert(n_2, 3.5);
		Q.insert(n_3, 4.9);
		Q.insert(n_4, 1.2);
		Q.insert(n_5, Double.MAX_VALUE);

		Q.print();

		Q.decrease_key(n_0, 0.0);

		Q.print();

		Q.extract_min();

		Q.print();

		return;
	}

	/**
	 * Function testing solution quality on graphs
	 */
	private static void quality_test_graph() {
		System.out.println("Starting the graph quality test ...");
		File f = new File("src/Benchmark_Graphs/larger_map.osm");
		Graph G = Graph_Util.read_osm_graph(f);
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());

		System.out.println("Reducing graph ...");

		G = Graph_Util.largest_connected_component(G);
		G.reduce_graph();

		// System.out.println(G.toString());
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
		System.out.println(Graph_Util.connected_components(G));

		String csv = "k, random gap, greedy gap, sequential addition gap, MGM gap, fm, 1-OPT gap\n";
		String csv_dist = "k, random, greedy, sequential addition, MGM dist, fm, 1-OPT\n";
		String csv_fm = "k, random, greedy, sequential addition, MGM dist, fm,  1-OPT\n";
		int stepsize = G.num_edges() / 25;
		boolean[] dist_sel = new boolean[G.num_edges()];
		boolean[] gap_sel = new boolean[G.num_edges()];
		boolean[] rd_sel = new boolean[G.num_edges()];
		boolean[] greedy_sel = new boolean[G.num_edges()];
		boolean[] one_opt_sel = new boolean[G.num_edges()];
		boolean[] fm_sel = new boolean[G.num_edges()];
		for (int i = stepsize; i < G.num_edges(); i += stepsize) {
			csv += i + ", ";
			csv_dist += i + ", ";
			csv_fm += i + ", ";
			System.out.println(i + "/" + G.num_edges() + " edges to be selected");
			double quality = 0.0;
			double quality_gap = 0.0;
			double quality_fm = 0.0;

			int rounding_runs = 25;
			for (int j = 0; j < rounding_runs; j++) {
				rd_sel = Utility.random_selection(G.num_edges(), i);
				quality += Utility.max(Graph_Eval.dijkstra_default(rd_sel, G)) / rounding_runs;
				quality_gap += Graph_Eval.evaluate_selection(rd_sel, G).second / rounding_runs;
				quality_fm += Graph_Eval.evaluate_selection_farthest_midpoint(rd_sel, G).second / rounding_runs;
			}
			csv_dist += quality + ", ";
			csv += quality_gap + ", ";
			csv_fm += quality_fm + ", ";

			greedy_sel = GraphSolver.greedy_heuristic(G, i);
			csv_dist += Utility.max(Graph_Eval.dijkstra_default(greedy_sel, G)) + ", ";
			csv += Graph_Eval.evaluate_selection(greedy_sel, G).second + ", ";
			csv_fm += Graph_Eval.evaluate_selection_farthest_midpoint(greedy_sel, G).second + ", ";

			dist_sel = GraphSolver.k_seq_add_dist(dist_sel, G, stepsize);
			csv_dist += Utility.max(Graph_Eval.dijkstra_default(dist_sel, G)) + ", ";
			csv += Graph_Eval.evaluate_selection(dist_sel, G).second + ", ";
			csv_fm += Graph_Eval.evaluate_selection_farthest_midpoint(dist_sel, G).second + ", ";

			gap_sel = GraphSolver.k_seq_add_gap(gap_sel, G, stepsize);
			csv_dist += Utility.max(Graph_Eval.dijkstra_default(gap_sel, G)) + ", ";
			csv += Graph_Eval.evaluate_selection(gap_sel, G).second + ", ";
			csv_fm += Graph_Eval.evaluate_selection_farthest_midpoint(gap_sel, G).second + ", ";

			fm_sel = GraphSolver.sequential_addition_farthest_midpoint(G, i, false);
			csv_dist += Utility.max(Graph_Eval.dijkstra_default(fm_sel, G)) + ", ";
			csv += Graph_Eval.evaluate_selection(fm_sel, G).second + ", ";
			csv_fm += Graph_Eval.evaluate_selection_farthest_midpoint(fm_sel, G).second + ", ";

			one_opt_sel = GraphSolver.one_opt_naive(G, gap_sel.clone());
			csv_dist += Utility.max(Graph_Eval.dijkstra_default(one_opt_sel, G)) + "\n";
			csv += Graph_Eval.evaluate_selection(one_opt_sel, G).second + "\n";
			csv_fm += Graph_Eval.evaluate_selection_farthest_midpoint(one_opt_sel, G).second + "\n";

		}
		System.out.println();
		System.out.println(csv);
		System.out.println();
		System.out.println(csv_dist);
		System.out.println();
		System.out.println(csv_fm);
		// System.out.println(csv_dist);
	}

	public static void lsp_eval_test() {
		File f = new File("src/Benchmark_Graphs/Rome.txt");
		Graph G = Graph_Util.read_dimacs_graph(f);
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());

		System.out.println("Reducing graph ...");
		String csv = "";
		G = Graph_Util.largest_connected_component(G);
		G.reduce_graph();
		try {
			ProgressBar pb = new ProgressBar("Progress", G.num_edges());
			int stepsize = 100;
			for (int k = stepsize; k < G.num_edges(); k += stepsize) {
				// System.out.println(G.toString());
				// System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
				csv += k + ", ";
				boolean[] S = Utility.random_selection(G.num_edges(), k);

				double[] distance = Graph_Eval.dijkstra_default(S, G);
				double max_dist = Utility.max(distance);
				// System.out.println("Max Dist: " + max_dist);
				double max_composite_gap = Graph_Eval.evaluate_selection(S, G, distance).second;
				// System.out.println("Max Composite Gap: " + max_composite_gap);
				double farthest_midpoint_dist = Graph_Eval.evaluate_selection_farthest_midpoint(S, G, distance).second;
				// System.out.println("Farthest_Midpoint: " + farthest_midpoint_dist);
				long start;
				long duration_1 = 0;
				long duration_2 = 0;
				for (int i = 0; i < 2; i++) {
					S = Utility.random_selection(G.num_edges(), k);
					start = System.nanoTime();
					double LSP = Graph_Eval.evaluate_LSP(G, S);
					duration_1 += System.nanoTime() - start;

					// System.out.println("LSP: " + LSP + " in " + duration + " ns");

					start = System.nanoTime();
					double LSP_optimised = Graph_Eval.evaluate_LSP_memory_optimised(G, S);
					duration_2 = System.nanoTime() - start;
					// System.out.println("LSP optimised: " + LSP_optimised + " in " + duration + "
					// ns");
				}
				csv += (duration_1 / 2) + ", ";
				csv += (duration_2 / 2) + ", ";
				csv += "\n";
				pb.stepBy(stepsize);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(csv);
	}

	/**
	 * Dont try this at home -- takes at least 30 min
	 */
	public static void graph_test() {
		System.out.println(
				"Starting the large graph quality test. -- Abort if this is taking too long for you. You can expect hours of waiting.");

		File f = new File("src/Benchmark_Graphs/100k_j_d.txt");
		Graph G = Graph_Util.read_graph(f);

		G = Graph_Util.largest_connected_component(G);
		System.out.println(G.num_nodes());
		G.reduce_graph();

		double[] result = new double[G.num_nodes()];
		boolean[] sel = new boolean[G.num_edges()];
		// boolean[] visited = Graph_Util.visit_breadth_first(G_2,
		// G.get_nodes().get(4));

		int k = 10000;

		sel = GraphSolver.sequential_addition_farthest_midpoint(G, k, true);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);
		Graph_Util.write_graph(G, "SA_fm_pajek.txt", sel);
		System.out.println("Maximum distance (sequential addition farthest_midpoint): "
				+ Graph_Eval.evaluate_selection(sel, G).second);

		sel = Utility.random_selection(G.num_edges(), k);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);
		Graph_Util.write_graph(G, "random_pajek.txt", sel);
		System.out.println("Maximum distance (random select): " + Graph_Eval.evaluate_selection(sel, G).second);

		sel = GraphSolver.greedy_heuristic(G, k);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);

		Graph_Util.write_graph(G, "greedy_pajek.txt", sel);
		System.out.println("Maximum distance (greedy select): " + Graph_Eval.evaluate_selection(sel, G).second);

		sel = GraphSolver.sequential_addition_dist(G, k, true);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);
		Graph_Util.write_graph(G, "SA_dist_pajek.txt", sel);
		System.out.println(
				"Maximum distance (sequential addition dist): " + Graph_Eval.evaluate_selection(sel, G).second);

		sel = GraphSolver.sequential_addition_gap(G, k, true);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);
		Graph_Util.write_graph(G, "SA_gap_pajek.txt", sel);
		System.out
				.println("Maximum distance (sequential addition gap): " + Graph_Eval.evaluate_selection(sel, G).second);

		sel = GraphSolver.one_opt_naive(G, sel);
		result = Graph_Eval.dijkstra_gap_eval(sel, G);
		Graph_Util.write_graph(G, "1-OPT_SA_pajek.txt", sel);
		System.out.println("Maximum distance (1-OPT from sequential addition gap): "
				+ Graph_Eval.evaluate_selection(sel, G).second);
	}

	/**
	 * Comparing polyline run time
	 */
	private static void runtime_test_polyline() {
		System.out.println("Starting the polyline run time test ...");

		double range = 10.0;
		String csv = "n; uniform; greedy; parametric search; parametric binary search, graph heuristic\n";
		for (int i = 0; i < 20; i++) {
			int n = 1000;
			int k = 140;
			double[] A = Utility.random_Array(n, range);
			PolyLineSolver.uniform_selection_heuristic(A, k);
			PolyLineSolver.greedy_heuristic(A, k);
			PolyLineSolver.parametric_search(A, k);
			PolyLineSolver.parametric_search_binary_search(A, k);
			PolyLineSolver.graph_heuristic(A, k);
		}

		for (int n = 5000; n < 70000; n += 2500) {
			System.out.println(n);
			csv += n + "; ";
			int k = (int) Math.round(Math.sqrt(n));
			double[] A = Utility.random_Array(n, range);

			long duration;
			long time;

			time = System.nanoTime();
			for (int i = 0; i < 10; i++) {
				PolyLineSolver.uniform_selection_heuristic(A, k);
			}
			duration = (System.nanoTime() - time) / 10;
			csv += duration + "; ";

			time = System.nanoTime();
			for (int i = 0; i < 10; i++) {
				PolyLineSolver.greedy_heuristic(A, k);
			}
			duration = (System.nanoTime() - time) / 10;
			csv += duration + "; ";

			time = System.nanoTime();
			for (int i = 0; i < 5; i++) {
				PolyLineSolver.parametric_search(A, k);
			}
			duration = (System.nanoTime() - time) / 5;
			csv += duration + "; ";

			time = System.nanoTime();
			for (int i = 0; i < 5; i++) {
			PolyLineSolver.parametric_search_binary_search(A, k);
			}
			duration = (System.nanoTime() - time) / 5;
			csv += duration + "; ";

			time = System.nanoTime();
			PolyLineSolver.graph_heuristic(A, k);
			duration = System.nanoTime() - time;
			csv += duration + "\n";
		}
		System.out.println(csv);
	}

	/**
	 * Function comparing polyline quality
	 */
	private static void quality_test_polyline() {
		System.out.println("Starting the polyline quality test ...");
		int size = 10000;
		double range = 10.0;
		System.out.println("Creating random array of double values -- size: " + size + ", value range: " + range);
		String csv = "k, random, uniform, greedy, parametric, farthest midpoint, graph, 1-OPT\n";
		for (int k = 250; k <= 10000; k += 250) {
			csv += k + ", ";
			System.out.println(k + "/" + 10000);

			boolean[][] sel = new boolean[5][];
			double[][] polylines = new double[5][];
			// System.out.println(Arrays.toString(A));
			double gap = 0.0;

			for (int i = 0; i < 5; i++) {
				polylines[i] = Utility.random_Array(size, range);
			}

			for (int i = 0; i < 5; i++) {
				sel[i] = Utility.random_selection(size, k);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				sel[i] = PolyLineSolver.uniform_selection_heuristic(polylines[i], k);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				sel[i] = PolyLineSolver.greedy_heuristic(polylines[i], k);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				sel[i] = PolyLineSolver.parametric_search_binary_search(polylines[i], k);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				Graph G = new Graph(polylines[i]);
				sel[i] = GraphSolver.sequential_addition_farthest_midpoint(G, k, false);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				sel[i] = PolyLineSolver.graph_heuristic(polylines[i], k);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + ", ";

			gap = 0.0;
			for (int i = 0; i < 5; i++) {
				sel[i] = PolyLineSolver.one_opt(polylines[i], sel[i]);
				gap += Utility.max_gap(polylines[i], sel[i]) / 5;
			}
			csv += gap + "\n";

		}

		System.out.println(csv);
	}

}
