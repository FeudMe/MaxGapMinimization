import java.io.File;
import java.text.NumberFormat;
import java.util.Locale;

/**
 * Driver Code for Maximum Gap Minimization testing
 * 
 * @author Oliver Leenders
 *
 */
public class Main {
	public static void main(String[] args) {
		System.out.println("This will take a while to finish.\nSeveral tests will run and print their results as csv:");
		System.out.println("Starting the polyline run time test ...");
		// runtime_test_polyline();
		System.out.println("Starting the polyline quality test ...");
		// quality_test_polyline();
		System.out.println("Starting the graph quality test ...");
		quality_test_graph();
		System.out.println(
				"Starting the large graph quality test. -- Abort if this is taking too long for you. You can expect at least several hours of waiting.");
		//graph_test();
	}

	/**
	 * Function testing solution quality on graphs
	 */
	private static void quality_test_graph() {
		File f = new File("src/map.osm");
		Graph G = Graph_Util.read_osm_graph(f);
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());

		System.out.println("Reducing graph ...");
		
		G = Graph_Util.largest_connected_component(G);
		G.reduce_graph();
		
		// System.out.println(G.toString());
		System.out.println("|V| = " + G.num_nodes() + ", |E| = " + G.num_edges());
		System.out.println(Graph_Util.connected_components(G));
		
		String csv = "k; random gap; greedy gap; sequential addition gap; MGM gap; 1-OPT fast gap; 1-OPT thorough gap\n";
		String csv_dist = "k; random; greedy; sequential addition; MGM dist; 1-OPT fast; 1-OPT thorough\n";
		int stepsize = G.num_edges() / 15;
		boolean[] dist_sel = new boolean[G.num_edges()];
		boolean[] gap_sel = new boolean[G.num_edges()];
		boolean[] rd_sel = new boolean[G.num_edges()];
		boolean[] greedy_sel = new boolean[G.num_edges()];
		boolean[] one_opt_sel = new boolean[G.num_edges()];
		for (int i = stepsize; i < G.num_edges(); i += stepsize) {
			csv += i + "; ";
			System.out.println(i + "/" + G.num_edges() + " edges to be selected");
			double quality = 0.0;
			double quality_gap = 0.0;

			int rounding_runs = 25;
			for (int j = 0; j < rounding_runs; j++) {
				rd_sel = Utility.random_selection(G.num_edges(), i);
				quality += Utility.max(Graph_Util.dijkstra_default(rd_sel, G)) / rounding_runs;
				quality_gap += Graph_Util.evaluate_selection(rd_sel, G).second / rounding_runs;
			}
			csv_dist += String.format("%1$,.7f", quality) + "; ";
			csv += String.format("%1$,.7f", quality_gap) + "; ";

			greedy_sel = GraphSolver.greedy_heuristic(G, i);
			csv_dist += String.format("%1$,.7f", Utility.max(Graph_Util.dijkstra_default(greedy_sel, G))) + "; ";
			csv += String.format("%1$,.7f", Graph_Util.evaluate_selection(greedy_sel, G).second) + "; ";

			dist_sel = GraphSolver.k_seq_add_dist(dist_sel, G, stepsize);
			csv_dist += String.format("%1$,.7f", Utility.max(Graph_Util.dijkstra_default(dist_sel, G))) + "; ";
			csv += String.format("%1$,.7f", Graph_Util.evaluate_selection(dist_sel, G).second) + "; ";

			gap_sel = GraphSolver.k_seq_add_gap(gap_sel, G, stepsize);
			csv_dist += String.format("%1$,.7f", Utility.max(Graph_Util.dijkstra_default(gap_sel, G))) + "; ";
			csv += String.format("%1$,.7f", Graph_Util.evaluate_selection(gap_sel, G).second) + "; ";

			
			long time = System.currentTimeMillis();
			one_opt_sel = GraphSolver.one_opt(G, greedy_sel.clone());
			long duration = System.currentTimeMillis() - time;
			
			
			csv_dist += String.format("%1$,.7f", Utility.max(Graph_Util.dijkstra_default(one_opt_sel, G))) + "; ";
			csv += String.format("%1$,.7f", Graph_Util.evaluate_selection(one_opt_sel, G).second) + "; ";
			
			time = System.currentTimeMillis();
			one_opt_sel = GraphSolver.one_opt_naive(G, greedy_sel.clone());
			long duration_naive = System.currentTimeMillis() - time;
			System.out.println("smart 1-OPT / naive 1-OPT " + duration / (float) duration_naive);
			csv_dist += String.format("%1$,.7f", Utility.max(Graph_Util.dijkstra_default(one_opt_sel, G))) + "\n";
			csv += String.format("%1$,.7f", Graph_Util.evaluate_selection(one_opt_sel, G).second) + "\n";
		}
		System.out.println(csv);
		//System.out.println(csv_dist);
	}

	/**
	 * Dont try this at home -- takes at least 30 min
	 */
	public static void graph_test() {

		File f = new File("src/100k_j_d.txt");
		Graph G_2 = Graph_Util.read_graph(f);
		G_2 = Graph_Util.largest_connected_component(G_2);
		System.out.println(G_2.num_nodes());
		G_2.reduce_graph();
		double[] result = new double[G_2.num_nodes()];
		boolean[] sel = new boolean[G_2.num_edges()];
		// boolean[] visited = Graph_Util.visit_breadth_first(G_2,
		// G.get_nodes().get(4));

		int k = 10000;

		sel = Utility.random_selection(G_2.num_edges(), k);
		result = Graph_Util.dijkstra_gap_eval(sel, G_2);
		Graph_Util.write_graph(G_2, "random_pajek.txt", sel);
		System.out.println("Maximum distance (random select): " + Graph_Util.evaluate_selection(sel, G_2).second);

		sel = GraphSolver.greedy_heuristic(G_2, k);
		result = Graph_Util.dijkstra_gap_eval(sel, G_2);

		Graph_Util.write_graph(G_2, "greedy_pajek.txt", sel);
		System.out.println("Maximum distance (greedy select): " + Graph_Util.evaluate_selection(sel, G_2).second);

		sel = GraphSolver.sequential_addition_dist(G_2, k, true);
		result = Graph_Util.dijkstra_gap_eval(sel, G_2);
		Graph_Util.write_graph(G_2, "SA_dist_pajek.txt", sel);
		System.out.println(
				"Maximum distance (sequential addition dist): " + Graph_Util.evaluate_selection(sel, G_2).second);

		sel = GraphSolver.sequential_addition_gap(G_2, k, true);
		result = Graph_Util.dijkstra_gap_eval(sel, G_2);
		Graph_Util.write_graph(G_2, "SA_gap_pajek.txt", sel);
		System.out.println(
				"Maximum distance (sequential addition gap): " + Graph_Util.evaluate_selection(sel, G_2).second);

		sel = GraphSolver.one_opt(G_2, sel);
		result = Graph_Util.dijkstra_gap_eval(sel, G_2);
		Graph_Util.write_graph(G_2, "1-OPT_SA_pajek.txt", sel);
		System.out.println("Maximum distance (1-OPT from sequential addition gap): "
				+ Graph_Util.evaluate_selection(sel, G_2).second);
	}

	/**
	 * Comparing polyline run time
	 */
	private static void runtime_test_polyline() {
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
			PolyLineSolver.parametric_search(A, k);
			duration = System.nanoTime() - time;
			csv += duration + "; ";

			time = System.nanoTime();
			PolyLineSolver.parametric_search_binary_search(A, k);
			duration = System.nanoTime() - time;
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
		int size = 10000;
		double range = 10.0;
		System.out.println("Creating random array of double values -- size: " + size + ", value range: " + range);
		double[] A = Utility.random_Array(size, range);
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.GERMAN);
		String csv = "k; random; uniform; greedy; parametric; graph; 1-OPT\n";
		for (int k = 500; k <= 10000; k += 500) {
			csv += k + "; ";
			System.out.println(k + "/" + 10000);
			// System.out.println(Arrays.toString(A));

			boolean[] sel = Utility.random_selection(size, k);
			double gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "; ";

			sel = PolyLineSolver.uniform_selection_heuristic(A, k);
			gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "; ";

			sel = PolyLineSolver.greedy_heuristic(A, k);
			gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "; ";

			sel = PolyLineSolver.parametric_search_binary_search(A, k);
			gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "; ";

			sel = PolyLineSolver.graph_heuristic(A, k);
			gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "; ";

			sel = PolyLineSolver.one_opt(A, sel);
			gap = Utility.max_gap(A, sel);
			csv += nf.format(gap) + "\n";

		}

		System.out.println(csv);
	}

}
