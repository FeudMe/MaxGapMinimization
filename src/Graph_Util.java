import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Scanner;

/**
 * Utility class for Graph stuph
 * 
 * @author Oliver Leenders
 *
 */
public abstract class Graph_Util {
	/**
	 * Private Constructor
	 */
	public Graph_Util() {
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
	 * Creates a Graph from an array representation of a polyline.
	 * 
	 * @param A array representation of a polyline
	 * @return Graph G
	 */
	public static Graph polyline_to_graph(double[] A) {
		Graph G = new Graph();
		G.add_node(new Node(0));
		for (int i = 0; i < A.length; i++) {
			G.add_node(new Node(i + 1));
			G.add_edge(new Edge(i, i + 1, A[i], i));
		}
		return G;
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
			if (n.degree() == 1) {
				distance[n.key] = 0.0;
			} else {
				distance[n.key] = Double.MAX_VALUE;
			}
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

	public static double[] dijkstra_default(boolean[] selection, Graph G) {
		double[] distance = new double[G.num_nodes()];
		for (Node n : G.get_nodes().values()) {
			distance[n.key] = Double.MAX_VALUE;
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

	public static double[] dijkstra_two_source_removal_upd(boolean[] selection, Graph G, double[] distance,
			Edge removed_edge, double nuke_radius) {
		if (!(G.get_nodes().entrySet().iterator().next().getValue() instanceof TwoDNode)) {
			return dijkstra_gap_eval(selection, G);
		} else {
			return dijkstra_gap_eval_rm_upd_subgraph(selection, G, distance, removed_edge, nuke_radius);
		}
	}

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
			if (Math.min(d_1, d_2) <= nuke_radius && n.degree() != 1) {
				distance[n.key] = Double.MAX_VALUE;
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

	private static double[] dijkstra_gap_eval_rm_upd_subgraph(boolean[] selection, Graph G, double[] distance,
			Edge removed_edge, double nuke_radius) {
		TwoDNode r_1 = (TwoDNode) G.get_nodes().get(removed_edge.first);
		TwoDNode r_2 = (TwoDNode) G.get_nodes().get(removed_edge.second);
		// create aubgraph
		TwoDGraph S = new TwoDGraph();
		double[] correct = dijkstra_gap_eval(selection, G);
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
				System.err.println("Error: distance[" + i + "] contains value " + distance[i] + " != correct " + correct[i]);
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
				System.err.println("Distance from removed edge: " + Math.min(d_1, d_2) + " vs radius " + (nuke_radius + removed_edge.length));
				System.err.println("Incident edges:");
				for(Edge e : n.incident_edges.values()) {
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

	/**
	 * Computes a distance array from a single source to all other vertices in a
	 * graph
	 * 
	 * @param S source
	 * @param G graph G
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

	private static void run_dijkstra(TwoDGraph G, double[] distance) {
		Comparator<Tupel<Double, TwoDNode>> comp = new Comparator<Tupel<Double, TwoDNode>>() {
			public int compare(Tupel<Double, TwoDNode> n_1, Tupel<Double, TwoDNode> n_2) {
				if (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] > 0) {
					return -1;
				} else if (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] < 0) {
					return 1;
				} else
					return 0;
			}
		};

		PriorityQueue<Tupel<Double, TwoDNode>> Q = new PriorityQueue<Tupel<Double, TwoDNode>>(G.num_nodes(), comp);
		for (TwoDNode n : G.get_nodes().values()) {
			Q.add(new Tupel<Double, TwoDNode>(distance[n.reduced_key], n));
		}

		while (Q.size() > 0) {
			Tupel<Double, TwoDNode> t = Q.remove();
			TwoDNode n = t.second;
			if (distance[n.reduced_key] < t.first) { // here lies the problem!!! or does it?
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.get_nodes().get(e.first).reduced_key;
				int index_second = G.get_nodes().get(e.second).reduced_key;
				double alt = distance[n.reduced_key] + e.length;

				if (e.first == n.key) {
					if (distance[index_second] > alt) {
						distance[index_second] = alt;
						Q.add(new Tupel<Double, TwoDNode>(distance[index_second], G.get_nodes().get(e.second)));
					}
				} else if (e.second == n.key) {
					if (distance[index_first] > alt) {
						distance[index_first] = alt;
						Q.add(new Tupel<Double, TwoDNode>(distance[index_first], G.get_nodes().get(e.first)));
					}
				} else {
					System.out.println("alert");
				}
			}
		}
	}

	private static void run_dijkstra(Graph G, double[] distance) {
		Comparator<Tupel<Double, Node>> comp = new Comparator<Tupel<Double, Node>>() {
			public int compare(Tupel<Double, Node> n_1, Tupel<Double, Node> n_2) {
				if (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] > 0) {
					return -1;
				} else if (distance[n_2.second.reduced_key] - distance[n_1.second.reduced_key] < 0) {
					return 1;
				} else
					return 0;
			}
		};

		PriorityQueue<Tupel<Double, Node>> Q = new PriorityQueue<Tupel<Double, Node>>(G.num_nodes(), comp);
		for (Node n : G.get_nodes().values()) {
			Q.add(new Tupel<Double, Node>(distance[n.reduced_key], n));
		}

		while (Q.size() > 0) {
			Tupel<Double, Node> t = Q.remove();
			Node n = t.second;
			if (distance[n.reduced_key] < t.first) {
				continue;
			}
			for (Edge e : n.incident_edges.values()) {
				int index_first = G.get_nodes().get(e.first).reduced_key;
				int index_second = G.get_nodes().get(e.second).reduced_key;
				double alt = distance[n.reduced_key] + e.length;

				if (e.first == n.key) {
					if (distance[index_second] > alt) {
						distance[index_second] = alt;
						Q.add(new Tupel<Double, Node>(distance[index_second], G.get_nodes().get(e.second)));
					}
				} else if (e.second == n.key) {
					if (distance[index_first] > alt) {
						distance[index_first] = alt;
						Q.add(new Tupel<Double, Node>(distance[index_first], G.get_nodes().get(e.first)));
					}
				} else {
					System.out.println("alert");
				}
			}
		}
	}

	/**
	 * Read a graph in the representation of the 100k graph
	 * 
	 * @param f file f
	 * @return graph G
	 */
	public static Graph read_graph(File f) {
		Graph G = new Graph();
		G.TWO_D_TYPE = 1;
		try {
			Scanner sc = new Scanner(f);
			int nodes = sc.nextInt();
			sc.nextLine();
			sc.nextLine();
			for (int i = 0; i < nodes; i++) {
				String[] s = sc.nextLine().split(" ");
				final double lat = Double.parseDouble(s[0]);
				final double lon = Double.parseDouble(s[1]);
				G.add_node(new TwoDNode(i, lon, lat));
			}
			int i = 0;
			while (sc.hasNextLine()) {
				String l = sc.nextLine();
				String[] values = l.split(" ");
				TwoDNode n_1 = (TwoDNode) G.get_nodes().get(Integer.parseInt(values[0]));
				TwoDNode n_2 = (TwoDNode) G.get_nodes().get(Integer.parseInt(values[1]));
				double weight = Utility.dist_lat_lon(n_1.y, n_1.x, n_2.y, n_2.x);
				Edge e = new Edge(n_1.key, n_2.key, weight, i);
				G.add_edge(e);
				i++;
			}
			sc.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return G;
	}

	public static Graph read_test_graph(File f) {
		Graph G = new Graph();
		G.TWO_D_TYPE = 2;
		try {
			Scanner sc = new Scanner(f);
			String l = sc.nextLine();
			int num_nodes = Integer.parseInt(l);
			l = sc.nextLine();
			int num_edges = Integer.parseInt(l);
			for (int i = 0; i < num_nodes; i++) {
				l = sc.nextLine();
				String[] split = l.split(" ");
				G.add_node(new TwoDNode(i, Double.parseDouble(split[0]), Double.parseDouble(split[1])));
			}
			for (int i = 0; i < num_edges; i++) {
				l = sc.nextLine();
				String[] split = l.split(" ");
				TwoDNode first = (TwoDNode) G.get_nodes().get(Integer.parseInt(split[0]));
				TwoDNode second = (TwoDNode) G.get_nodes().get(Integer.parseInt(split[1]));
				double length = Math.sqrt(
						(first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
				G.add_edge(new Edge(first.key, second.key, length, i));
			}
			sc.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return G;
	}

	/**
	 * Read a graph in the DIMACS challenge 9 format
	 * 
	 * @param f file to read from
	 * @return Graph
	 */
	public static Graph read_dimacs_graph(File f) {
		Graph G = new Graph();
		G.TWO_D_TYPE = 0;
		try {
			Scanner sc = new Scanner(f);
			String l = sc.nextLine();
			String[] words = l.split(" ");
			int n = Integer.parseInt(words[2]);
			int j = 0;
			for (int i = 0; i < n; i++) {
				sc.nextLine();
				G.add_node(new Node(i));
			}
			while (sc.hasNextLine()) {
				l = sc.nextLine();
				String[] split = l.split(" ");
				Edge e = new Edge(Integer.parseInt(split[1]) - 1, Integer.parseInt(split[2]) - 1,
						Integer.parseInt(split[3]), j);
				G.add_edge(e);
				j++;
			}
			sc.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return G;
	}

	public static Graph read_usa_graph(File f) {
		Graph G = new Graph();
		G.TWO_D_TYPE = 2;
		try {
			Scanner sc = new Scanner(f);
			String l = sc.nextLine();
			int num_nodes = Integer.parseInt(l);
			for (int i = 0; i < num_nodes; i++) {
				l = sc.nextLine();
				String[] split = l.split(" ");
				G.add_node(new TwoDNode(i, Double.parseDouble(split[1]), Double.parseDouble(split[2])));
			}
			l = sc.nextLine();
			int num_edges = Integer.parseInt(l);
			for (int i = 0; i < num_edges; i++) {
				l = sc.nextLine();
				String[] split = l.split(" ");
				l = sc.nextLine();
				String[] split_two = l.split(" ");
				G.add_edge(new Edge(Integer.parseInt(split[0]), Integer.parseInt(split[1]),
						Double.parseDouble(split_two[1]), i));
			}
			sc.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return G;
	}

	/**
	 * Get the largest connected component of a graph
	 * 
	 * @param G graph
	 * @return largest connected component as a graph
	 */
	public static Graph largest_connected_component(Graph G) {
		ArrayList<Tupel<Node, Integer>> sizes = connected_components(G);

		Tupel<Node, Integer> max = new Tupel<Node, Integer>(null, 0);
		for (Tupel<Node, Integer> t : sizes) {
			if (t.second > max.second) {
				max.first = t.first;
				max.second = t.second;
			}
		}
		boolean[] visited = new boolean[G.num_nodes()];
		iterative_BFS(G, max.first, visited);
		Graph G_2 = new Graph();
		G_2.TWO_D_TYPE = G.TWO_D_TYPE;
		for (int x = 0; x < visited.length; x++) {
			if (visited[x]) {
				G_2.add_node(G.get_nodes().get(x));
			}
		}

		for (Edge e : G.get_edges()) {
			if (visited[e.first] && visited[e.second]) {
				G_2.add_edge(e);
			}
		}

		return G_2;
	}

	/**
	 * Get the connected component containing a specific vertex
	 * 
	 * @param G graph G
	 * @param n vertex
	 * @return graph
	 */
	public static Graph connected_component(Graph G, Node n) {
		boolean[] visited = new boolean[G.num_nodes()];
		iterative_BFS(G, n, visited);
		Graph G_2 = new Graph();
		G_2.TWO_D_TYPE = G.TWO_D_TYPE;
		for (int x = 0; x < visited.length; x++) {
			if (visited[x]) {
				G_2.add_node(G.get_nodes().get(x));
			}
		}

		for (Edge e : G.get_edges()) {
			if (visited[e.first] && visited[e.second]) {
				G_2.add_edge(e);
			}
		}

		return G_2;

	}

	/**
	 * Get a list of nodes in disjoint connected components and the sizes of said
	 * connected components
	 * 
	 * @param G graph G
	 * @return list of tuples of nodes and sizes
	 */
	public static ArrayList<Tupel<Node, Integer>> connected_components(Graph G) {
		ArrayList<Tupel<Node, Integer>> sizes = new ArrayList<>();
		boolean[] checked = new boolean[G.num_nodes()];

		int i = 0;
		while (i < checked.length) {
			if (checked[i]) {
				break;
			}
			boolean[] visited = new boolean[G.num_nodes()];

			iterative_BFS(G, G.get_nodes().get(i), visited);
			for (int x = 0; x < visited.length; x++) {
				if (visited[x]) {
					checked[x] = true;
				}
			}
			int c = Utility.num_selected(visited);

			sizes.add(new Tupel<Node, Integer>(G.get_nodes().get(i), c));

			int j = i;
			while (j < visited.length) {
				if (!checked[i]) {
					break;
				}
				j++;
				i++;
			}

		}
		return sizes;
	}

	/**
	 * Implementation of breadth-first search found on Geeksforgeeks.com
	 * 
	 * @param G       graph G
	 * @param S       start node
	 * @param visited array of visited vertices
	 */
	public static void visit_breadth_first(Graph G, Node S, boolean[] visited) {
		visited[S.key] = true;
		for (Edge e : S.incident_edges.values()) {
			if (e.first == S.key && (!visited[e.second])) {
				visit_breadth_first(G, G.get_nodes().get(e.second), visited);
			} else if (!visited[e.first]) {
				visit_breadth_first(G, G.get_nodes().get(e.first), visited);
			}
		}
	}

	public static void iterative_BFS(Graph G, Node v, boolean[] visited) {
		Queue<Node> Q = new ArrayDeque<Node>();
		visited[v.key] = true;
		Q.add(v);

		while (!Q.isEmpty()) {
			v = Q.poll();
			for (Edge e : v.incident_edges.values()) {
				Node w;
				if (e.first == v.key) {
					w = G.get_nodes().get(e.second);
				} else {
					w = G.get_nodes().get(e.first);
				}
				if (!visited[w.key]) {
					// mark it as discovered and enqueue it
					visited[w.key] = true;
					Q.add(w);
				}
			}
		}
	}

	/**
	 * Write out a graph in pajek format with only selected edges
	 * 
	 * @param G        graph G
	 * @param fileName name of the desired file
	 * @param sel      selection
	 */
	public static void write_graph(Graph G, String fileName, boolean[] sel) {
		ArrayList<TwoDNode> nodes = new ArrayList<TwoDNode>();
		for (Node n : G.get_nodes().values()) {
			if (n instanceof TwoDNode) {
				nodes.add((TwoDNode) n);
			}
		}
		try {
			BufferedWriter bf = new BufferedWriter(new FileWriter("src/" + fileName));
			bf.append("*Network " + fileName + "\n*vertices " + (nodes.size()) + "\n");
			for (int i = 0; i < nodes.size(); i++) {
				StringBuilder sb = new StringBuilder();
				sb.append(nodes.get(i).key + 1);
				sb.append(" ");
				sb.append("\"" + (nodes.get(i).key + 1) + "\"");
				sb.append(" ");
				sb.append(nodes.get(i).x);
				sb.append(" ");
				sb.append(nodes.get(i).y);
				sb.append("\n");

				bf.append(sb.toString());
			}
			bf.append("*arcs\n");
			for (Edge e : G.get_edges()) {
				if (sel[e.index]) {
					StringBuilder sb = new StringBuilder();
					sb.append(e.first + 1);
					sb.append(" ");
					sb.append(e.second + 1);
					sb.append("\n");
					bf.append(sb.toString());
				}
			}

			bf.flush();
			bf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static Graph read_osm_graph(File f) {
		Graph G = new Graph();
		G.TWO_D_TYPE = 1;
		try {
			Scanner sc = new Scanner(f, "UTF-8");
			HashMap<Long, Integer> key_reductor = new HashMap<Long, Integer>();
			int line_num = 0;
			while (sc.hasNextLine()) {
				String l = sc.nextLine();
				line_num++;
				if (!l.isBlank()) {
					l = l.trim();
					// System.out.println(l);
					String[] split = l.split(" ");
					// System.out.println(split[0]);
					if (split[0].equals("<node")) {
						// System.out.println(l);
						String[] id_splt = split[1].split("\"");
						if (id_splt.length < 2) {
							System.out.println(l + " <- l.:" + line_num + " hasNextLine: " + sc.hasNextLine());
						}
						long id = Long.parseLong(id_splt[1]);
						int i = 1;
						while (i < split.length && !split[i].startsWith("lat=")) {
							i++;
						}
						String[] lat_split = split[i].split("\"");
						double lat = Double.parseDouble(lat_split[1]);
						String[] lon_split = split[i + 1].split("\"");
						double lon = Double.parseDouble(lon_split[1]);
						key_reductor.put(id, G.num_nodes());
						G.add_node(new TwoDNode(G.num_nodes(), lon, lat));
					} else if (split[0].equals("<way")) {
						int prev_node_index = -1;
						l = sc.nextLine().trim();
						line_num++;
						String[] w_splt = l.split(" ");
						while ((!l.equals("</way>")) && sc.hasNextLine()) {
							if (w_splt[0].equals("<tag")) {
								l = sc.nextLine().trim();
								line_num++;
								w_splt = l.split(" ");
								continue;
							} else if (w_splt[0].equals("<nd")) {
								String[] nd_splt = w_splt[1].split("\"");
								long nd_id_l = Long.parseLong(nd_splt[1]);
								int nd_id = key_reductor.get(nd_id_l);
								TwoDNode curr_node = (TwoDNode) G.get_nodes().get(nd_id);
								if (prev_node_index != -1) {
									TwoDNode prev_node = (TwoDNode) G.get_nodes().get(prev_node_index);
									// System.out.println(prev_node_index + ", " + curr_node.key + ", " + nd_id_l);
									G.add_edge(new Edge(prev_node_index, curr_node.key,
											Utility.dist_lat_lon(prev_node.y, prev_node.x, curr_node.y, curr_node.x),
											G.num_edges()));
								}
								prev_node_index = curr_node.key;
							}
							l = sc.nextLine().trim();
							line_num++;
							w_splt = l.split(" ");
						}
					} else {
						continue;
					}
				}
			}
			sc.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		return G;
	}

}
