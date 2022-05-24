import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
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
		int[] checked = new int[G.num_nodes()];
		for (int i = 0; i < checked.length; i++) {
			checked[i] = -1;
		}
		for (int i = 0; i < checked.length; i++) {
			if (checked[i] != -1) {
				continue;
			}
			int c = iterative_BFS(G, G.get_nodes().get(i), checked);
				
			sizes.add(new Tupel<Node, Integer>(G.get_nodes().get(i), c));
		}
		return sizes;
	}

	public static Graph subgraph_dimacs() {
		return null;
	}
	
	public static Graph subgraph_100_k(double percentage) {
		File f = new File("src/Benchmark_Graphs/100k_j_d.txt");
		Graph G = Graph_Util.read_graph(f);

		double x_min = Double.MAX_VALUE;
		double x_max = -Double.MAX_VALUE;
		double y_min = Double.MAX_VALUE;
		double y_max = -Double.MAX_VALUE;
		
		
		for (Node n : G.get_nodes().values()) {
			TwoDNode n_2d = (TwoDNode) n;
			x_max = Math.max(x_max, n_2d.x);
			x_min = Math.min(x_min, n_2d.x);
			y_max = Math.max(y_max, n_2d.y);
			y_min = Math.min(y_min, n_2d.y);
		}
		// System.out.println(x_min + " " + x_max + " " + y_min + " " + y_max);
		double side_length_factor = Math.sqrt(percentage);
		double new_x_length = side_length_factor * (x_max - x_min);
		double new_y_length = side_length_factor * (y_max - y_min);
		
		double new_max_x_cutoff = x_min + new_x_length;
		double new_max_y_cutoff = y_min + new_y_length;
		// System.out.println(new_max_x_cutoff + " " + new_max_y_cutoff);
		
		Graph G_2 = new Graph();
		Iterator<Node> it = G.get_nodes().values().iterator();
		
		while (it.hasNext()) {
			TwoDNode n_2d = (TwoDNode) it.next();
			if (n_2d.x <= new_max_x_cutoff || n_2d.y <= new_max_y_cutoff) {
				G_2.add_node(new TwoDNode(n_2d.key, n_2d.x, n_2d.y));
			}
		}
		for (Edge e : G.get_edges()) {
			if (G_2.get_nodes().containsKey(e.first) && G_2.get_nodes().containsKey(e.second)) {
				G_2.add_edge(e);
			}
		}
		G_2.reduce_graph();
		G_2 = largest_connected_component(G_2);
		G_2.reduce_graph();
		return G_2;
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

	public static int iterative_BFS(Graph G, Node v, boolean[] visited) {
		Queue<Node> Q = new ArrayDeque<Node>();
		visited[v.key] = true;
		Q.add(v);
		
		int c = 0; 
		
		while (!Q.isEmpty()) {
			v = Q.poll();
			c++;
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
		return c;
	}

	public static int iterative_BFS(Graph G, Node v, int[] checked) {
		Queue<Node> Q = new ArrayDeque<Node>();
		checked[v.key] = v.key;
		Q.add(v);
		
		int c = 0; 
		
		while (!Q.isEmpty()) {
			Node next = Q.poll();
			c++;
			for (Edge e : next.incident_edges.values()) {
				Node w;
				if (e.first == next.key) {
					w = G.get_nodes().get(e.second);
				} else {
					w = G.get_nodes().get(e.first);
				}
				if (checked[w.key] == -1) {
					// mark it as discovered and enqueue it
					checked[w.key] = v.key;
					Q.add(w);
				}
			}
		}
		return c;
	}
	
	public static int iterative_BFS(Graph G, Node v, boolean[] visited, boolean[] checked) {
		Queue<Node> Q = new ArrayDeque<Node>();
		visited[v.key] = true;
		checked[v.key] = true;
		Q.add(v);
		
		int c = 0; 
		
		while (!Q.isEmpty()) {
			v = Q.poll();
			c++;
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
					checked[w.key] = true;
					Q.add(w);
				}
			}
		}
		return c;
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
			BufferedWriter bf = new BufferedWriter(new FileWriter("src/Results/" + fileName));
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
	
	public static Tree read_tree(File f) {
		Tree T = new Tree();
		try {
			Scanner sc = new Scanner(f);
			while (sc.hasNextLine()) {
				String[] split = sc.nextLine().split(" ");
				int id_1 = Integer.parseInt(split[0]);
				int id_2 = Integer.parseInt(split[1]);
				double l = Double.parseDouble(split[2]);
				Node n_1 = new Node(id_1);
				Node n_2 = new Node(id_2);
				if (T.nodes.size() == 0) {
					T.add_root(n_1);
					T.add_child(n_2, n_1.key, l);
				} else {
					if (id_1 < id_2) {
						T.add_child(n_2, n_1.key, l);
					} else {
						T.add_child(n_1, n_2.key, l);
					}
				}
			}
			sc.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return T;
	}

}
