import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
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

}
