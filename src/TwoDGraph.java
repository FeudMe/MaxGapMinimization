import java.util.ArrayList;
import java.util.HashMap;

/**
 * Simple graph datastructure allowing indexation of edges and vertices
 * 
 * @author Oliver Leenders
 *
 */
public class TwoDGraph {

	/**
	 * The edges of the graph
	 */
	private ArrayList<Edge> edges;
	/**
	 * The vertices of the graph
	 */
	private HashMap<Integer, TwoDNode> nodes;

	/**
	 * Constructor for an empty graph
	 */
	public TwoDGraph() {
		this.edges = new ArrayList<Edge>();
		this.nodes = new HashMap<Integer, TwoDNode>();
	}

	/**
	 * Add a node to the graph
	 * 
	 * @param n node to be added
	 */
	public void add_node(TwoDNode n) {
		this.nodes.put(n.key, n);
	}

	/**
	 * Add an edge to the graph. If the endpoints are not contained in the graph,
	 * add them as well.
	 * 
	 * @param e edge to be added
	 */
	public void add_edge(Edge e) {
		if ((!this.nodes.containsKey(e.first)) || (!this.nodes.containsKey(e.second))) {
			System.err.println("Nodes not found in 2D-Graph!");
			System.err.println("Edge " + e.toString());
			System.err.println("Node " + e.first + " found: " + this.nodes.containsKey(e.first));
			System.err.println("Node " + e.second + " found: " + this.nodes.containsKey(e.second));
			return;
		}
		this.nodes.get(e.first).add_incident_edge(e);
		this.nodes.get(e.second).add_incident_edge(e);
		this.edges.add(e);
	}

	/**
	 * Remove an edge from the graph
	 * 
	 * @param e edge to be removed
	 */
	public void remove_edge(Edge e) {
		this.edges.remove(e.index);
	}

	/**
	 * Remove a node from the graph
	 * 
	 * @param n node to be removed
	 */
	public void removeNode(TwoDNode n) {
		for (Edge e : n.incident_edges.values()) {
			this.edges.remove(e.index);
		}
		this.nodes.remove(n.key);
	}

	/**
	 * Get the number of nodes in the graph
	 * 
	 * @return number of nodes
	 */
	public int num_nodes() {
		return this.nodes.size();
	}

	/**
	 * Get the number of edges in the graph
	 * 
	 * @return number of edges
	 */
	public int num_edges() {
		return this.edges.size();
	}

	/**
	 * Get the nodes of the graph
	 * 
	 * @return hashmap of integers to nodes
	 */
	public HashMap<Integer, TwoDNode> get_nodes() {
		return this.nodes;
	}

	/**
	 * Get the edges of the graph
	 * 
	 * @return ArrayList of edges
	 */
	public ArrayList<Edge> get_edges() {
		return this.edges;
	}

	/**
	 * Update the indices of the vertices and edges such that all indices start from
	 * zero and are continously increasing
	 */
	public void reduce_graph() {
		int i = 0;
		for (Node n : this.nodes.values()) {
			n.reduce_key(i);
			i++;
		}
		i = 0;
		for (Edge e : this.edges) {
			e.reduce_edge(i);
			e.first = this.nodes.get(e.first).reduced_key;
			e.second = this.nodes.get(e.second).reduced_key;
			i++;
		}
		HashMap<Integer, TwoDNode> replace = new HashMap<>();
		for (TwoDNode n : this.nodes.values()) {
			n.key = n.reduced_key;
			replace.put(n.key, n);
		}
		this.nodes = replace;

	}
	/**
	 * Return a string representation of the graph
	 */
	@Override
	public String toString() {
		String s = "";
		for (Node n : this.nodes.values()) {
			s += n.to_long_string();
		}
		return s;
	}
}
