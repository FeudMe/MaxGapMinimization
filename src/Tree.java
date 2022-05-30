import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class Tree {
	// if two vertices are connected by an edge the one with lower index should be
	// parent
	public ArrayList<Node> nodes;
	public ArrayList<Edge> edges;

	public Tree() {
		nodes = new ArrayList<>();
		edges = new ArrayList<>();
	}
	
	public Tree(Graph G) {
		nodes = new ArrayList<>(G.num_nodes());
		edges = new ArrayList<>(G.num_nodes() - 1);
		for (int i = 0; i < G.num_nodes(); i++) {
			this.nodes.add(new Node(i));
		}
		boolean[] visited = new boolean[this.nodes.size()];
		// Use BFS to generate a tree from the graph
		Queue<Edge> Q = new LinkedList<Edge>();
		Q.add(G.get_edges().get(0));
		while (Q.size() > 0) {
			Edge e = Q.poll();
			if (visited[e.first] && visited[e.second]) {
				continue;
			} 
			Edge e_new = new Edge(e.first, e.second, e.length, edges.size());
			edges.add(e_new);
			nodes.get(e.first).incident_edges.put(e_new.index, e_new);
			nodes.get(e.second).incident_edges.put(e_new.index, e_new);
			if (!visited[e.first] && !visited[e.second]) {
				Node v = G.get_nodes().get(e.first);
				Node w = G.get_nodes().get(e.second);
				visited[v.key] = true;
				visited[w.key] = true;
				for (Edge e_i : v.incident_edges.values()) {
					
					Q.add(e_i);
				}
				for (Edge e_i : w.incident_edges.values()) {
					Q.add(e_i);
				}
			} else {
				Node v = (visited[e.first]) ? G.get_nodes().get(e.second) : G.get_nodes().get(e.first); 
				visited[v.key] = true;
				for (Edge e_i : v.incident_edges.values()) {
					Q.add(e_i);
				}
			}
		}
	}

	public Node get_root() {
		return (nodes.size() > 0) ? nodes.get(0) : null;
	}

	public Node get(int i) {
		return (nodes.size() > i) ? nodes.get(i) : null;
	}

	/**
	 * Adds a new root to the tree
	 * @param n
	 */
	public void add_root(Node n) {
		this.nodes.add(n);
	}

	/**
	 * Adds a new child to the graph
	 * @param c new child 
	 * @param p key of parent
	 * @param length length of new edge
	 */
	public void add_child(Node c, int p, double length) {
		this.nodes.add(c);
		Edge e = new Edge(p, c.key, length, edges.size());
		c.incident_edges.put(e.index, e);
		nodes.get(p).incident_edges.put(e.index, e);
		edges.add(e);
	}

	/**
	 * Returns whether first argument node is direct parent of second argument node
	 * 
	 * @param p candidate parent
	 * @param c candidate child
	 * @return whether first is parent of second
	 */
	public static boolean is_parent(Node p, Node c) {
		return p.key < c.key && (p.incident_edges.containsKey(c.key));
	}

	/**
	 * Returns whether first argument node is direct child of second argument node
	 * 
	 * @param c candidate child
	 * @param p candidate parent
	 * @return whether first is child of second
	 */
	public static boolean is_child(Node c, Node p) {
		return c.key > p.key && (c.incident_edges.containsKey(p.key));
	}

	@Override
	public String toString() {
		String s = "";
		for (Node n : this.nodes) {
			s += n.to_long_string() + " ";
		}
		s += "\n";
		for (Edge e : this.edges) {
			s += e.toString();
		}
		return s;
	}
}
