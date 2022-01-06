import java.util.HashMap;

/**
 * Class representing a node
 * 
 * @author Oliver Leenders
 *
 */
public class Node {
	/**
	 * Key/index of the node
	 */
	public int key;
	/**
	 * reduced key/index of the node
	 */
	public int reduced_key;
	/**
	 * Incident edges of the node
	 */
	public HashMap<Integer, Edge> incident_edges;

	/**
	 * Default constructor
	 * 
	 * @param set_key index/key of the node
	 */
	public Node(final int set_key) {
		this.key = set_key;
		this.reduced_key = set_key;
		this.incident_edges = new HashMap<Integer, Edge>();
	}

	/**
	 * Constructor of a node with given incident edges
	 * 
	 * @param set_key            index/key of the node
	 * @param set_incident_edges incident edges of the node
	 */
	public Node(final int set_key, HashMap<Integer, Edge> set_incident_edges) {
		this.key = set_key;
		this.reduced_key = set_key;
		this.incident_edges = set_incident_edges;
	}

	/**
	 * Get the degree of the node
	 * 
	 * @return degree
	 */
	public int degree() {
		return this.incident_edges.size();
	}

	/**
	 * Set the reduced key of the node
	 * 
	 * @param set_reduced_key the reduced key
	 */
	public void reduce_key(int set_reduced_key) {
		this.reduced_key = set_reduced_key;
	}

	/**
	 * Add an incident edge to the node
	 * 
	 * @param e edge to be added
	 */
	public void add_incident_edge(Edge e) {
		this.incident_edges.put(e.index, e);
	}

	/**
	 * Remove an incident edge of the node
	 * 
	 * @param e edge to be removed
	 */
	public void remove_incident_edge(Edge e) {
		this.incident_edges.remove(e.index);
	}

	/**
	 * Clear the incident edges of the node
	 */
	public void clear_incident_edges() {
		this.incident_edges.clear();
	}

	/**
	 * Change the index of the node
	 * 
	 * @param i new index
	 */
	public void change_index(int i) {
		this.key = i;
	}

	/**
	 * Overriding hashcode such that nodes can be hashed without collisions
	 */
	@Override
	public int hashCode() {
		return this.key;
	}

	/**
	 * Override equals such that keys are compared
	 */
	@Override
	public boolean equals(Object o) {
		if (o instanceof Node) {
			return ((Node) o).key == this.key;
		}
		return false;
	}

	/**
	 * Generates a string representation of a node
	 * @return string representation
	 */
	@Override
	public String toString() {
		String s = "(" + this.key + ")";
		return s;
	}

	/**
	 * Generates a longer string representation 
	 * @return string representation
	 */
	public String to_long_string() {
		String s = "(" + this.key + ")\n";
		for (Edge e : this.incident_edges.values()) {
			s += "\t" + e.toString() + "\n";
		}
		return s;
	}
}
