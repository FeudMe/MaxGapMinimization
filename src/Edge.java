/**
 * Class representing undirected weighted edges supporting indexation
 * 
 * @author Oliver Leenders
 *
 */
public class Edge {
	/**
	 * First key
	 */
	public int first;
	/**
	 * Second key
	 */
	public int second;
	/**
	 * length of the edge
	 */
	public final double length;
	/**
	 * index of the edge
	 */
	public int index;

	/**
	 * Default constructor
	 * 
	 * @param set_first  set first key
	 * @param set_second set second key
	 * @param set_length set edge length
	 * @param set_index  set index of edge
	 */
	public Edge(final int set_first, final int set_second, final double set_length, final int set_index) {
		this.first = set_first;
		this.second = set_second;
		this.length = set_length;
		this.index = set_index;
	}

	/**
	 * Overriding equals such that undirected edge behaviour is mimicked
	 */
	@Override
	public boolean equals(Object o) {
		if (o instanceof Edge) {
			return ((Edge) o).first == this.first && ((Edge) o).second == this.second
					&& ((Edge) o).length == this.length && ((Edge) o).index == this.index;

		}
		return false;
	}

	/**
	 * Update index of edge
	 * 
	 * @param set_index new index
	 */
	public void reduce_edge(int set_index) {
		this.index = set_index;
	}

	/**
	 * Override hashcode such that collision free hashing is possible
	 */
	@Override
	public int hashCode() {
		return this.index;
	}

	/**
	 * Generate a string representation of the edge
	 */
	@Override
	public String toString() {
		return "{" + first + ", " + second + "} weight: " + length;
	}

}
