import java.util.HashMap;

/**
 * Class representing a node with 2d coordinates
 * 
 * @author Oliver Leenders
 *
 */
public class TwoDNode extends Node {
	/**
	 * x coordinate
	 */
	public final double x;
	/**
	 * y coordinate
	 */
	public final double y;

	/**
	 * standard constructor
	 * 
	 * @param set_key index/key of the node
	 */
	public TwoDNode(int set_key) {
		super(set_key);
		this.x = 0.0;
		this.y = 0.0;
	}

	/**
	 * constructor, with given incident edges
	 * 
	 * @param set_key            index/key of the node
	 * @param set_incident_edges oncident edges of the node
	 */
	public TwoDNode(int set_key, HashMap<Integer, Edge> set_incident_edges) {
		super(set_key, set_incident_edges);
		this.x = 0.0;
		this.y = 0.0;
	}

	/**
	 * Constructor for 2D - nodes
	 * 
	 * @param set_key index/key of the node
	 * @param set_x   set x coordinate
	 * @param set_y   set y coordinate
	 */
	public TwoDNode(int set_key, final double set_x, final double set_y) {
		super(set_key);
		this.x = set_x;
		this.y = set_y;
	}

	/**
	 * Create a string representation of a 2d node
	 * 
	 * @return string representation of a 2d node
	 */
	@Override
	public String toString() {
		return "(" + this.key + " (" + this.reduced_key + ")) -- " + this.x + ", " + this.y;
	}

}
