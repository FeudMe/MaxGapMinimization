/**
 * Class representing a tupel with two elements
 * @author Oliver Leenders
 *
 * @param <T> first type
 * @param <U> second type
 */
public class Tupel<T, U> {
	/**
	 * First element
	 */
    public T first;
    /**
     * Second element
     */
    public U second;
    /**
     * Default constructor
     * @param set_first first element
     * @param set_second second element
     */
    public Tupel (T set_first, U set_second){
        this.first = set_first;
        this.second = set_second;
    }   
    
    /**
     * Creates a string representation of a tupel
     */
    @Override
    public String toString() {
    	return "(" + first.toString() + ", " + second.toString() + ")";
    }
}