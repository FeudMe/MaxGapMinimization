import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

/**
 * Different Utility functions
 * 
 * @author Oliver Leenders
 */
public abstract class Utility {
	/**
	 * Private Constructor
	 */
	private Utility() {
	}

	/**
	 * Removes duplicates from a sorted {@code ArrayList} of {@code double}-s
	 * 
	 * @param L list to remove duplicates from
	 */
	public static void removeDuplicates(ArrayList<Double> L) {
		for (int i = 0; i < L.size() - 1; i++) {
			double here = L.get(i);
			double next = L.get(i + 1);
			if (here == next) {
				L.remove(i);
				i--;
			}
		}
	}
	
	public static void removeDuplicates(LinkedList<Double> L) {
		double prev = -1.0d;
		Iterator<Double> it = L.iterator();
		while (it.hasNext()) {
			double curr = it.next();
			if (curr == prev) {
				it.remove();
			} else {
				prev = curr;
			}
		}
	}

	/**
	 * Creates a array of random numbers within a range
	 * 
	 * @param size  size of the array
	 * @param range upper range limit
	 * @return random array
	 */
	public static double[] random_Array(int size, double range) {
		double[] A = new double[size];
		Random rd = new Random();
		for (int i = 0; i < size; i++) {
			A[i] = rd.nextDouble() * range;
		}
		return A;
	}

	/**
	 * Generate a list of possible gap lengths.
	 * 
	 * @param A       polyline
	 * @param delta_k delta_k
	 * @return list of possible gap lengths smaller than delta_k
	 */
	public static ArrayList<Double> generate_possible_gaps(double[] A, double delta_k) {
		ArrayList<Double> gaps = new ArrayList<Double>();
		// create array of all possible gap lengths
		for (int i = 0; i < A.length; i++) {
			int j = i;
			double gap = 0.0;
			while (j < A.length && gap < delta_k) {
				gap += A[j];
				if (gap < delta_k) {
					gaps.add(gap);
				}
				j++;
			}
		}
		gaps.add(0.0);
		return gaps;
	}

	/**
	 * generate a selection realizing a certaing gap length
	 * 
	 * @param A          polyline
	 * @param gap_length gap length to be realized
	 * @return selection
	 */
	public static boolean[] selection_from_gap_length(double[] A, double gap_length) {
		boolean[] selection = new boolean[A.length];
		double gap_filler = 0.0;
		for (int i = 0; i < A.length; i++) {
			double a = A[i];
			if (gap_filler + a > gap_length) {
				selection[i] = true;
				gap_filler = 0.0;
			} else {
				selection[i] = false;
				gap_filler += A[i];
			}
		}
		return selection;
	}

	/**
	 * Calculate delta_k for a polyline
	 * 
	 * @param A polyline
	 * @param k number of segments to be selected
	 * @return delta_k
	 */
	public static double calculate_delta_k(double[] A, int k) {
		double delta_k = 0.0;
		for (double i : A) {
			delta_k += i;
		}
		delta_k /= k + 1;
		return delta_k;
	}

	/**
	 * Evaluate a selection and return the maximum gap
	 * 
	 * @param A         polyline
	 * @param selection selection
	 * @return length of the maximum gap
	 */
	public static double max_gap(double[] A, boolean[] selection) {
		double max_gap = 0.0;
		double gap = 0.0;
		for (int i = 0; i < A.length; i++) {
			if (selection[i]) {
				max_gap = Math.max(gap, max_gap);
				gap = 0.0;
			} else {
				gap = gap + A[i];
			}
		}
		return Math.max(max_gap, gap);
	}

	public static double[] min(double[] first, double[] second) {
		double[] res = new double[Math.min(first.length, second.length)];
		for (int i = 0; i < res.length; i++) {
			res[i] = Math.min(first[i], second[i]);
		}
		return res;
	}

	public static double dist_lat_lon(double lat_1, double lon_1, double lat_2, double lon_2) {
		double p = 0.017453292519943295; // Math.PI / 180
		double a = 0.5 - Math.cos((lat_2 - lat_1) * p) / 2
				+ Math.cos(lat_1 * p) * Math.cos(lat_2 * p) * (1 - Math.cos((lon_2 - lon_1) * p)) / 2;

		return 12742 * Math.asin(Math.sqrt(a)); // 2 * R; R = 6371 km
	}

	/**
	 * Check how many segments are selected in a selection
	 * 
	 * @param sel selection
	 * @return number of selected segments
	 */
	public static int num_selected(boolean[] sel) {
		int i = 0;
		for (boolean b : sel) {
			if (b)
				i++;
		}
		return i;
	}

	/**
	 * Return the maximum value of a double array
	 * 
	 * @param D double array
	 * @return maximum value
	 */
	public static double max(double[] D) {
		double max = 0.0;
		for (double d : D) {
			if (d > max) {
				max = d;
			}
		}
		return max;
	}

	/**
	 * Get the index of the maximum element in a double array
	 * 
	 * @param D double array
	 * @return index of maximum element
	 */
	public static int index_of_max(double[] D) {
		int index = 0;
		for (int i = 1; i < D.length; i++) {
			if (D[i] > D[index]) {
				index = i;
			}
		}
		return index;
	}

	/**
	 * Generate a random selection of k elements
	 * 
	 * @param n number of elements
	 * @param k number of elements to be selected
	 * @return selection
	 */
	public static boolean[] random_selection(int n, int k) {
		ArrayList<Boolean> sel = new ArrayList<Boolean>(n);
		for (int i = 0; i < n; i++) {
			if (i < k) {
				sel.add(true);
			} else {
				sel.add(false);
			}
		}
		Random rnd = new Random(System.currentTimeMillis());
		Collections.shuffle(sel, rnd);

		boolean[] sel_arr = new boolean[n];
		for (int i = 0; i < n; i++) {
			sel_arr[i] = sel.get(i);
		}

		return sel_arr;
	}

	/**
	 * Get the length of the maximum gap along with the index of the last edge in
	 * the gap
	 * 
	 * @param A         polyline
	 * @param selection selection
	 * @return Tupel of max gap length and index of last gap edge
	 */
	public static Tupel<Double, Integer> max_gap_with_index(double[] A, boolean[] selection) {
		double max_gap = 0.0;
		double gap = 0.0;
		int index = 0;
		for (int i = 0; i < A.length; i++) {
			if (selection[i]) {
				if (gap > max_gap) {
					max_gap = gap;
					index = i - 1;
				}
				gap = 0.0;
			} else {
				gap += A[i];
			}
		}
		if (gap > max_gap) {
			max_gap = gap;
			index = A.length - 1;
		}
		return new Tupel<Double, Integer>(max_gap, index);
	}

	/**
	 * Find the middle segment of the longest gap
	 * 
	 * @param A   polyline
	 * @param sel selection
	 * @return index of middle segment of the longest gap
	 */
	public static int middle_segment_of_longest_gap(double[] A, boolean[] sel) {
		Tupel<Double, Integer> t = max_gap_with_index(A, sel);
		double max_gap = t.first;
		int index = t.second;
		double gap = 0.0;

		while (gap + A[index] < max_gap / 2.0) {
			gap += A[index];
			index--;
		}

		return index;
	}
}