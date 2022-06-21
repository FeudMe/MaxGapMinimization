import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

/**
 * Class containing functions for solving MGMP on polylines
 * 
 * @author Oliver Leenders
 *
 */
public abstract class PolyLineSolver {
	/**
	 * Private constructor
	 */
	public PolyLineSolver() {
	}

	/**
	 * Implementation of Stankov and Storandt's uniform selection heuristic
	 * 
	 * @param A polyline
	 * @param k number of segments to be selected
	 * @return selection
	 */
	public static boolean[] uniform_selection_heuristic(double[] A, int k) {
		double[] sums = new double[A.length];
		sums[0] = A[0];
		for (int i = 1; i < A.length; i++) {
			sums[i] = sums[i - 1] + A[i];
		}

		double L = sums[sums.length - 1];
		double delta_k = L / ((double) (k + 1));
		double[] delta_k_points = new double[k];
		for (int i = 0; i < delta_k_points.length; i++) {
			delta_k_points[i] = delta_k * (i + 1);
		}

		boolean[] selection = new boolean[A.length];
		int sel_index = 0;
		for (int i = 0; i < A.length; i++) {
			double curr = sums[i];

			if (sel_index > k - 1) {
				selection[i] = false;
			} else if (curr > delta_k_points[sel_index]) {
				selection[i] = true;
				sel_index++;
			} else {
				selection[i] = false;
			}
		}

		return selection;
	}

	/**
	 * Implementation of the greedy approximation selecting the k longest segments
	 * 
	 * @param A polyline
	 * @param k k
	 * @return selection
	 */
	public static boolean[] greedy_heuristic(double[] A, int k) {
		@SuppressWarnings("unchecked")
		Tupel<Double, Integer>[] sel = new Tupel[k];
		for (int i = 0; i < sel.length; i++) {
			sel[i] = new Tupel<Double, Integer>(0.0, -1);
		}
		for (int i = 0; i < A.length; i++) {
			int j = 0;
			if (sel[j].first < A[i]) {
				sel[j].first = A[i];
				sel[j].second = i;
				j++;
			}
			while (sel[j].first < A[i]) {
				sel[j - 1].first = sel[j].first;
				sel[j - 1].second = sel[j].second;
				sel[j].first = A[i];
				sel[j].second = i;
				j++;
				if (j >= sel.length) {
					break;
				}
			}
		}
		boolean[] selection = new boolean[A.length];
		for (int i = 0; i < selection.length; i++) {
			selection[i] = false;
		}
		for (int i = 0; i < sel.length; i++) {
			selection[sel[i].second] = true;
		}

		return selection;
	}

	/**
	 * Implementation of Stankov and Storand's parametric search algorithm for MGMP
	 * (exact)
	 * 
	 * @param A polyline
	 * @param k number of segments to be selected
	 * @return selection
	 */
	public static boolean[] parametric_search(double[] A, int k) {
		if (k == A.length) {

			Utility.selection_from_gap_length(A, 0);

		}
		// calculate delta_k
		double delta_k = Utility.calculate_delta_k(A, k);
		// generate all possible gaps < delta_k
		LinkedList<Double> gaps = Utility.generate_possible_gaps(A, delta_k);
		// sort array descendingly
		Collections.sort(gaps, Collections.reverseOrder());
		// filter duplicates
		Utility.removeDuplicates(gaps);
		double[] gaps_arr = new double[gaps.size()];
		Iterator<Double> it = gaps.iterator();
		{
			int i = 0;
			while (it.hasNext()) {
				gaps_arr[i] = it.next();
				i++;
			}
		}
		// for each gap length
		// compute minimal gap
		int index_of_min_gap = -1;
		int i = 0;
		for (; i < gaps.size(); i++) {
			LinkedList<Integer> selected_indices = new LinkedList<>();
			double current_gap = gaps_arr[i];
			double gap_filler = 0.0;
			for (int j = 0; j < A.length && selected_indices.size() <= k; j++) {
				if (gap_filler + A[j] > current_gap) {
					selected_indices.add(j);
					gap_filler = 0.0;
				} else {
					gap_filler += A[j];
				}
			}
			if (selected_indices.size() > k) {
				index_of_min_gap = i - 1;
				break;
			}
		}
		double min_gap = gaps_arr[index_of_min_gap];
		boolean[] selection = Utility.selection_from_gap_length(A, min_gap);
		return selection;
	}

	/**
	 * Implementation of Stankov and Storandt's parametric search algorithm for MGMP
	 * using binary search (exact)
	 * 
	 * @param A polyline
	 * @param k number of segments to be selected
	 * @return selection
	 */
	public static boolean[] parametric_search_binary_search(double[] A, int k) {
		// calculate delta_k
		double delta_k = Utility.calculate_delta_k(A, k);
		// generate all possible gaps < delta_k
		LinkedList<Double> gaps = Utility.generate_possible_gaps(A, delta_k);
		// sort array descendingly
		Collections.sort(gaps, Collections.reverseOrder());
		// filter duplicates
		Utility.removeDuplicates(gaps);
		double[] gaps_arr = new double[gaps.size()];
		Iterator<Double> it = gaps.iterator();
		{
			int i = 0;
			while (it.hasNext()) {
				gaps_arr[i] = it.next();
				i++;
			}
		}
		// for each gap length
		// compute minimal gap
		int index_to_check = gaps.size() / 2;
		int left = 0;
		int right = gaps.size() - 1;
		double prev_opt_gap = Double.MAX_VALUE;
		boolean optimal = false;
		while (!optimal) {
			LinkedList<Integer> selected_indices = new LinkedList<>();
			double current_gap = gaps_arr[index_to_check];
			double gap_filler = 0.0;
			for (int j = 0; j < A.length && selected_indices.size() <= k; j++) {
				if (gap_filler + A[j] > current_gap) {
					selected_indices.add(j);
					gap_filler = 0.0;
				} else {
					gap_filler += A[j];
				}
			}
			if (selected_indices.size() > k) {
				right = index_to_check - 1;
			} else if (selected_indices.size() < k) {
				left = index_to_check + 1;
			} else if (current_gap < prev_opt_gap) {
				prev_opt_gap = current_gap;
				left = index_to_check + 1;
			} else if (current_gap > prev_opt_gap) {
				right = index_to_check - 1;
			} else {
				optimal = true;
			}
			index_to_check = (left + right) / 2;
		}
		return Utility.selection_from_gap_length(A, prev_opt_gap);
	}

	/**
	 * Efficient implementation of the Graph gap minimzation approach for polylines
	 * 
	 * @param A polyline
	 * @param k
	 * @return selection
	 */
	public static boolean[] graph_heuristic(double[] A, int k) {
		Random rd = new Random(System.currentTimeMillis());
		boolean[] sel = new boolean[A.length];
		sel[rd.nextInt(A.length)] = true;
		for (int i = 1; i < k; i++) {
			int j = Utility.middle_segment_of_longest_gap(A, sel);
			if (sel[j]) {
				System.out.println("here" + j);
			}
			sel[j] = true;
		}
		return sel;
	}

	/**
	 * Computes a 1-OPT solution for polylines based of a given selection
	 * 
	 * @param A   polyline
	 * @param sel
	 * @return 1-OPT selection
	 */
	public static boolean[] one_opt(double[] A, boolean[] sel) {
		boolean has_changed = true;
		// System.out.println(prev_max_gap);
		double curr_max_gap_length = Utility.max_gap(A, sel);
		while (has_changed) {
			has_changed = false;
			for (int i = 0; i < A.length; i++) {
				if (sel[i]) {
					sel[i] = false;
					int j = Utility.middle_segment_of_longest_gap(A, sel);
					if (j != i) {
						sel[j] = true;
						double new_max_gap_length = Utility.max_gap(A, sel);
						if (new_max_gap_length < curr_max_gap_length) {
							has_changed = true;
							curr_max_gap_length = new_max_gap_length;
							// System.out.println(new_max_gap_length);
						} else {
							sel[j] = false;
							sel[i] = true;
						}
					} else {
						sel[i] = true;
					}
				}
			}
		}
		return sel;
	}

}
