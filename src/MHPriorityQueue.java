public class MHPriorityQueue {
	Object MinHeap[];
	int map[];
	int size = 0;
	int capacity;

	public MHPriorityQueue(int set_capacity) {
		capacity = set_capacity;
		map = new int[capacity];
		for (int i = 0; i < capacity; i++) {
			map[i] = -1;
		}
		MinHeap = new Object[capacity];
	}

	@SuppressWarnings("unchecked")
	public void insert(Node n, double key) {
		if (size == capacity) {
			return;
		}
		int i = size;
		// sift up
		while (i > 0) {
			int parent = (i - 1) >>> 1;
			
			Tupel<Double, Node> p = (Tupel<Double, Node>) MinHeap[parent];
			// if el is at correct position
			if (key >= p.first) {
				break;
			}
			// else swap target with origin
			MinHeap[i] = p;
			map[p.second.reduced_key] = i; 
			i = parent;
		}
		MinHeap[i] = new Tupel<Double, Node>(key, n);
		map[n.reduced_key] = i;
		size++;
		return;
	}

	@SuppressWarnings("unchecked")
	public Tupel<Double, Node> extract_min() {
		Tupel<Double, Node> t = (Tupel<Double, Node>) MinHeap[0];
		MinHeap[0] = MinHeap[size - 1];
		map[t.second.reduced_key] = -1;
		Tupel<Double, Node> sift = (Tupel<Double, Node>) MinHeap[size - 1];
		if (size > 0) {
			int i = 0;
			int half = size >>> 1;
			while (i < half) {
				int child = (i << 1) + 1;
				Tupel<Double, Node> c = (Tupel<Double, Node>) MinHeap[child];
				int right = child + 1;
				if (right < size && ((Tupel<Double, Node>)MinHeap[right]).first < c.first) {
					c = (Tupel<Double, Node>) MinHeap[child = right];
				}
				if (sift.first <= c.first) {
					break;
				}
				MinHeap[i] = c;
				map[c.second.reduced_key] = i;
				i = child;
			}
			MinHeap[i] = sift;
			map[sift.second.reduced_key] = i;
		}
		size--;
		return t;
	}

	@SuppressWarnings("unchecked")
	public void decrease_key(Node n, double key) {
		int i = map[n.reduced_key];
		Tupel<Double, Node> t = (Tupel<Double, Node>) MinHeap[i];
		t.first = key;
		while (i > 0) {
			int parent = (i - 1) >>> 1;
			
			Tupel<Double, Node> p = (Tupel<Double, Node>) MinHeap[parent];
			// if el is at correct position
			if (key >= p.first) {
				break;
			}
			// else swap target with origin
			MinHeap[i] = p;
			map[p.second.reduced_key] = i; 
			i = parent;
		}
		MinHeap[i] = new Tupel<Double, Node>(key, n);
		map[n.reduced_key] = i;
		return;
	}

	@SuppressWarnings("unchecked")
	public void print() {
		for (int i = 0; i < size; i++) {
			System.out.print(((Tupel<Double, Node>)MinHeap[i]).first + ", ");
		}
		System.out.println();
		for (int i = 0; i < size; i++) {
			System.out.print(map[i] + " ");
		}
		System.out.println();
	}

	public boolean contains(Node n) {
		return map[n.reduced_key] != -1;
	}
	
	public boolean contains(int key) {
		if (key < capacity) {
			return map[key] != -1;
		}
		return false;
	}

	public int size() {
		return size;
	}
}
