
public class Tour {
	private class Node {
		private Point p;
		private Node next;

		public Node(Point p) {
			this.p = p;
			next = null;
		}

		public void setNext(Node next) {
			this.next = next;
		}

		public Node getNext() {
			return next;
		}

		public Point getPoint() {
			return p;
		}
	}

	private Node start;
	private int size;
	private double totalDistance;
	private double weightedDistance;
	private double load;
	private static final double capacity = 1000;
	private static final double sleighWeight = 10;
	
	public double getCapacity(){
		return capacity;
	}
	
	public Tour(){
		Point p = new Point(-1, 90, 0, sleighWeight);
		start = new Node(p);
		start.setNext(start);
		size = 1;
		totalDistance = 0;
		weightedDistance = 0;
		load = sleighWeight;
	}
	
	public int size() {
		return size;
	}
	
	public void show(){

		Node curr = start;
		for(int i = 0; i < size; i++){
			StdOut.printf("%s", curr.getPoint().toString());
			curr = curr.getNext();
		}
		//should be circular
		if(curr.getNext()!=null)
			StdOut.printf("%s", curr.getNext().getPoint().toString());
	}
	
	public double distance(){
		
		Node curr = start;
		double dist = 0;

		for(int i = 0; i < size; i++){
			dist = dist + curr.getPoint().distanceTo(curr.getNext().getPoint());
			curr = curr.getNext();
		}
		return dist;
	}
	
    public void insertNearest(Point p){
    	Node node = new Node(p);
    	
    	Node nearestNode = null;
    	Node curr = start;
    	double nearest = -1;

    	for(int i = 0; i < size; i++){
    		double dist = curr.getPoint().distanceTo(p);
    		if(dist < nearest || nearest < 0){
    			nearest = dist;
    			nearestNode = curr;
    		}
    		curr = curr.getNext();
    	}
    	
    	Node oldNext = nearestNode.getNext();
    	nearestNode.setNext(node);
    	node.setNext(oldNext);
    	size++;
    }
    
    public void insertSmallest(Point p){
    	Node node = new Node(p);
    	
    	Node smallestNode = null;
    	Node curr = start;
    	double smallest = -1;
    	double currDistanceTo = -1;
    	double nextDistanceTo = -1;
    	for(int i = 0; i < size; i++){
    		if(nextDistanceTo!=-1){
    			currDistanceTo = nextDistanceTo;
    			nextDistanceTo = curr.getNext().getPoint().distanceTo(p);
    		}else{
    			currDistanceTo = curr.getPoint().distanceTo(p);
    			nextDistanceTo = curr.getNext().getPoint().distanceTo(p);
    		}
    		double dist = totalDistance - curr.getPoint().distanceTo(curr.getNext().getPoint());
    		dist += currDistanceTo + nextDistanceTo;
    		
    		if(dist < smallest || smallest < 0){
    			smallest = dist;
    			smallestNode = curr;
    		}
    		curr = curr.getNext();   		
    	}
    	
    	Node oldNext = smallestNode.getNext();
    	smallestNode.setNext(node);
    	node.setNext(oldNext);
    	size++;
    	totalDistance = smallest;
    }
    
    public void insertSmallestByWeightedDistance(Point p){
    	Node node = new Node(p);
    	
    	Node smallestNode = null;
    	Node curr = start;
    	double smallest = -1;
    	double currDistanceTo = -1;
    	double nextDistanceTo = -1;
    	double oldDistance = 0;
    	double distanceSoFar = 0;
    	double unloaded = 0;
    	
    	double currentVal =  weightedDistance;
    	
    	for(int i = 0; i < size; i++){ 
    		
    		currDistanceTo = curr.getPoint().haversineDistanceTo(p);
    		nextDistanceTo = curr.getNext().getPoint().haversineDistanceTo(p);
    		distanceSoFar += oldDistance;
    		oldDistance = curr.getPoint().haversineDistanceTo(curr.getNext().getPoint());
    		
    		//distanceSoFar += oldDistance;
    	
    		//System.out.println(currDistanceTo +":" + nextDistanceTo + ":" + oldDistance + ":" + distanceSoFar + ":" + getWeightFromPointOn(i));
    		double currWeightedDistance = currentVal + p.getGift() * (distanceSoFar + currDistanceTo);

    		if(i>0)
    			unloaded+=curr.getPoint().getGift();
       		currWeightedDistance += (nextDistanceTo+currDistanceTo-oldDistance) * (load-unloaded);
       		
       		if(currWeightedDistance<smallest || smallest == -1){
       			smallest = currWeightedDistance;
       			smallestNode = curr;
       		}
    		curr = curr.getNext();
    	}

    	Node oldNext = smallestNode.getNext();
    	smallestNode.setNext(node);
    	node.setNext(oldNext);
    	size++;
    	weightedDistance = smallest;
    	load += p.getGift();
    }
    
    
    public double getLoad(){
    	return load;
    }
    
    public double getWeightedDistance(){
    	return weightedDistance;
    }
    
    public void draw(){
    	if(start == null)
    		return;
    	
    	Node curr = start;
    	
    	for(int i = 0; i < size; i++){
    		curr.getPoint().draw();
    		curr.getPoint().drawTo(curr.getNext().getPoint());
    		curr = curr.getNext();
    	}

    }
    
    public Point getStart(){
    	return start.getPoint();
    }
    
    public double calculateWeightedDistance(){
    	Node curr = start;
    	double total = 0;
    	double currLoad = load;
    	for(int i = 0; i < size-1; i++){
    		total += currLoad*curr.getPoint().haversineDistanceTo(curr.getNext().getPoint());
    		
    		currLoad -= curr.getNext().getPoint().getGift();
    		curr = curr.getNext();
    	}
    	total += sleighWeight * curr.getPoint().haversineDistanceTo(start.getPoint());
    	return total;
    }
    
    public void drawHaversine(int width, int height){
    	if(start == null)
    		return;
    	
    	Node curr = start;
    	
    	for(int i = 0; i < size; i++){
    		curr.getPoint().drawHaversinePoint(width, height);
    		//curr.getPoint().drawToHaversinePoint(curr.getNext().getPoint(), width, height);
    		curr = curr.getNext();
    	}
    }
    
    public void dump(int tripId){
    	if(start == null)
    		return;
    	
    	Node curr = start;
    	
    	for(int i = 0; i < size; i++){
    		if(i > 0)
    			System.out.println(tripId + "," + curr.getPoint().getId());
    		curr = curr.getNext();
    	}
    }
}
