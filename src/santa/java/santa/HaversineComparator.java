
public class HaversineComparator implements java.util.Comparator<Point> {
	public static final Point northPole = new Point(-1, 90, 0, 0);
	    @Override
	    public int compare(Point o1, Point o2) {
	        // write comparison logic here like below , it's just a sample
	    	
	    	double number1 = o1.haversineDistanceTo(northPole);
	    	double number2 = o2.haversineDistanceTo(northPole);
	    	int compareTo = number2 > number1 ? -1 : (number2 < number1 ? 1 : 0) ;
	    	return compareTo;
	    }
}