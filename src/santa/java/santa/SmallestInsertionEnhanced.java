import java.util.Arrays;

/*************************************************************************
 * YOU DO NOT NEED TO MODIFY THIS FILE
 *
 * Compilation: javac SmallestInsertion.java Execution: java SmallestInsertion <
 * file.txt Dependencies: Tour.java Point.java StdIn.java StdDraw.java
 *
 * Run smallest insertion heuristic for traveling salesperson problem and plot
 * results.
 *
 * % java SmallestInsertion < tsp1000.txt
 *
 *************************************************************************/

public class SmallestInsertionEnhanced {

	public static int findNext(Point[] gifts) {
		for (int i = 0; i < gifts.length; i++) {
			if (!gifts[i].isDelivered())
				return i;
		}
		return -1;
	}

	public static void shuffleArray(Point[] ar) {
		// If running on Java 6 or older, use `new Random()` on RHS here
		java.util.Random rnd = new java.util.Random(System.currentTimeMillis());
		for (int i = ar.length - 1; i > 0; i--) {
			int index = rnd.nextInt(i + 1);
			// Simple swap
			Point a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
	}
	

	public static void main(String[] args) {

		// get dimensions
		// int w = StdIn.readInt();
		// int h = StdIn.readInt();
		int w = 720;
		int h = 720;
		StdDraw.setCanvasSize(w, h);
		StdDraw.setXscale(0, w);
		StdDraw.setYscale(0, h);

		// turn on animation mode
		StdDraw.show(0);

		// run smallest insertion heuristic

		Point[] gifts = new Point[100000];
		int i = 0;

		while (!StdIn.isEmpty()) {
			int id = StdIn.readInt();
			double x = StdIn.readDouble();
			double y = StdIn.readDouble();
			double gift = StdIn.readDouble();
			Point p = new Point(id, x, y, gift);
			gifts[i++] = p;
		}
		
		/*
		java.util.List<Point> list = Arrays.asList(gifts);
		
		java.util.Collections.sort(list, new HaversineComparator());

		gifts = (Point[])list.toArray();
		*/
		
		int counter = 0;
		int tourCounter = 1;
		double totalWeightedDistance = 0;

		Point lastPoint = null;
		double firstDistance = -1;
		Point firstPoint = null;
		double searchRadius = 6371;

		int doneIndex = 0;
		Tour tour = new Tour();
		
		while ((doneIndex = findNext(gifts)) >= 0) {
			int tripCounter = 0;
			for (int j = doneIndex; j < gifts.length; j++) {
				if (gifts[j].isDelivered())
					continue;
				
				
				if (firstDistance != -1)
					searchRadius = firstDistance * tour.getStart().getGift() / 50;
				
				if(searchRadius > 500)
					searchRadius = 500;
				
				
				if (gifts[j].getGift() + tour.getLoad() < tour.getCapacity()
						&& ((lastPoint != null
						&& lastPoint.haversineDistanceTo(gifts[j]) > firstDistance && firstDistance != -1)
						|| (firstPoint != null && firstPoint.haversineDistanceTo(gifts[j]) > searchRadius))
						) {
					continue;
				} else {
					if (gifts[j].getGift() + tour.getLoad() > tour.getCapacity()) {
						break;
						/*
						tour.drawHaversine(w, h);
						StdDraw.show(0);
						System.out.println(
								"tour #:" + tourCounter + " Tour WeightedDistance:" + tour.getWeightedDistance()
										+ " # of Gifts:" + tour.size() + " load:" + tour.getLoad());
						counter += tour.size() - 1;
						totalWeightedDistance += tour.getWeightedDistance();
						// tour.dump(tourCounter);
						tour = new Tour();
						tourCounter++;
						firstDistance = -1;
						lastPoint = null;
						firstPoint = null;
						*/
					}
				}

				tour.insertSmallestByWeightedDistance(gifts[j]);
				if (firstDistance == -1) {
					firstDistance = tour.getStart().haversineDistanceTo(gifts[j]);
				}

				if (firstPoint == null) {
					firstPoint = gifts[j];
				}

				lastPoint = gifts[j];
				gifts[j].setDelivered();
				tripCounter++;
				//if(tripCounter >= 67)
				//	break;
			}
			tour.drawHaversine(w, h);
			StdDraw.show(0);
			System.out.println("tour #:" + tourCounter + " Tour WeightedDistance:" + tour.getWeightedDistance()
					+ " # of Gifts:" + tour.size() + " load:" + tour.getLoad());
			counter += tour.size() - 1;
			totalWeightedDistance += tour.getWeightedDistance();
			// tour.dump(tourCounter);
			tour = new Tour();
			tourCounter++;
			firstDistance = -1;
			lastPoint = null;
			firstPoint = null;
		}

		// tour.dump(tourCounter);
		tour.drawHaversine(w, h);
		StdDraw.show(0);
		System.out.println("tour #:" + tourCounter + " Tour WeightedDistance:" + tour.getWeightedDistance()
				+ " # of Gifts:" + tour.size() + " load:" + tour.getLoad());
		counter += tour.size() - 1;
		totalWeightedDistance += tour.getWeightedDistance();
		System.out.println("total # of gifts:" + counter);
		System.out.println("total weighted distance:" + totalWeightedDistance);

	}

}
