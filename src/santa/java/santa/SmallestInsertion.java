/*************************************************************************
 *  YOU DO NOT NEED TO MODIFY THIS FILE
 *
 *  Compilation:  javac SmallestInsertion.java
 *  Execution:    java SmallestInsertion < file.txt
 *  Dependencies: Tour.java Point.java StdIn.java StdDraw.java
 *
 *  Run smallest insertion heuristic for traveling salesperson problem
 *  and plot results.
 *
 *  % java SmallestInsertion < tsp1000.txt
 *
 *************************************************************************/

public class SmallestInsertion {

    public static void main(String[] args) {

        // get dimensions
        //int w = StdIn.readInt();
        //int h = StdIn.readInt();
    	int w = 720;
    	int h = 720;
        StdDraw.setCanvasSize(w, h);
        StdDraw.setXscale(0, w);
        StdDraw.setYscale(0, h);

        // turn on animation mode
        StdDraw.show(0);

        // run smallest insertion heuristic
        Tour tour = new Tour();
        int counter = 0;
        int tourCounter = 1;
        double totalWeightedDistance = 0;
        Point lastPoint = null;
        double firstDistance = -1;
        
        while (!StdIn.isEmpty()) {
        	int id = StdIn.readInt();
            double x = StdIn.readDouble();
            double y = StdIn.readDouble();
            double gift = StdIn.readDouble();
            Point p = new Point(id, x, y, gift);

            if(p.getGift()+tour.getLoad()>tour.getCapacity() || 
            		(lastPoint!=null && lastPoint.haversineDistanceTo(p)>firstDistance && firstDistance != -1) ||
            		(lastPoint!=null && lastPoint.haversineDistanceTo(p)>500) ){
            	tour.drawHaversine(w, h);
                StdDraw.show(0);
            	System.out.println("tour #:" + tourCounter + " Tour WeightedDistance:" + tour.getWeightedDistance() + " # of Gifts:" + tour.size() + " load:" + tour.getLoad());
            	counter+= tour.size()-1;
            	totalWeightedDistance += tour.getWeightedDistance();
            	//tour.dump(tourCounter);
            	tour = new Tour();
            	tourCounter++;
            	firstDistance = -1;
            	lastPoint = null;
            	//break;
            }
            
            tour.insertSmallestByWeightedDistance(p);
            if(firstDistance == -1){
            	firstDistance = tour.getStart().haversineDistanceTo(p);
            }
             //uncomment the 4 lines below to animate
             //StdDraw.clear();
             //tour.draw();
             //StdDraw.show(0);
             //StdDraw.text(100, 0, "" + tour.distance());
             //StdDraw.show(50);
            
            //counter++;
            
            lastPoint = p;
            //if(counter%1000==0)
            //	System.out.println(counter);
        }
        //tour.dump(tourCounter);
        tour.drawHaversine(w, h);
        StdDraw.show(0);
        System.out.println("tour #:" + tourCounter + " Tour WeightedDistance:" + tour.getWeightedDistance() + " # of Gifts:" + tour.size() + " load:" + tour.getLoad());
        counter+= tour.size()-1;
        totalWeightedDistance += tour.getWeightedDistance();
        System.out.println("total # of gifts:" + counter);
        System.out.println("total weighted distance:" + totalWeightedDistance);
        
        // draw to standard draw 
        //StdDraw.show(0);
        
        // print tour to standard output
        //tour.show();
        //StdOut.printf("Tour distance =  %.4f\n", tour.distance());
        //StdOut.printf("Number of points = %d\n", tour.size());
    }
    

}
