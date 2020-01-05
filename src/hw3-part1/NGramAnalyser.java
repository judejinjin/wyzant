import java.io.*;
import java.lang.Exception;
import java.util.*;

// Main class to calculate counts of ngram in a binary file
public class NGramAnalyser{
    private int lengthOfNgram;
    private int lengthOfSlide;
    private HashMap<NGram, Integer> counter;

	//constructor
	//initialize the ngram HashMap counter
    public NGramAnalyser(int n, int slide){
		lengthOfNgram = n;
		lengthOfSlide = slide;
		counter = new HashMap<NGram, Integer>();
    }

	//add a ngram to HashMap counter
    public void add(NGram ngram){
		Integer n = counter.get(ngram);
		if(n != null){
	    	//System.out.println("found:" +ngram.toString());
	    	n = new Integer(n.intValue()+1);
	    	counter.put(ngram, n);
		}else{
	    	//System.out.println("not found:" +ngram.toString());
	    	counter.put(ngram, new Integer(1));
		}
    }

	//sort the NGram counter in decreasing order
	//and print it out into PrintStream
    public void print(PrintStream ps){
		Set<Map.Entry<NGram, Integer>> set = counter.entrySet();
        List<Map.Entry<NGram, Integer>> list = new ArrayList<Map.Entry<NGram, Integer>>(
                set);
        Collections.sort(list, new Comparator<Map.Entry<NGram, Integer>>() {
            public int compare(Map.Entry<NGram, Integer> o1,
                    Map.Entry<NGram, Integer> o2) {
				//when there is a tie in counts
				//compare ngram, smaller first
				if(o2.getValue().compareTo(o1.getValue())==0){
					return o2.getKey().compareTo(o1.getKey());
				}else{
					return o2.getValue().compareTo(o1.getValue());
				}
            }
        });
	
		for (Map.Entry<NGram, Integer> entry : list) {
	    	NGram key = entry.getKey();
	    	Integer value = entry.getValue();
	    	ps.print(key.toString() + "\t" + value.toString() + "\n");
		}

    }

	//slide the ngram buffer
    public static void slide(byte[] ngram, byte[] slide){
		//slide the ngram to the left
		for(int i = 0; i < ngram.length-slide.length; i++){
	    	ngram[i] = ngram[slide.length+i];
		}
		//put the slide into the ngram
		for(int i = 0; i < slide.length; i++){
		    ngram[ngram.length-slide.length+i] = slide[i];
		}
    }

    public static void main(String args[]){
		int lengthOfNgram=0, lengthOfSlide=0;
		String inputFilename, outputFilename;
		if(args.length != 4){
	    	System.out.println("usage: java NGramAnalyser [length of ngram] [length of sliding window] [input file] [output file]");
	    	System.exit(1);
		}

		try{
	    	lengthOfNgram = Integer.parseInt(args[0]);
		}catch(Exception e){
	    	System.out.println("usage: java NGramAnalyser [length of ngram] [length of sliding window] [input file] [output file]");
	    	System.exit(1);
		}

		try{
	    	lengthOfSlide = Integer.parseInt(args[1]);
		}catch(Exception e){
	    	System.out.println("usage: java NGramAnalyser [length of ngram] [length of sliding window] [input file] [output file]");
	    	System.exit(1);
		}

		if(lengthOfNgram > 3 || lengthOfNgram < 0){
	    	System.out.println("length of ngram should be 1-3(inclusive).");
	    	System.exit(1);
		}
	
		if(lengthOfSlide > lengthOfNgram){
	    	System.out.println("length of sliding window should not be longer than length of ngram");
	    	System.exit(1);
		}
		inputFilename = args[2];
		outputFilename = args[3];
	    	FileInputStream is = null;
		PrintStream ps = null;
		try{
		    is = new FileInputStream(inputFilename);
		}catch(Exception e){
		    System.out.println("error opening input file:" + inputFilename);
		    System.exit(1);
		}

		try{
		    ps = new PrintStream(new FileOutputStream(outputFilename));		    
		}catch(Exception e){
		    System.out.println("error opening output file:" + outputFilename);
		    System.exit(1);
		}
		
		try{

		    NGramAnalyser analyser = new NGramAnalyser(lengthOfNgram, lengthOfSlide);
		    //PrintStream os = new PrintStream(new FileOutputStream(inputFilename + "_" + lengthOfNgram + "_" + lengthOfSlide +".csv"));
		    
	    	byte[] ngram = new byte[lengthOfNgram];
	    	byte[] slide = new byte[lengthOfSlide];
	    	int counter = 0;
	    	long ngramCounter = 0;
	    
	    	int read = is.read(ngram);
	    
	    	while(true) {

				counter += read;
				ngramCounter++;
				analyser.add(new NGram(ngram));
		
				//for(int i = 0; i < ngram.length; i++){
		    		//os.print(String.format("%02x", ngram[i]));
				//}
				//os.print("\n");
				//System.out.println(new NGram(ngram).toString());

				read = is.read(slide);
				if ( read < lengthOfSlide ) {
		    		break;
				}

				slide(ngram, slide);
	    	}

	    	System.out.println(counter + " bytes read.");
	    	System.out.println(ngramCounter + " " + lengthOfNgram + "-gram read.");
	    	//os.close();
	    	is.close();

	    	analyser.print(ps);
	    	ps.close();
		}catch(Exception e){
		    e.printStackTrace();
		    System.out.println("usage: java NGramAnalyser [length of ngram] [length of sliding window] [input file] [output file]");
		}

	}
}
