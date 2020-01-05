public class NGram{

	//use 4 byte integer to store n-gram where n <= 3
    private int ngram;
	//length of n-gram
    private int n;

    public NGram(byte[] b){
		n = b.length;
		ngram = 0;
		//shift the bits left by one byte when i > 0
		for(int i = 0; i < n; i++){
	 	   ngram = ngram << (8*(i==0?0:1)) | (b[i] & 0xFF);
		}
    }
    
    public int n(){
	return n;
    }

    @Override
    public int hashCode(){
	return ngram;
    }

    @Override
    public boolean equals(Object b){
		if(b instanceof NGram){
	    	return (n()==((NGram)b).n()) && (hashCode() == ((NGram)b).hashCode());
		}else{
	    	return false;
		}
    }

	//special compareTo logic to enable sorting in decreasing order of ngram
    public int compareTo(NGram b){
		if(hashCode() > b.hashCode())
	    	return -1;
		else if(hashCode() == b.hashCode()){
	    	return 0;
		}else
	    	return 1;
    }
    
    @Override
    public String toString(){
		String s = "";
		for(int i = 0; i < n; i++){
	    	if(i== 0)
		    				s = String.format("%02x", (ngram >> (n-1-i)*8 & 0xFF) );
	    	else
		    				s += " " + String.format("%02x", (ngram >> (n-1-i)*8 & 0xFF) );
		}
		return s;
    }
}
