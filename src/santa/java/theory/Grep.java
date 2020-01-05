/*************************************************************************
 *  Compilation:  javac Grep.java
 *  Execution:    java Grep pattern file
 *
 *  % java Grep "..oo..oo." words.txt
 *  bloodroot
 *  schoolbook
 *  schoolroom
 * 
 *  %  java Grep "a.*e.*i.*o.*u" words.txt
 *  adventitious
 *  facetious
 *  sacrilegious

 *
 *************************************************************************/

public class Grep { 

    public static void main(String[] args) { 
        String re = ".*" + args[0] + ".*";
        In in = new In(args[1]);
        while (!in.isEmpty()) {
            String line = in.readLine();
            if (line.matches(re)) System.out.println(line);
        }
   }

}
