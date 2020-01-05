/*************************************************************************
 *  Compilation:  javac Harvester.java
 *  Execution:    java Harvester pattern url
 *  Dependencies: In.java
 *  
 *  Downloads the web page and prints out all strings that
 *  match the specified pattern.
 *
 *  % java Harvester "[a-z]+@([a-z]+\.)+(edu|com|biz|tv)"  http://www.cs.princeton.edu/courses/archive/fall06/cos126/people.html
 * doug@cs.princeton.edu
 * nadiah@cs.princeton.edu
 * ajfeldma@cs.princeton.edu
 * gareis@princeton.edu
 * astoler@princeton.edu
 * maia@cs.princeton.edu
 * elliottk@cs.princeton.edu
 * vivek@cs.princeton.edu
 * troconno@cs.princeton.edu
 *
 *  %  java Harvester "GCG(CGG|AGG)*CTG" chromosomeX.txt
 *  GCGCGGAGGCGGCGGCGGCGGCGGAGGCGGCGGCTG
 *  GCGCGGCGGAGGCTG
 *
 *************************************************************************/

import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class Harvester { 

    public static void main(String[] args) { 
        String regexp = args[0];
        In in = new In(args[1]);
        String input = in.readAll();

        Pattern pattern = Pattern.compile(regexp);
        Matcher matcher = pattern.matcher(input);
    
        while (matcher.find()) {
            String s = matcher.group();
            System.out.println(s);
        }

    
   }

}
