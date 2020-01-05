
/**********************************************************************************
 *
 *  Compilation:  javac -cp .:stdplayer.jar A.java
 *  Execution:    java  -cp .:stdplayer.jar A
 *  Dependencies: Wave.java
 *
 *  Note:  under Windows, replace the : with a ; when specifying the classpath.
 *  Note:  the file stdplayer.jar must be in the current directory
 *
 * Generates and plays an A for two seconds.
 *
 **********************************************************************************/

import javazoom.jl.player.StdPlayer;

public class A {
    public static void main(String[] args) {
        StdPlayer.open();
        Wave A = new Wave(440.0, 2.0, .8);
        A.play();
        StdPlayer.close();
        System.exit(0);
    }
}
