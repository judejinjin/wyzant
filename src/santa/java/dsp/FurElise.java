
/**********************************************************************************
 *
 *  Compilation:  javac -cp .:stdplayer.jar FurElise.java
 *  Execution:    java  -cp .:stdplayer.jar FurElise
 *  Dependencies: Wave.java
 *
 *  Note:  under Windows, replace the : with a ; when specifying the classpath.
 *  Note:  the file stdplayer.jar must be in the current directory
 *
 *  Generates the first 9 notes of Fur Elise using the Wave class and plays them.
 *
 **********************************************************************************/

import javazoom.jl.player.StdPlayer;

public class FurElise {
    public static void main(String[] args) {
        double amplitude = 0.8;
        StdPlayer.open();

        // Setup the notes - multiply freq by 2 to raise an octave
        Wave E  = new Wave(659.25 * 2, 0.3, amplitude);
        Wave DS = new Wave(622.25 * 2, 0.3, amplitude);
        Wave B  = new Wave(493.88 * 2, 0.3, amplitude);
        Wave D  = new Wave(587.33 * 2, 0.3, amplitude);
        Wave C  = new Wave(523.25 * 2, 0.3, amplitude);
        Wave A  = new Wave(440.00 * 2, 0.9, amplitude);

        // e, d#, e, d#, e, b, d, c, a
        E.play();
        DS.play();
        E.play();
        DS.play();
        E.play();
        B.play();
        D.play();
        C.play();
        A.play();

        StdPlayer.close();
        System.exit(0);
    }
}
