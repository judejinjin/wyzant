
/**********************************************************************************
 *
 *  Compilation:  javac -classpath .:stdplayer.jar Scale.java Wave.java
 *  Execution:    java  -classpath .:stdplayer.jar Scale note
 *
 *  Note:  under Windows, replace the : with a ; when specifying the classpath.
 *  Note:  the file player.jar must be in the current directory
 *
 *  Plays a major scale. Legal values are listed below
 *
 *             A, A#, B, C,C#, D, D#, E, F, F#, G, G#
 *
 *
 *  % java -cp .:stdplayer.jar Scale A
 *  % java -cp .:stdplayer.jar Scale A#
 *  % java -cp .:stdplayer.jar Scale C#
 *
 **********************************************************************************/

import javazoom.jl.player.StdPlayer;

public class Scale {
    public static void main(String[] args) {
        int[] majorScale = { 0, 2, 4, 5, 7, 9, 11, 12 };     // scale
        double octave    = 1.0;                              // octave multiplier (a power of 2)
        double seconds   = 1.0;                              // duration of each note
        double amplitude = 0.8;                              // volume

        // two octaves worth
        double[] freqs = {  440.00,  466.16,  493.88,  523.25,  554.37,  587.33,
                            622.25,  659.26,  698.46,  739.99,  783.99,  830.61,
                            880.00,  932.32,  987.76, 1046.50, 1108.74, 1174.66,
                           1244.50, 1318.52, 1396.92, 1479.98, 1567.98, 1661.22  };

        String[] notes = {  "A",    "A#",   "B",    "C",    "C#",   "D",
                            "D#",   "E",    "F",    "F#",   "G",    "G#"  };


        // find the index of the first note
        int offset;
        for (offset = 0; offset < notes.length; offset++)
           if (notes[offset].equals(args[0].toUpperCase())) break;
        if (offset == notes.length) {
           System.out.println("Invalid scale name");
           System.exit(0);
        }

        StdPlayer.open();

        // play the scale up
        for (int i = 0; i < 8; i++) {
            double freq = freqs[offset + majorScale[i]];
            Wave w = new Wave(freq * octave, seconds, amplitude);
            w.play();
        }

        // play the scale down
        for (int i = 7; i >= 0; i--) {
            double freq = freqs[offset + majorScale[i]];
            Wave w = new Wave(freq * octave, seconds, amplitude);
            w.play();
        }

        StdPlayer.close();
        System.exit(0);
    }
}
