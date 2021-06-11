
// package essential;
import java.util.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.TargetDataLine;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

/**
 * A utility class provides general functions for recording sound.
 * 
 * @author www.codejava.net
 *
 */
public class SoundRecordingUtil {
    private static final int BUFFER_SIZE = 4096;
    private ByteArrayOutputStream recordBytes;
    private TargetDataLine audioLine;
    private AudioFormat format;

    private boolean isRunning;

    public final int UPPER_LIMIT = 300;
    public final int LOWER_LIMIT = 40;

    public final int[] RANGE = new int[] { 40, 80, 120, 180, UPPER_LIMIT + 1 };

    // Find out in which range
    public int getIndex(int freq) {
        int i = 0;
        while (RANGE[i] < freq)
            i++;
        return i;
    }

    /**
     * Defines a default audio format used to record
     */
    AudioFormat getAudioFormat() {
        float sampleRate = 44100;
        int sampleSizeInBits = 8;
        int channels = 1;
        boolean signed = true;
        boolean bigEndian = true;
        return new AudioFormat(sampleRate, sampleSizeInBits, channels, signed, bigEndian);
    }

    class Complex {
        private final double re; // the real part
        private final double im; // the imaginary part

        // create a new object with the given real and imaginary parts
        public Complex(double real, double imag) {
            re = real;
            im = imag;
        }

        // return a string representation of the invoking Complex object
        public String toString() {
            if (im == 0)
                return re + "";
            if (re == 0)
                return im + "i";
            if (im < 0)
                return re + " - " + (-im) + "i";
            return re + " + " + im + "i";
        }

        // return abs/modulus/magnitude
        public double abs() {
            return Math.hypot(re, im);
        }

        // return angle/phase/argument, normalized to be between -pi and pi
        public double phase() {
            return Math.atan2(im, re);
        }

        // return a new Complex object whose value is (this + b)
        public Complex plus(Complex b) {
            Complex a = this; // invoking object
            double real = a.re + b.re;
            double imag = a.im + b.im;
            return new Complex(real, imag);
        }

        // return a new Complex object whose value is (this - b)
        public Complex minus(Complex b) {
            Complex a = this;
            double real = a.re - b.re;
            double imag = a.im - b.im;
            return new Complex(real, imag);
        }

        // return a new Complex object whose value is (this * b)
        public Complex times(Complex b) {
            Complex a = this;
            double real = a.re * b.re - a.im * b.im;
            double imag = a.re * b.im + a.im * b.re;
            return new Complex(real, imag);
        }

        // return a new object whose value is (this * alpha)
        public Complex scale(double alpha) {
            return new Complex(alpha * re, alpha * im);
        }

        // return a new Complex object whose value is the conjugate of this
        public Complex conjugate() {
            return new Complex(re, -im);
        }

        // return a new Complex object whose value is the reciprocal of this
        public Complex reciprocal() {
            double scale = re * re + im * im;
            return new Complex(re / scale, -im / scale);
        }

        // return the real or imaginary part
        public double re() {
            return re;
        }

        public double im() {
            return im;
        }

        // return a / b
        public Complex divides(Complex b) {
            Complex a = this;
            return a.times(b.reciprocal());
        }

        // return a new Complex object whose value is the complex exponential of this
        public Complex exp() {
            return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
        }

        // return a new Complex object whose value is the complex sine of this
        public Complex sin() {
            return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
        }

        // return a new Complex object whose value is the complex cosine of this
        public Complex cos() {
            return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re) * Math.sinh(im));
        }

        // return a new Complex object whose value is the complex tangent of this
        public Complex tan() {
            return sin().divides(cos());
        }

        // a static version of plus
        public Complex plus(Complex a, Complex b) {
            double real = a.re + b.re;
            double imag = a.im + b.im;
            Complex sum = new Complex(real, imag);
            return sum;
        }

        // See Section 3.3.
        public boolean equals(Object x) {
            if (x == null)
                return false;
            if (this.getClass() != x.getClass())
                return false;
            Complex that = (Complex) x;
            return (this.re == that.re) && (this.im == that.im);
        }

        // See Section 3.3.
        public int hashCode() {
            return Objects.hash(re, im);
        }
    }

    class FFT {

        // compute the FFT of x[], assuming its length n is a power of 2
        public Complex[] fft(Complex[] x) {
            int n = x.length;

            // base case
            if (n == 1)
                return new Complex[] { x[0] };

            // radix 2 Cooley-Tukey FFT
            if (n % 2 != 0) {
                throw new IllegalArgumentException("n is not a power of 2");
            }

            // compute FFT of even terms
            Complex[] even = new Complex[n / 2];
            for (int k = 0; k < n / 2; k++) {
                even[k] = x[2 * k];
            }
            Complex[] evenFFT = fft(even);

            // compute FFT of odd terms
            Complex[] odd = even; // reuse the array (to avoid n log n space)
            for (int k = 0; k < n / 2; k++) {
                odd[k] = x[2 * k + 1];
            }
            Complex[] oddFFT = fft(odd);

            // combine
            Complex[] y = new Complex[n];
            for (int k = 0; k < n / 2; k++) {
                double kth = -2 * k * Math.PI / n;
                Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
                y[k] = evenFFT[k].plus(wk.times(oddFFT[k]));
                y[k + n / 2] = evenFFT[k].minus(wk.times(oddFFT[k]));
            }
            return y;
        }
    }

    /**
     * Start recording sound.
     * 
     * @throws LineUnavailableException if the system does not support the specified
     *                                  audio format nor open the audio data line.
     */
    public void start() throws LineUnavailableException {
        format = getAudioFormat();
        DataLine.Info info = new DataLine.Info(TargetDataLine.class, format);

        // checks if system supports the data line
        if (!AudioSystem.isLineSupported(info)) {
            throw new LineUnavailableException("The system does not support the specified format.");
        }

        audioLine = AudioSystem.getTargetDataLine(format);

        audioLine.open(format);
        audioLine.start();

        byte[] buffer = new byte[BUFFER_SIZE];
        int bytesRead = 0;

        recordBytes = new ByteArrayOutputStream();// out in toptal
        isRunning = true;

        while (isRunning) {
            bytesRead = audioLine.read(buffer, 0, buffer.length);
            recordBytes.write(buffer, 0, bytesRead);
            System.out.print(recordBytes);
        }
        // return recordBytes;
    }

    void getHash(ByteArrayOutputStream out) {
        byte audio[] = out.toByteArray();

        final int totalSize = audio.length;

        int amountPossible = totalSize / 4096;

        // When turning into frequency domain we'll need complex numbers:
        Complex[][] results = new Complex[amountPossible][];

        // For all the chunks:
        for (int times = 0; times < amountPossible; times++) {
            Complex[] complex = new Complex[4096];
            for (int i = 0; i < 4096; i++) {
                // Put the time domain data into a complex number with imaginary
                // part as 0:
                complex[i] = new Complex(audio[(times * 4096) + i], 0);
            }
            // Perform FFT analysis on the chunk:
            FFT f = new FFT();
            results[times] = f.fft(complex);
        }

        double[][] highscores = new double[results.length][5];
        for (int i = 0; i < results.length; i++) {
            for (int j = 0; j < 5; j++) {
                highscores[i][j] = 0;
            }
        }

        double[][] recordPoints = new double[results.length][UPPER_LIMIT];
        for (int i = 0; i < results.length; i++) {
            for (int j = 0; j < UPPER_LIMIT; j++) {
                recordPoints[i][j] = 0;
            }
        }

        long[][] points = new long[results.length][5];
        for (int i = 0; i < results.length; i++) {
            for (int j = 0; j < 5; j++) {
                points[i][j] = 0;
            }
        }

        for (int t = 0; t < results.length; t++) {
            for (int freq = LOWER_LIMIT; freq < UPPER_LIMIT - 1; freq++) {
                // Get the magnitude:
                double mag = Math.log(results[t][freq].abs() + 1);

                // Find out which range we are in:
                int index = getIndex(freq);

                // Save the highest magnitude and corresponding frequency:
                if (mag > highscores[t][index]) {
                    highscores[t][index] = mag;
                    recordPoints[t][freq] = 1;
                    points[t][index] = freq;
                }
            }

            // for (int k = 0; k < 5; k++) {
            //     System.out.println("" + highscores[t][k] + ";" + recordPoints[t][k] + "\t");
            // }
            long h = hash(points[t][0], points[t][1], points[t][2], points[t][3]);
            System.out.print(h);
        }

        // outFile.write("\n");

    }

    private static final int FUZ_FACTOR = 2;

    private long hash(long p1, long p2, long p3, long p4) {
        return (p4 - (p4 % FUZ_FACTOR)) * 100000000 + (p3 - (p3 % FUZ_FACTOR)) * 100000 + (p2 - (p2 % FUZ_FACTOR)) * 100
                + (p1 - (p1 % FUZ_FACTOR));
    }

    // public void getHash() {
    // int chunkSize = 4096;

    // byte audio[] = recordBytes.toByteArray();
    // int totalSize = audio.length;
    // int sampledChunkSize = totalSize / chunkSize;
    // Complex[][] result = new Complex[sampledChunkSize][];

    // double[][] highscores = new double[350][];
    // int[][] points = new int[350][];

    // // SoundRecordingUtil c = new SoundRecordingUtil();
    // // Complex c1 = new Complex();

    // for (int j = 0; j < sampledChunkSize; j++) {
    // Complex[] complexArray = new Complex[chunkSize];

    // for (int i = 0; i < chunkSize; i++) {
    // complexArray[i] = new Complex(audio[(j * chunkSize) + i], 0);
    // }
    // FFT f = new FFT();
    // result[j] = f.fft(complexArray);
    // }
    // final int[] RANGE = new int[] { 40, 80, 120, 180, 300};

    // // find out in which range is frequency
    // final int getIndex(int freq) {
    // int i = 0;
    // while (RANGE[i] < freq)
    // i++;
    // return i;
    // }

    // // result is complex matrix obtained in previous step
    // for(int t = 0;t<result.length;t++)
    // {
    // for (int freq = 40; freq < 300; freq++) {
    // // Get the magnitude:
    // double mag = Math.log(result[t][freq].abs() + 1);

    // // Find out which range we are in:
    // int index = getIndex(freq);

    // // Save the highest magnitude and corresponding frequency:
    // if (mag > highscores[t][index]) {
    // points[t][index] = freq;
    // }
    // }
    // }
    // // form hash tag
    // long h = hash(points[t][0], points[t][1], points[t][2], points[t][3]);
    // }

    // private static final int FUZ_FACTOR = 2;

    // private long hash(long p1, long p2, long p3, long p4) {
    // return (p4 - (p4 % FUZ_FACTOR)) * 100000000 + (p3 - (p3 % FUZ_FACTOR))
    // * 100000 + (p2 - (p2 % FUZ_FACTOR)) * 100
    // + (p1 - (p1 % FUZ_FACTOR));
    // }

    // }}

    // int blockSizeY = 10;
    // int blockSizeX = 10;
    // boolean logModeEnabled = false;
    // int size = 50;

    // @Override
    // protected void paintComponent(Graphics g) {
    // super.paintComponent(g);

    // Graphics2D g2d = (Graphics2D) g;
    // }
    // for (int i = 0; i < result.length; i++) {
    // int freq = 1;
    // for (int line = 1; line < size; line++) {
    // // To get the magnitude of the sound at a given frequency slice
    // // get the abs() from the complex number.
    // // In this case I use Math.log to get a more managable number (used for
    // color)
    // double magnitude = Math.log(results[i][freq].abs() + 1);

    // // The more blue in the color the more intensity for a given frequency point:
    // g2d.setColor(new Color(0, (int) magnitude * 10, (int) magnitude * 20));
    // // Fill:
    // g2d.fillRect(i * blockSizeX, (size - line) * blockSizeY, blockSizeX,
    // blockSizeY);

    // // I used a improviced logarithmic scale and normal scale:
    // if (logModeEnabled && (Math.log10(line) * Math.log10(line)) > 1) {
    // freq += (int) (Math.log10(line) * Math.log10(line));
    // } else {
    // freq++;
    // }
    // }
    // }

    /**
     * Stop recording sound.
     * 
     * @throws IOException if any I/O error occurs.
     */
    public void stop() throws IOException {
        isRunning = false;

        if (audioLine != null) {
            audioLine.drain();
            audioLine.close();
        }
    }

    /**
     * Save recorded sound data into a .wav file format.
     * 
     * @param wavFile The file to be saved.
     * @throws IOException if any I/O error occurs.
     */
    public void save(File wavFile) throws IOException {
        byte[] audioData = recordBytes.toByteArray();
        for(int i=0; i<audioData.length; i++) {
            System.out.print(audioData[i]);
        }
        ByteArrayInputStream bais = new ByteArrayInputStream(audioData);
        AudioInputStream audioInputStream = new AudioInputStream(bais, format,
                audioData.length / format.getFrameSize());

        AudioSystem.write(audioInputStream, AudioFileFormat.Type.WAVE, wavFile);
        getHash(recordBytes);

        audioInputStream.close();
        recordBytes.close();

        // return wavFile;
    }

    public static void main(String[] args) {
        final int RECORD_TIME = 10000;
        // ByteArrayOutputStream recordBytes = new ByteArrayOutputStream();
        File wavFile = new File("/Users/divita/Downloads/audio/RecordAudio.wav");
        wavFile.setWritable(true);

        final SoundRecordingUtil recorder = new SoundRecordingUtil();

        // create a separate thread for recording
        Thread recordThread = new Thread(new Runnable() {
            @Override
            public void run() {
                try {
                    System.out.println("Start recording...");
                    // final ByteArrayOutputStream recordBytes = new ByteArrayOutputStream();
                    recorder.start();
                    // recorder.getHash(recordBytes);
                } catch (LineUnavailableException ex) {
                    ex.printStackTrace();
                    System.exit(-1);
                }
            }
        });

        recordThread.start();

        try {
            Thread.sleep(RECORD_TIME);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }

        try {
            // ByteArrayOutputStream recordBytes = new ByteArrayOutputStream();
            recorder.stop();
            recorder.save(wavFile);
            // ByteArrayOutputStream recordBytes = recorder.start();
            // recorder.getHash(recordBytes);
            System.out.println("STOPPED");
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        System.out.println("DONE");
    }
}