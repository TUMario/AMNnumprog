package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
        int n = c.length;
        Complex[] retComp = new Complex[c.length];
        if (n == 1) {
            retComp[0]=c[0];
        }
        else {
            int m = n>>1;
            Complex[] placeholder1 = new Complex[m];
            Complex[] placeholder2 = new Complex[m];
            for (int i = 0; i<=(n-2) ; i+=2) {
                placeholder1[i/2] = c[i];
            }
            for (int i = 1; i<=(n-1) ; i+=2) {
                placeholder2[i/2] = c[i];
            }
            Complex[] z1 = ifft(placeholder1);
            Complex[] z2 = ifft(placeholder2);
            Complex omega = Complex.fromPolar(1, 2 * Math.PI / n);
            for (int j = 0; j <= m-1; j++) {
                retComp[j] = z1[j].add((omega.power(j)).mul(z2[j])) ;
                retComp[m+j] = z1[j].sub((omega.power(j)).mul(z2[j]));
            }
        }
        return retComp;
    }
}
