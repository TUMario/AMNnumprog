import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Die Klasse LinearInterpolation beschreibt ein stueckweises
 * Interpolationsverfahren. In jedem Intervall zwischen zwei benachbarten
 * Stuetzstellen werden die entsprechenden Stuetzwerte mit einer Geraden
 * verbunden. Liegt eine zu auswertende Stelle z ausserhalb der Stuetzgrenzen,
 * werden die aeussersten Werte y[0] bzw. y[n] zurueckgegeben.
 *
 * @author braeckle
 *
 */
public class LinearInterpolation implements InterpolationMethod {

    /** Die Stuetzstellen x_i */
    double[] x;
    /** Die Stuetzwerte y_i */
    double[] y;

    @Override
    public void init(double a, double b, int n, double[] y) {
        this.y = y;
        x = new double[n + 1];
        double h = (b - a) / n;

        for (int i = 0; i < n + 1; i++) {
            x[i] = a + i * h;
        }
    }

    /**
     * Initialisierung des Interpolationsverfahrens mit beliebigen Stuetzstellen
     * und ordnet die Stuetzstellen der Groesse nach. Die Faelle
     * "x und y sind unterschiedlich lang" oder "eines der beiden Arrays ist
     * leer" werden nicht beachtet.
     *
     * @param x
     *            Stuetzstellen
     * @param y
     *            Stuetzwerte
     */
    public void init(final double[] x, double[] y) {
        if (x.length != y.length || x.length == 0)
            return;

        int n = x.length;

        /* die Stuetzstellen werden der Groesse nach geordnet */
        List<Integer> indices = new ArrayList<Integer>(n);
        for (int i = 0; i < n; i++)
            indices.add(i);

        Comparator<Integer> comparator = new Comparator<Integer>() {

            @Override
            public int compare(Integer i, Integer j) {
                if (x[i] - x[j] < 0)
                    return -1;
                if (x[i] - x[j] > 0)
                    return 1;
                else
                    return 0;
            }
        };

        Collections.sort(indices, comparator);

        this.x = new double[n];
        this.y = new double[n];

        for (int i = 0; i < n; i++) {
            int index = indices.get(i);
            this.x[i] = x[index];
            this.y[i] = y[index];
        }
    }

    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird eine Gerade mit den Punkten (x_i,y_i)
     * und (x_i+1, y_i+1) gebildet und in z ausgewertet.
     *
     * Die Stuetzstellen liegen der Groesse nach geordnet vor. Der Fall
     * ungeordneter Stuetzstellen oder leerer Stuetzstellen muss nicht extra
     * behandelt werden.
     */
    @Override
    public double evaluate(double z) {
        if (z <= x[0]) return y[0];
        if (z >= x[x.length-1]) return y[y.length-1];
        int binSearchUp = (x.length-1);
        int binSearchLow = 0;
        int pivotIndex = 0;
        double[] p1xy = {0.0,0.0};
        double[] p2xy = {0.0,0.0};
        while (binSearchUp-binSearchLow > 0) {
            pivotIndex = (binSearchUp+binSearchLow)/2;
            if (x[pivotIndex] == z) {
                    return y[pivotIndex];
            }
            if (x[pivotIndex] < z) {
                    if (x[pivotIndex + 1] > z) {
                        p1xy[0] = x[pivotIndex]; /// [0] = x1, [1] = y1
                        p1xy[1] = y[pivotIndex];
                        p2xy[0] = x[pivotIndex + 1];
                        p2xy[1] = y[pivotIndex + 1];
                        break;
                    }
                    binSearchLow = binSearchLow + (binSearchUp-binSearchLow)/2 + 1;
            }
            if (x[pivotIndex] > z) {
                if (x[pivotIndex-1] < z) {
                    p1xy[0] = x[pivotIndex-1];
                    p1xy[1] = y[pivotIndex-1];
                    p2xy[0] = x[pivotIndex];
                    p2xy[1] = y[pivotIndex];
                    break;
                }
                binSearchUp = binSearchUp - (binSearchUp-binSearchLow)/2 + 1;
            }
        }
        return p1xy[1] + ((p2xy[1]-p1xy[1]) / (p2xy[0]-p1xy[0])) * (z - p1xy[0]); //formel Lin Interp.
    }

}
