package asg.ion;

import asg.cliche.Command;
import asg.cliche.Param;
import asg.cliche.ShellFactory;
import geopt.geometry.Vector;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.ZoomableChart;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.views.ChartPanel;
import java.awt.Dimension;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Scanner;
import javax.swing.JFrame;

import static java.lang.Math.*;

/**
 *
 * @author asg
 */
public class TrajectorySolver {

    private int dimX;
    private int dimY;
    private int dimZ;

    private double[][][] potential;
    private double[][][] eX;
    private double[][][] eY;
    private double[][][] eZ;

    private void setDimensions(int x, int y, int z) {
        this.dimX = x;
        this.dimY = y;
        this.dimZ = z;
        potential = new double[x][y][z];
        eX = new double[x][y][z];
        eY = new double[x][y][z];
        eZ = new double[x][y][z];
    }

    @Command
    public void loadField(
            @Param(name="dimX") int dimX,
            @Param(name="dimY") int dimY,
            @Param(name="dimZ") int dimZ,
            @Param(name="file-name") String fileName) throws IOException {

        setDimensions(dimX, dimY, dimZ);

        Scanner s = new Scanner(new File(fileName));

        int x, y, z;
        boolean electrode;
        double value;

        while (s.hasNext()) {
            x = s.nextInt();
            y = s.nextInt();
            z = s.nextInt();
            electrode = s.nextInt() == 1;
            value = s.nextDouble();

            potential[x][y][z] = electrode ? Double.NaN : value;
        }
    }

    // ix < dimX - 1, etc.
    private double dx(int ix, int iy, int iz) {
        return (potential[ix + 1][iy][iz] - potential[ix][iy][iz]);
    }
    private double dy(int ix, int iy, int iz) {
        return (potential[ix][iy + 1][iz] - potential[ix][iy][iz]);
    }
    private double dz(int ix, int iy, int iz) {
        return (potential[ix][iy][iz + 1] - potential[ix][iy][iz]);
    }

    @Command
    public void computeE() {
        for (int ix = 0; ix < dimX - 1; ix++) {
            for (int iy = 0; iy < dimY - 1; iy++) {
                for (int iz = 0; iz < dimZ - 1; iz++) {
                    eX[ix][iy][iz] = -dx(ix, iy, iz);
                    eY[ix][iy][iz] = -dy(ix, iy, iz);
                    eZ[ix][iy][iz] = -dz(ix, iy, iz);
                }
                eX[ix][iy][dimZ - 1] = dx(ix, iy, dimZ - 2);
            }
            for (int iz = 0; iz < dimZ - 1; iz++) {
                eX[ix][dimY - 1][iz] = dx(ix, dimY - 2, iz);
            }
        }
        for (int iy = 0; iy < dimY - 1; iy++) {
            for (int iz = 0; iz < dimZ - 1; iz++) {
                eX[dimX - 1][iy][iz] = dx(dimX - 2, iy, iz);
            }
        }
    }

    public static interface Callback {
        void params(double dt, double m, double e,
                double x0, double y0, double z0,
                double kineticE, double dirX, double dirY, double dirZ);
        void point(double t, double x, double y, double z,
                double vx, double vy, double vz);
        void close();
    }

    private static final double interpolate(double x1, double x2, double y1, double y2, double x) {
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
    }

    private static final int X = 0;
    private static final int Y = 1;
    private static final int Z = 2;


    // returns null if hitting an electrode
    @Command
    public final double[] getIntensity(double x, double y, double z) {
        final int ix = (int) Math.round(x);
        final int iy = (int) Math.round(y);
        final int iz = (int) Math.round(z);

        if (Double.isNaN(potential[ix][iy][iz])) {
            return null;
        }

        return new double[] {
            interpolate(eX, x, y, z),
            interpolate(eY, x, y, z),
            interpolate(eZ, x, y, z)
        };
    }

    // lattice with space 1
    private static final double interpolate(double[][][] v, double x, double y, double z) {
        // http://en.wikipedia.org/wiki/Trilinear_interpolation
        final int
                xF = (int)floor(x), xC = (int)ceil(x),
                yF = (int)floor(y), yC = (int)ceil(y),
                zF = (int)floor(z), zC = (int)ceil(z);
        final double
                xd = x - xF,
                yd = y - yF,
                zd = z - zF;
        final double
                i1 = v[xF][yF][zF] * (1 - zd) + v[xF][yF][zC] * zd,
                i2 = v[xF][yC][zF] * (1 - zd) + v[xF][yC][zC] * zd,
                j1 = v[xC][yF][zF] * (1 - zd) + v[xF][yF][zC] * zd,
                j2 = v[xC][yC][zF] * (1 - zd) + v[xF][yC][zC] * zd;
        final double
                w1 = i1 * (1 - yd) + i2 * yd,
                w2 = j1 * (1 - yd) + j2 * yd;

        return w1 * (1 - xd) + w2 * xd;
    }

    @Command
    public final double getPotential(double x, double y, double z) {
        return interpolate(potential, x, y, z);
    }

    private final double getEnergy(double e, double m,
            double x, double y, double z,
            double vx, double vy, double vz) {
        return e * getPotential(x, y, z) + m * (vx*vx + vy*vy + vz*vz) / (2 * K);
    }

    /**
     * On the unit system
     *
     * The user specifies quantities in following units:
     *  e  |e| units
     *  E  eV
     *  m  amu
     *  V  V
     *  x  arb.u.
     *  t  arb.u.
     * Where arbitrary units for x and t are in agreement:
     * if x is in mm, then t in ms, x in um, then t in us, etc.
     *
     * The Newton's equation is
     *       e/m dV/dx = d2x/dt2, quantities in SI.
     * Transforming the eq. to our unit system:
     *       K e/m dV/dx = d2x/dt2, where K = |elemCharge|/1amu = |elemCharge|*NAvogadro*1000/1kg
     * The same K appears in v(E): v = sqrt(2E/m) SI => v = sqrt(2E/m K),
     *   Ekinetic = mv2/2K for the same reason.
     */
    private static final double K = 9.648534147e6;

    public void solve(double dt, double m, double e,
                double x0, double y0, double z0,
                double kineticE, double dirX, double dirY, double dirZ,
                int maxSteps,
                Callback... callbacks) {

        for (Callback c : callbacks) {
            c.params(dt, m, e, x0, y0, z0, kineticE, dirX, dirY, dirZ);
        }

        Vector dir = new Vector(dirX, dirY, dirZ);
        dir.normalizeIt();
        Vector v0 = dir.multiply(Math.sqrt(kineticE * 2 / m * K));

        int steps = 0;
        double t = 0;
        double x = x0, y = y0, z = z0;
        double vx = v0.getX(), vy = v0.getY(), vz = v0.getZ();

        double[] es;
        while ((maxSteps <= 0 || steps < maxSteps) &&
                x >= 0 && x <= dimX &&
                y >= 0 && y <= dimY &&
                z >= 0 && z <= dimZ) {
            steps++;

            es = getIntensity(x, y, z);
            if (es == null) {
                break;
            }
            t+= dt;

            vx += dt * es[X] * e / m * K;
            vy += dt * es[Y] * e / m * K;
            vz += dt * es[Z] * e / m * K;

            x += dt * vx;
            y += dt * vy;
            z += dt * vz;

            for (Callback c : callbacks) {
                c.point(t, x, y, z, vx, vy, vz);
            }
        }
        for (Callback c : callbacks) {
            c.close();
        }
    }

    private JFrame chartFrame = null;

    @Command
    public void showGraph() {
        if (chartFrame == null) {
            JFrame frame = new JFrame("Graph");
            frame.setSize(new Dimension(600, 400));
            ZoomableChart chart = new ZoomableChart();
            chart.addTrace(energyTrace);
            chart.addTrace(xCoordTrace);
            chart.addTrace(yCoordTrace);
            chart.addTrace(zCoordTrace);
            frame.getContentPane().add(new ChartPanel(chart));
            chartFrame = frame;
        }
        chartFrame.setVisible(true);
    }

    private static final int TRACE_SIZE_LIMIT = 10000;
    
    private ITrace2D energyTrace = new Trace2DLtd(TRACE_SIZE_LIMIT, "energy");
    private Callback energyPlottingCallback = new EnergyPlottingCallback();

    private ITrace2D xCoordTrace = new Trace2DLtd(TRACE_SIZE_LIMIT, "x");
    private ITrace2D yCoordTrace = new Trace2DLtd(TRACE_SIZE_LIMIT, "y");
    private ITrace2D zCoordTrace = new Trace2DLtd(TRACE_SIZE_LIMIT, "z");
    private Callback coordPlottingCallback = new XYZPlottingCallback();

    private abstract class PlottingCallback implements Callback {

        protected double m;
        protected double e;

        protected abstract void clearTraces();

        public void params(double dt, double m, double e, double x0, double y0, double z0, double kineticE, double dirX, double dirY, double dirZ) {
            clearTraces();
            this.m = m;
            this.e = e;
        }

        public abstract void point(double t, double x, double y, double z, double vx, double vy, double vz);

        public void close() {}
    }
    
    private class EnergyPlottingCallback extends PlottingCallback {

        protected void clearTraces() {
            energyTrace.removeAllPoints();
        }

        public void point(double t, double x, double y, double z, double vx, double vy, double vz) {
            energyTrace.addPoint(t, getEnergy(e, m, x, y, z, vx, vy, vz));
        }        
    }
    
    private class XYZPlottingCallback extends PlottingCallback {

        protected void clearTraces() {
            xCoordTrace.removeAllPoints();
            yCoordTrace.removeAllPoints();
            zCoordTrace.removeAllPoints();
        }

        public void point(double t, double x, double y, double z, double vx, double vy, double vz) {
            xCoordTrace.addPoint(t, x);
            yCoordTrace.addPoint(t, y);
            zCoordTrace.addPoint(t, z);
        }
    }

    private class OutputCallback implements Callback {

        private PrintStream s;

        public OutputCallback(PrintStream s) {
            this.s = s;
        }

        public void params(double dt, double m, double e, double x0, double y0, double z0, double kineticE, double dirX, double dirY, double dirZ) {
            s.println(String.format("# dt=%f, m=%f, e=%f, r0=(%f %f %f), K=%f, dirV0=(%f %f %f)",
                    dt, m, e, x0, y0, z0, kineticE, dirX, dirY, dirZ));
            s.println("t\tx\ty\tz\tvx\tvy\tvz");
        }

        public void point(double t, double x, double y, double z, double vx, double vy, double vz) {
            s.print(t); s.print('\t');

            s.print(x); s.print('\t');
            s.print(y); s.print('\t');
            s.print(z); s.print('\t');

            s.print(vx); s.print('\t');
            s.print(vy); s.print('\t');
            s.print(vz);

            s.println();
        }

        public void close() {
            s.close();
        }
    }

    private int maxSteps = 10000;
    private double dt = 0.001;
    private Vector r0;
    private Vector dirV0;
    private double k0;
    private double charge;
    private double mass;

    @Command
    public void setDt(double dt) {
        this.dt = dt;
    }

    @Command
    public void setMaxSteps(int maxSteps) {
        this.maxSteps = maxSteps;
    }

    @Command
    public void setParticle(double m, double e) {
        this.mass = m;
        this.charge = e;
    }

    @Command
    public void setV0(double energy, double dirX, double dirY, double dirZ) {
        this.k0 = energy;
        this.dirV0 = new Vector(dirX, dirY, dirZ);
    }

    @Command
    public void setR0(double x, double y, double z) {
        this.r0 = new Vector(x, y, z);
    }

    @Command
    public void fly(String outFileBase) throws FileNotFoundException {
        final Callback writer = new OutputCallback(new PrintStream(
                new BufferedOutputStream(new FileOutputStream(outFileBase + ".txt"))));
        solve(dt, mass, charge, 
                r0.getX(), r0.getY(), r0.getZ(),
                k0, dirV0.getX(), dirV0.getY(), dirV0.getZ(),
                maxSteps,
                energyPlottingCallback,
//                coordPlottingCallback,
                writer);
    }

    @Command
    public void exit() {
        System.exit(0);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        ShellFactory.createConsoleShell("solver", "Ion Trajectory Solver", new TrajectorySolver())
                .commandLoop();
    }

}
