/*
Snyder Projection Implementation

Copyright (c) 2012-2015, APIS Point, LLC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
package gis.proj.drivers;

import gis.proj.Azimuthal;
import gis.proj.Datum;
import gis.proj.Ellipsoid;
import gis.proj.Projection;
import gis.proj.SnyderMath;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Path2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.Semaphore;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.ListCellRenderer;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.event.MouseInputListener;

import com.beust.jcommander.IParameterValidator;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

public final class SnyderVis {

    private static class MapCanvas extends JComponent implements MouseInputListener, MouseWheelListener {
        private static final long serialVersionUID = 1L;

        private final static Color INDIAN_RED = new Color(205,  92,  92, 255);
        private final static Color ROYAL_BLUE = new Color( 65, 105, 225, 255);
        private final static Color PERU_BROWN = new Color(205, 133,  63, 255);

        private final static Color TRNS_BCKG = new Color(  0,   0,   0, 192);
        private final static Color TRNS_FONT = new Color(255, 255, 255, 128);

        private final static Font   HUD_FONT = new Font(Font.MONOSPACED, Font.PLAIN, 25);
        private final static String HUD_STRG = "pps\u2245%,15.3f";

        private double scale = 100;
        private double dt;
        private double projDivisor;

        private boolean contrastMode = false;

        //
        // Animation variables
        //
        private Thread rotateThread = null;

        private boolean rotLon = false;
        private boolean rotLat = false;

        private double rotation    = 0;
        private double rotationInc = StrictMath.toRadians(1.0);

        private ArrayList<double[][]> renderLines = new ArrayList<double[][]>();
        private double[] env;
        private double[] envStart = {Double.MAX_VALUE, Double.MAX_VALUE, Double.MIN_VALUE, Double.MIN_VALUE};

        // The animation mutex dynamically throttles the drawing
        private Semaphore animationMutex = new Semaphore(1, true);

        //
        // Generalized drawing variables
        //
        private boolean invertedY = true;
        private double  xoff      = 0;
        private double  yoff      = 0;
        private double  xoffs     = 0;
        private double  yoffs     = 0;

        private Point sLoc;

        private Path2D path   = new Path2D.Float();
        private Stroke stroke = new BasicStroke(0.0f); // Force a 1px stroke regardless of scale

        private boolean drawAsInv = false;

        private Graphics2D g2;
        private double[][] projPnt;

        //
        // Line data
        //
        private ArrayList<double[][]> worldLines        = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleLonLines = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleLatLines = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleEquator  = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleMeridian = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleDateLine = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleNPoleCir = new ArrayList<double[][]>();
        private ArrayList<double[][]> graticuleSPoleCir = new ArrayList<double[][]>();

        //
        // Graticule variables
        //
        private boolean drawLonGraticule = false;
        private boolean drawLatGraticule = false;

        //
        // Projection variables
        //
        private Projection proj      = null;
        private Datum      datum     = null;
        private Ellipsoid  ellipsoid = null;

        private boolean inverseProjection = false;

        private String lonRotStr = "lon0";
        private String latRotStr = "lat1";
        private double lonRot    = 0.0;
        private double latRot    = 0.0;

        //
        // Points Per Second (PPS)
        //
        private double pps   = 0;
        private long   ppsst = 0; // Start Time;
        private double ppset = 0; // Elapse Time

        public MapCanvas(
                String dataSetFileName,
                double dt) {

            setMinimumSize(new Dimension(600, 400));
            setPreferredSize(getMinimumSize());

            this.dt = dt > 1.0 ? 1.0 : dt;

            loadDataSet(dataSetFileName);
            setEllipsoid("WGS84");

            addMouseListener(this);
            addMouseMotionListener(this);
            addMouseWheelListener(this);

            //
            // Graticule construction
            //
            int    step;
            int    ticks;
            int    hlf_ticks;
            int    degStep = 15;
            double degPace = 00.1;

            double fillVal;
            double[] vals;

            //
            // Latitudes
            //
            ticks     = (int) StrictMath.ceil(360.0 / degPace) + 1;
            hlf_ticks = ticks >> 1;

            double[][] equator = new double[ticks][ticks];
            for(step = 1; step <= hlf_ticks; step++) {
                equator[0][hlf_ticks + step] = step * degPace * SnyderMath.DEG_TO_RAD;
                equator[0][hlf_ticks - step] = -equator[0][hlf_ticks + step];
            }
            graticuleEquator.add(equator);

            for(step = degStep; step < 90; step += degStep) {
                fillVal = step * SnyderMath.DEG_TO_RAD;

                vals = new double[ticks];
                Arrays.fill(vals, fillVal);
                graticuleLatLines.add(new double[][]{equator[0], vals});

                vals = new double[ticks];
                Arrays.fill(vals, -fillVal);
                graticuleLatLines.add(new double[][]{equator[0], vals});
            }

            fillVal = (66.0 + 33.0 / 60.0 + 44.0 / 3600.0) * SnyderMath.DEG_TO_RAD;

            vals = new double[ticks];
            Arrays.fill(vals,  fillVal);
            graticuleNPoleCir.add(new double[][]{equator[0], vals});

            vals = new double[ticks];
            Arrays.fill(vals, -fillVal);
            graticuleSPoleCir.add(new double[][]{equator[0], vals});

            //
            // Longitudes
            //
            ticks     = (int) StrictMath.ceil(180 / degPace) + 1;
            hlf_ticks = ticks >> 1;

            double[][] meridian = new double[ticks][ticks];
            for(step = 1; step <= hlf_ticks; step++) {
                meridian[1][hlf_ticks + step] = step * degPace * SnyderMath.DEG_TO_RAD;
                meridian[1][hlf_ticks - step] = -meridian[1][hlf_ticks + step];
            }
            graticuleMeridian.add(meridian);

            for(step = degStep; step < 180; step += degStep) {
                fillVal = step * SnyderMath.DEG_TO_RAD;

                vals = new double[ticks];
                Arrays.fill(vals, fillVal);
                graticuleLonLines.add(new double[][] {vals, meridian[1]});

                vals = new double[ticks];
                Arrays.fill(vals, -fillVal);
                graticuleLonLines.add(new double[][] {vals, meridian[1]});
            }

            vals = new double[ticks];
            Arrays.fill(vals, SnyderMath.PI);
            graticuleDateLine.add(new double[][] {vals, meridian[1]});
        }

        private double[] envelope(double[][] input, double env[]) {
            double minx = env[0],
                    miny = env[1],
                    maxx = env[2],
                    maxy = env[3];

            for(int i = 0; i < input[0].length; ++i) {
                if(input[0][i] < minx) minx = input[0][i];
                if(input[0][i] > maxx) maxx = input[0][i];
                if(input[1][i] < miny) miny = input[1][i];
                if(input[1][i] > maxy) maxy = input[1][i];
            }

            return new double[] {
                minx, miny, maxx, maxy
            };
        }

        private void drawPoly(Graphics2D g, double[][] pnts, double[] extent) {
            path.reset();
            path.moveTo(pnts[0][0], pnts[1][0]);

            for(int i = 1, j = 0; i < pnts[0].length; ++i, ++j) {
                if(
                    StrictMath.abs(pnts[0][i] - pnts[0][j]) / extent[0] < dt &&
                    StrictMath.abs(pnts[1][i] - pnts[1][j]) / extent[1] < dt)
                    path.lineTo(pnts[0][i], pnts[1][i]);
                else
                    path.moveTo(pnts[0][i], pnts[1][i]);
            }

            g.draw(path);
        }

        private void projectDraw(ArrayList<double[][]> data, Graphics2D g) {
            // Zeroize variables used for actual BogoPPS calculations
            pps   = 0.0;
            ppset = 0.0;

            renderLines.clear();
            env = envStart;

            for(double[][] points : data) {
                pps += points[0].length;
                ppsst = System.nanoTime();

                projPnt = proj.forward(points[0], points[1], ellipsoid, datum);
                ppset += System.nanoTime() - ppsst;

                if(inverseProjection) {
                    ppsst = System.nanoTime();

                    projPnt = proj.inverse(projPnt[0], projPnt[1],
                            ellipsoid,
                            datum);

                    ppset += System.nanoTime() - ppsst;
                }

                env = envelope(projPnt, env);
                renderLines.add(projPnt);
            }

            // Do not including the drawing time
            for(double[][] p : renderLines)
                try {
                    drawPoly(g, p, new double[] {
                            StrictMath.abs(env[2] - env[0]),
                            StrictMath.abs(env[3] - env[1])});
                } catch(Exception e) {
                    System.err.println("Error drawing to screen!!");
                }
        }

        public void paintComponent(Graphics grphx) {
            // Defaults for drawing without a selected projections
            pps   = 0.0;
            ppset = 1.0;

            super.paintComponent(grphx);

            if(contrastMode == true) {
                grphx.setColor(Color.DARK_GRAY);
                grphx.fillRect(0, 0, getWidth(), getHeight() - 1);
                grphx.setColor(Color.GRAY);
            }

            if(proj != null && datum != null && ellipsoid != null) {
                g2 = (Graphics2D) grphx.create();

                drawAsInv =
                        Azimuthal.class.isAssignableFrom(proj.getClass()) &
                        inverseProjection == false &
                        rotation != 0.0;

                g2.setClip(0, 0, getWidth(), getHeight());
                g2.setStroke(stroke);

                if(drawAsInv)
                    g2.rotate(rotation, getWidth() * 0.5, getHeight() * 0.5 );

                g2.translate(
                        getWidth()  * 0.5 + (drawAsInv ? -1 : 1) * xoff,
                        getHeight() * 0.5 + (drawAsInv ? -1 : 1) * (invertedY ? yoff : -yoff));

                projDivisor = inverseProjection
                    ? 1.0 : ellipsoid.isSphere() ? ellipsoid.getProperty("R") : ellipsoid.getProperty("a");

                g2.scale(scale / projDivisor, invertedY ? -scale / projDivisor : scale / projDivisor);

                // Draw the world data
                projectDraw(worldLines, g2);

                if(drawLonGraticule)
                    projectDraw(graticuleLonLines, g2);

                if(drawLatGraticule)
                    projectDraw(graticuleLatLines, g2);

                g2.setColor(ROYAL_BLUE);
                projectDraw(graticuleEquator, g2);

                g2.setColor(PERU_BROWN);
                projectDraw(graticuleNPoleCir, g2);
                projectDraw(graticuleSPoleCir, g2);

                g2.setColor(INDIAN_RED);
                projectDraw(graticuleMeridian, g2);
                projectDraw(graticuleDateLine, g2);
            }

            ((Graphics2D) grphx).setRenderingHint(
                    RenderingHints.KEY_TEXT_ANTIALIASING,
                    RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

            grphx.setColor(TRNS_BCKG);
            grphx.fillRect(0, 0, getWidth(), 38);
            grphx.setColor(TRNS_FONT);
            grphx.setFont(HUD_FONT);

            grphx.drawString(String.format(HUD_STRG, pps * (1.0E9 / ppset)), 10, 25);

            animationMutex.release();
        }

        public void mousePressed(MouseEvent e) {
            sLoc = e.getPoint();
            xoffs = xoff;
            yoffs = yoff;
        }

        public void mouseClicked(MouseEvent e) { }
        public void mouseReleased(MouseEvent e) { }
        public void mouseEntered(MouseEvent e) { }
        public void mouseExited(MouseEvent e) { }

        public void mouseDragged(MouseEvent e) {
            xoff = xoffs + e.getPoint().x - sLoc.x;
            yoff = yoffs + e.getPoint().y - sLoc.y;
            repaint();
        }

        public void mouseMoved(MouseEvent e) { }

        public void mouseWheelMoved(MouseWheelEvent e) {
            scale += e.getWheelRotation() < 0 ? -1.0 : 1.0;
            if(scale < 1.0)
                scale = 1.0;
            repaint();
        }

        private void rotate() {
            if(rotLon) {
                lonRot -= rotationInc;
                datum.setUserOverrideProperty(lonRotStr, SnyderMath.normalizeLonRad(lonRot));
            }
            if(rotLat) {
                latRot += rotation > SnyderMath.PI_DIV_2 ? -rotationInc : rotationInc;
                if(StrictMath.abs(latRot) >= SnyderMath.POLE_RAD) {
                    // If the animation starts on either pole correct for a stuck earth
                    // caused by incrementing from 90 to 91
                    latRot = rotation > SnyderMath.PI_DIV_2 ? SnyderMath.N_PI_DIV_2 + rotationInc : SnyderMath.PI_DIV_2 - rotationInc;

                    rotation = rotation > SnyderMath.PI_DIV_2 ? 0 : SnyderMath.PI;
                    if(Azimuthal.class.isAssignableFrom(proj.getClass())) {
                        lonRot += SnyderMath.PI;
                        datum.setUserOverrideProperty(lonRotStr, SnyderMath.normalizeLonRad(lonRot));
                    }
                }
                datum.setUserOverrideProperty(latRotStr, latRot);
            }
            repaint();
        }

        public Projection getProjection() {
            return this.proj;
        }

        public void setProjection(Projection proj) {
            this.proj = proj;
            this.datum = getDatum();

            repaint();
            startRotation();
        }

        public void setLongitudeRotation(String str) {
            this.lonRotStr = str;
        }

        public void setLatitudeRotation(String str) {
            this.latRotStr = str;
        }

        public String getLongitudeRotation() {
            return this.lonRotStr;
        }

        public String getLatitudeRotation() {
            return this.latRotStr;
        }

        public void drawGraticule(boolean lat, boolean b) {
            if(lat)
                drawLatGraticule = b;
            else
                drawLonGraticule = b;
            repaint();
        }

        private void startRotation() {
            if(proj != null && (rotateThread == null || rotateThread.isAlive() == false) && (rotLat || rotLon)) {
                rotateThread = new Thread(new Runnable() {
                    public void run() {
                        do {
                            animationMutex.acquireUninterruptibly();
                            if(datum != null)
                                rotate();
                        } while(rotLat || rotLon);
                    }
                });
                rotateThread.setDaemon(true);
                rotateThread.setName("Rotation");
                rotateThread.start();
            }
        }

        public void rotate(boolean lat, boolean b) {
            if(lat)
                rotLat = b;
            else
                rotLon = b;

            startRotation();
        }

        public void doInverseProjection(boolean b) {
            inverseProjection = b;
            repaint();
        }

        public void resetMapView() {
            xoff  = 0;
            yoff  = 0;
            xoffs = 0;
            yoffs = 0;
            scale = 100;
            repaint();
        }

        public String getEllipsoid() {
            return ellipsoid.getName();
        }

        public void setEllipsoid(String name) {
            if(EllipsoidFactory.getInstance().getEllipsoidNames().contains(name)) {
                ellipsoid = EllipsoidFactory.getInstance().getEllipsoid(name);
                repaint();
            }
        }

        public Datum getDatum() {
            return proj != null ? DatumFactory.getInstance().getDatum(proj.getClass().getName()) : null;
        }

        public void setContrastMode(boolean b) {
            this.contrastMode = b;
            repaint();
        }

        public void loadDataSet(String filename) {
            Scanner s = null;

            try {
                s = new Scanner(new File(filename));
                double[] lon, lat;
                int pCount;

                worldLines.clear();
                System.gc();

                while(s.hasNext()) {
                    s.next();
                    s.next();
                    pCount = s.nextInt();
                    s.nextLine();

                    lon = new double[pCount];
                    lat = new double[pCount];

                    for(int i = 0; i < pCount; i++) {
                        lon[i] = s.nextDouble() * SnyderMath.DEG_TO_RAD;
                        lat[i] = s.nextDouble() * SnyderMath.DEG_TO_RAD;
                    }
                    worldLines.add(new double[][] {lon, lat});
                }
            } catch(Exception e) {
                e.printStackTrace();
            } finally {
                try {
                    if(s != null)
                        s.close();
                } catch(Exception e) {}
            }

            repaint();
        }

    } // class

    private static final StringBuilder sb = new StringBuilder();

    private static final String WRN_TITLE = "Uh-Oh!";
    private static final String INF_TITLE = "Interesting...";

    @Parameter(names = "-data", description = "Data set filename")
    private String dataSet = Snyder.WMAP_DIR + "gshhg_l.txt";

    @Parameter(names = "-dt", description = "Positive drawing threshold, percent [0, 1.0]", validateWith = PositiveDouble.class)
    private double dt = 0.05;

    private SnyderVis() {}

    public static class PositiveDouble implements IParameterValidator {
        public void validate(String name, String value) throws ParameterException {
            double n = Double.parseDouble(value);
            if (n <= 0.0) {
                throw new ParameterException("Parameter " + name + " should be positive (found " + value +")");
            }
        }
    }

    private static String getEllipsoidInformation(String ellipsoid) {
        String[] keys = EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getPropertyNames().toArray(
                new String[EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getPropertyNames().size()]);
        Arrays.sort(keys);

        // This is bad, but it's a way quick to get the first
        // three elements in this order.
        keys[0] = "a";
        keys[1] = "1/f";
        keys[2] = "R";

        sb.delete(0, sb.length());
        sb.append("<html><b>Ellipsoid Parameters:</b><br><br><table border='0' cellpadding='0' cellspacing='0'>");
        for(String p : keys) {
            sb.append("<tr><td style='width:10px'/><td style='text-align:right'>");
            sb.append(p);
            sb.append("</td><td style='width:10px'/><td>");
            sb.append(EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getProperty(p));
            sb.append("</td></tr>");
        }
        sb.append("</table></html>");

        return sb.toString();
    }

    private static String getEllipsoidNameInformation(String ellipsoid) {
        sb.delete(0, sb.length());

        sb.append("<html>");
        sb.append(EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getName());
        sb.append(" [ID: ");
        sb.append(EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getId());
        sb.append("]<br>");
        sb.append(EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).getDescription());
        sb.append("<br>");
        sb.append(EllipsoidFactory.getInstance().getEllipsoid(ellipsoid).isSphere() ? "Spherical" : "Ellipsoidal");
        sb.append("</html>");

        return sb.toString();
    }

    private void realizeGuiAndDisplay() {
        final JFrame f = new JFrame("// Snyder Vis");
        final MapCanvas mapCanvas = new MapCanvas(dataSet, dt);

        //
        // PROJECTION PROPERTIES
        //
        ArrayList<Projection> projClasses = new ArrayList<Projection>();
        Set<String>           posPClasses = DatumFactory.getInstance().getDatumNames();

        for(String prj : posPClasses) {
            try {
                Class<?> cls = Class.forName(prj);
                if(Projection.class.isAssignableFrom(cls) == true)
                    projClasses.add((Projection) cls.newInstance());
            } catch (Exception e) {
                System.err.println("Unable to load projection: " + prj);
            }
        }

        //
        // Sort the projections by name
        //
        Comparator<Projection> comparator = new Comparator<Projection>() {
            public int compare(Projection c1, Projection c2) {
                return c1.getName().compareTo(c2.getName());
            }
        };

        Collections.sort(projClasses, comparator);

        // Add none to top of sorted list
        projClasses.add(0, null);

        JLabel projectionsLabel = new JLabel("Projections");

        JComboBox<Projection> projComboBox = new JComboBox<Projection>(
                projClasses.toArray(new Projection[projClasses.size()]));

        projComboBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Projection p = (Projection) ((JComboBox<?>) e.getSource()).getSelectedItem();
                mapCanvas.setProjection(p);
                if(p != null && DatumFactory.getInstance().getDatum(p.getClass().getName()) == null)
                    JOptionPane.showMessageDialog(
                            f,
                            "Projection does not define a datum; it cannot render.",
                            WRN_TITLE,
                            JOptionPane.WARNING_MESSAGE);
            }
        });
        projComboBox.setRenderer(new ListCellRenderer<Projection>() {
            public Component getListCellRendererComponent(
                    JList<? extends Projection> list, Projection value, int index,
                    boolean isSelected, boolean cellHasFocus) {
                return new JLabel(value != null ? value.getName() : "None");
            }
        });

        MouseListener ProjExpansionML = new MouseListener() {
            public void mouseReleased(MouseEvent e) { }
            public void mousePressed(MouseEvent e)  { }
            public void mouseExited(MouseEvent e)   { }
            public void mouseEntered(MouseEvent e)  { }
            public void mouseClicked(MouseEvent e) {
                Projection proj = mapCanvas.getProjection();

                if(SwingUtilities.isRightMouseButton(e))
                    if(proj != null) {
                        Datum d = mapCanvas.getDatum();

                        JPanel paramsPanel = new JPanel(new GridLayout(0, 1));
                        JPanel paramPanel  = null;
                        JTextField jtf = null;

                        ArrayList<JTextField> jtfl = new ArrayList<JTextField>();

                        try {
                            for(String param : proj.getDatumProperties()) {
                                paramPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
                                paramPanel.add(new JLabel("<html><b>" + param + " </b></html>"));
                                jtf = new JTextField(
                                        Double.isNaN(d.getProperty(param)) ?
                                                "" :
                                                    Double.toString(d.getProperty(param)), 20);
                                jtf.setName(param);
                                jtfl.add(jtf);
                                paramPanel.add(jtf);
                                paramsPanel.add(paramPanel);
                            }
                        } catch(Exception ex) { }

                        if(jtfl.isEmpty() == false) {
                            int result = JOptionPane.showConfirmDialog(
                                    f,
                                    paramsPanel,
                                    "Projection parameters:",
                                    JOptionPane.OK_CANCEL_OPTION);

                            if (result == JOptionPane.OK_OPTION) {
                                for(JTextField tf : jtfl){
                                    if(tf.getText().isEmpty() == false)
                                        d.setUserOverrideProperty(
                                                tf.getName(),
                                                Snyder.parseDatumVal(tf.getText()));
                                }
                                mapCanvas.repaint();
                            }
                        } else {
                            JOptionPane.showMessageDialog(
                                    f,
                                    "This projection is not reporting\ndatum parameters.",
                                    INF_TITLE,
                                    JOptionPane.INFORMATION_MESSAGE);
                        }
                    }
                    else
                        JOptionPane.showMessageDialog(f, "Please select a projection.", WRN_TITLE, JOptionPane.WARNING_MESSAGE);
            }
        };
        JLabel projExpansion = new JLabel(">");
        projExpansion.addMouseListener(ProjExpansionML);

        JCheckBox inverseProj = new JCheckBox("Inv");
        inverseProj.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.doInverseProjection(e.getStateChange() == ItemEvent.SELECTED);
            }
        });

        JPanel projPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        projPanel.add(projectionsLabel);
        projPanel.add(projComboBox);
        projPanel.add(projExpansion);
        projPanel.add(inverseProj);

        //
        // ELLIPSOID PROPERTIES
        //
        JLabel ellipsoidsLabel = new JLabel("Ellipsoids");
        final JLabel ellipsoidNameLabel = new JLabel(getEllipsoidNameInformation(mapCanvas.getEllipsoid()));
        final JLabel ellipsoidsInfo = new JLabel(getEllipsoidInformation(mapCanvas.getEllipsoid()));

        // Sort ellipsoids
        String[] ellipsoidsList = EllipsoidFactory.getInstance().getEllipsoidNames().toArray(new String[EllipsoidFactory.getInstance().getEllipsoidNames().size()]);
        Arrays.sort(ellipsoidsList);

        JComboBox<String> ellipsoidsComboBox = new JComboBox<String>(ellipsoidsList);
        ellipsoidsComboBox.setSelectedItem(mapCanvas.getEllipsoid());
        ellipsoidsComboBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String ellipsoid = (String) ((JComboBox<?>) e.getSource()).getSelectedItem();
                mapCanvas.setEllipsoid(ellipsoid);
                ellipsoidNameLabel.setText(getEllipsoidNameInformation(ellipsoid));
                ellipsoidsInfo.setText(getEllipsoidInformation(ellipsoid));
            }
        });

        JPanel ellDropPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        ellDropPanel.add(ellipsoidsLabel);
        ellDropPanel.add(ellipsoidsComboBox);
        ellDropPanel.add(ellipsoidNameLabel);

        //
        // DRAWING PROPERTIES
        //
        JPanel drawingPropPanel = new JPanel(new GridLayout(0, 2));

        JCheckBox graticuleLon = new JCheckBox("Graticule Longitude");
        JCheckBox graticuleLat = new JCheckBox("Graticule Latitude");
        graticuleLon.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.drawGraticule(false, e.getStateChange() == ItemEvent.SELECTED);
            }
        });
        graticuleLat.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.drawGraticule(true, e.getStateChange() == ItemEvent.SELECTED);
            }
        });
        drawingPropPanel.add(graticuleLon);
        drawingPropPanel.add(graticuleLat);

        MouseListener checkBoxML = new MouseListener() {
            public void mouseReleased(MouseEvent e) { }
            public void mousePressed(MouseEvent e)  { }
            public void mouseExited(MouseEvent e)   { }
            public void mouseEntered(MouseEvent e)  { }
            public void mouseClicked(MouseEvent e) {
                if(SwingUtilities.isRightMouseButton(e)) {
                    JCheckBox src = (JCheckBox) e.getSource();
                    boolean lat = src.getText().matches(".*Latitude.*");
                    String seed = lat ?
                            mapCanvas.getLatitudeRotation() :
                            mapCanvas.getLongitudeRotation();

                    String tmp = JOptionPane.showInputDialog(
                           f,
                           src.getText().substring(0, src.getText().length() - 1)
                           + "projection parameter:",
                           seed);

                    if(tmp != null && tmp.trim().isEmpty() == false) {
                        if(lat)
                            mapCanvas.setLatitudeRotation(tmp);
                        else
                            mapCanvas.setLongitudeRotation(tmp);
                    }
                }
            }
        };

        JCheckBox rotLon = new JCheckBox("Rotate Longitude >");
        rotLon.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.rotate(false, e.getStateChange() == ItemEvent.SELECTED);
            }
        });
        rotLon.addMouseListener(checkBoxML);
        drawingPropPanel.add(rotLon);

        JCheckBox rotLat = new JCheckBox("Rotate Latitude >");
        rotLat.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.rotate(true, e.getStateChange() == ItemEvent.SELECTED);
            }
        });
        rotLat.addMouseListener(checkBoxML);
        drawingPropPanel.add(rotLat);

        JCheckBox cmDrawing = new JCheckBox("Contrast mode");
        cmDrawing.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                mapCanvas.setContrastMode(e.getStateChange() == ItemEvent.SELECTED);
            }
        });
        drawingPropPanel.add(cmDrawing);

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));

        JButton viewResetButton = new JButton("Reset Map View");
        viewResetButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                mapCanvas.resetMapView();
            }
        });
        buttonPanel.add(viewResetButton);

        JButton wdbButton = new JButton("Select World Data");
        wdbButton.addActionListener(new ActionListener() {
            private JFileChooser fc = new JFileChooser(Snyder.WMAP_DIR);

            public void actionPerformed(ActionEvent e) {
                ((JButton) e.getSource()).setEnabled(false);

                if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
                {
                    File file = fc.getSelectedFile();
                    if (file.canRead() && file.exists())
                    {
                        mapCanvas.loadDataSet(file.getAbsolutePath());
                    }
                }

                ((JButton) e.getSource()).setEnabled(true);
            }
        });
        buttonPanel.add(wdbButton);

        JPanel panel = new JPanel();
        Border drawingBorder = BorderFactory.createTitledBorder(
                BorderFactory.createLineBorder(Color.BLACK),
                "Properties");
        BoxLayout bl = new BoxLayout(panel, BoxLayout.Y_AXIS);

        JPanel ellipsoidInfoPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        ellipsoidInfoPanel.add(ellipsoidsInfo);

        panel.setLayout(bl);
        panel.setBorder(drawingBorder);
        panel.add(projPanel);
        panel.add(ellDropPanel);
        panel.add(ellipsoidInfoPanel);
        panel.add(drawingPropPanel);
        panel.add(buttonPanel);

        JLabel info = new JLabel(Snyder.licenseInformation("SnyderVis"));

        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(mapCanvas, BorderLayout.CENTER);
        f.getContentPane().add(panel, BorderLayout.WEST);
        f.getContentPane().add(info, BorderLayout.SOUTH);
        f.pack();
        f.setLocationRelativeTo(null);
        f.setVisible(true);
    }

    public static void main(String... args) {
        final SnyderVis snyderVis = new SnyderVis();
        JCommander jc = new JCommander(snyderVis);

        try {
            jc.parse(args);
        } catch(Exception e) {
            jc.usage();
            System.exit(-10);
        }

        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                snyderVis.realizeGuiAndDisplay();
            }
        });
    }

}