/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.visualisation;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import umcg.genetica.containers.Pair;
import umcg.genetica.math.Goniometry;

/**
 *
 * @author harmjan
 */
public class CircularPlot {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            CircularPlot p = new CircularPlot();
            p.run("/Volumes/iSnackHD/circ.png");
        } catch (IOException e) {
        }
    }
    private BufferedImage bimage;
    private Graphics2D g2d;

    public void run(String loc) throws IOException {
        int width = 500;
        int height = 500;

        int marginx = 50;
        int rx = 100;
        int ry = 100;

        int dx = rx * 2;
        int dy = ry * 2;

        bimage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        g2d = bimage.createGraphics();
//        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, width, height);

        g2d.setColor(Color.black);


        int originX = width / 2 - rx;
        int originY = height / 2 - ry;
        g2d.fillOval(originX, originY, dx, dy);

        
        int startX = originX + rx + 5;
        int startY = originY + ry + 5;
        
        Pair<Integer, Integer> startX1 = Goniometry.calcPosOnCircle(rx + 5, width / 2 - rx, (height / 2) - (int) Math.floor(ry), 0); 
        Pair<Integer, Integer> startX2 = Goniometry.calcPosOnCircle(rx + 10, width / 2 - rx, (height / 2) - (int) Math.floor(ry), 0);
        
        g2d.dispose();
        javax.imageio.ImageIO.write(bimage, "png", new File(loc));
    }
}
