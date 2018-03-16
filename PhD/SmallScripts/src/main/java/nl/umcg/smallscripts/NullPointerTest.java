/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.util.Primitives;

/**
 *
 * @author harmjan
 */
public class NullPointerTest {

    public static void main(String[] args) {
//        double x = -5.29;
//        double unit = 0.1;
//        double z = 2.3;
//        System.out.println((x % unit) + "\t" + (unit - (x % unit)) + "\t" + (x - unit - (x % unit)));
//
//        double[] x1 = new double[1000];
//        double[] y1 = new double[1000];
//        for (int i = -500; i < 500; i++) {
//            double p1 = Math.random();
//            double p2 = Math.random();
//
//            double p3 = Math.random();
//            double p4 = Math.random();
//
//            if (p1 > 0.5) {
//                p3 *= -1;
//            }
//
//            if (p2 > 0.5) {
//                p4 *= -1;
//            }
//            
//            double q = (2*Math.PI)/1000*i;
//            
//            x1[i+500] = Math.sin(q);
//            y1[i+500] = Math.cos(q);
////            y1[i] -= 300d;
//        }
//
//
//        ScatterPlot p = new ScatterPlot(500, 500, 0.1, x1, y1, "/Volumes/iSnackHD/tmp/fig.png");

        double x = 0.1;
        double z = 0.01;
        double f = 101;
        System.out.println(Math.log10(x));
        System.out.println(Math.log10(z));
        System.out.println(Math.log10(f));
        
        System.out.println(Math.pow(10, -0.1));

        //        Double[][] test = new Double[100][48000];
//
//        for(int i=0; i<test.length; i++){
//            for(int j=0; j<test[i].length; j++){
//                test[i][j] = Math.random();
//            }
//        }
//        try {
//            System.out.println("Waiting");
//            Thread.sleep(30000);
//        } catch (InterruptedException e) {
//            System.out.println("got interrupted!");
//        }
    }
}
