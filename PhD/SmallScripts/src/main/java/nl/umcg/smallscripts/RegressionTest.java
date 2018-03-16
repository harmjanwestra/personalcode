/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.math.stats.Regression;

/**
 *
 * @author harmjan
 */
public class RegressionTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        double[] x = new double[]{0.2, 5.6, 9.8, 2.5, 4.5};
        double[] y = new double[]{5.7, 9.2, 7.8, 5.6, 1.3};

        double[] coeff1 = new double[]{0,1,2};
        double[] coeff2 = Regression.getLinearRegressionCoefficients(x, y);

        System.out.println("a\t"+coeff1[0]+"\t"+coeff2[0]);
        System.out.println("b\t"+coeff1[1]+"\t"+coeff2[1]);
        System.out.println("se\t"+coeff1[2]+"\t"+coeff2[2]);
//        System.out.println("se\t"+coeff1[2]+"\t"+coeff2[2]);
    }
}
