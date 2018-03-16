/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author harmjan
 */
public class InverseTTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        // 9.81E-198	0	
        System.out.println(ZScores.zScoreToCorrelation(-42.01369283875257, 2893));
    }
}
