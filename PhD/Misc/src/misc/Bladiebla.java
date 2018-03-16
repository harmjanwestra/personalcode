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
public class Bladiebla {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        for(double i=0; i<30; i++){
            System.out.println(i+"\t"+ZScores.zToP(i));
        }
    }
}
