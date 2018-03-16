/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harmjan
 */
public class ZScoreToPValueCheck {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        Descriptives.zScoreToPValue();
        for (int i = 0; i < Descriptives.m_zScoreToPValue.length; i++) {
            System.out.println(i + "\t" + Descriptives.m_zScoreToPValue[i]);
        }
    }
}
