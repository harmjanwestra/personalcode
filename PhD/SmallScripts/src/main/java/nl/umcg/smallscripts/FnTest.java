/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class FnTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        
        System.out.println(Gpio.getFileName("/Volumes/iSnackHD/SoftwareTest/ExpressionData.txt"));
    }
}
