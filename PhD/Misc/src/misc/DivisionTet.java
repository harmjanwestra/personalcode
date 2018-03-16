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
public class DivisionTet {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        int nr = 1024 * 3 + 512 + 120;
        
        System.out.println(Gpio.humanizeFileSize(nr));
    }
}
