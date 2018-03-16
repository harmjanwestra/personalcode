/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package dashavcfconverter;

import java.io.IOException;
import umcg.genetica.io.trityper.converters.VCFToTriTyper;

/**
 *
 * @author harmjan
 */
public class DashaVCFConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            VCFToTriTyper t = new VCFToTriTyper();
            if (args.length == 2) {
                t.parse(args[0], args[1]);
            } else if (args.length == 3) {
                t.parse(args[0], args[1], args[2]);
            } else {
                System.out.println("Usage: DashaVCFConverter.jar indir outdir [snppattern]");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
