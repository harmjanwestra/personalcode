/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package pedandmaptotrityperconverter;

import java.io.IOException;
import umcg.genetica.io.trityper.converters.PedAndMapToTriTyper;

/**
 *
 * @author harmjan
 */
public class PedAndMapToTriTyperConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        if (args.length < 2) {
            System.out.println("Usage: indir outdir");
        } else {
            try {
                PedAndMapToTriTyper pd = new PedAndMapToTriTyper();
                pd.importPEDFile(args[0], args[1]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }
}
