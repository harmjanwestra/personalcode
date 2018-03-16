/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class SNPFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            TextFile in = new TextFile("", TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab);
            while (elems != null) {
                
                elems = in.readLineElems(TextFile.tab);
            }
            
            in.close();
            
            TextFile in2 = new TextFile("", TextFile.R);
            String[] elems2 = in2.readLineElems(TextFile.tab);
            while (elems != null) {
                
                elems2 = in2.readLineElems(TextFile.tab);
            }
            
            in.close();
        } catch (IOException e){
            
        }
    }
}
