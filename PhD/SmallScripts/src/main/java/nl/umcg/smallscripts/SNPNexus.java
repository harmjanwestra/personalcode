/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class SNPNexus {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try{
            TextFile in = new TextFile("/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/Real/AllSNPs.txt-WithProxies.txt", TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab);
            HashSet<String> snps = new HashSet<String>();
            while (elems != null) {
                snps.add(elems[0]);
                snps.add(elems[1]);
                elems = in.readLineElems(TextFile.tab);
            }
            in.close();
            
            TextFile out = new TextFile("/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/Real/AllSNPs.txt-WithProxies-SNPNexus.txt", TextFile.W);
            TextFile out2 = new TextFile("/Volumes/iSnackHD/SkyDrive/SNPFunctionalAnnotation/Trans-SNPs/Real/AllSNPs.txt-WithProxies-List.txt", TextFile.W);
            for(String s: snps){
                out.writeln("dbsnp\t"+s);
                out2.writeln(s);
            }
            out.close();
            out2.close();
            
        } catch (IOException e){
            e.printStackTrace();
        }
    }
}
