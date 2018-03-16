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
public class CreateRScripts {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
//            for (int chr = 1; chr < 23; chr++) {
//                String fileIn  = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.Rdata";
//                String fileOut = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.csv";
//                TextFile out   = new TextFile("/Volumes/iSnackHD/ImpData/rscripts/r" + chr + "-1240.R", TextFile.W);
//
//                String script  = "load(\"" + fileIn + "\");\nwrite.csv(file=\"" + fileOut + "\", as(genotypes, 'numeric'));\nq()\n";
//                out.writeln(script);
//                out.close();
//                
//                // set_229_chr11_imputed
//                fileIn  = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.Rdata";
//                fileOut = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.csv";
//                out   = new TextFile("/Volumes/iSnackHD/ImpData/rscripts/r" + chr + "-229.R", TextFile.W);
//
//                script  = "load(\"" + fileIn + "\");\nwrite.csv(file=\"" + fileOut + "\", as(genotypes, 'numeric'));\nq()\n";
//                out.writeln(script);
//                out.close();
//            }
            
            for (int chr = 1; chr < 23; chr++) {
                String fileIn  = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.Rdata";
                String fileOut = "/Volumes/iSnackHD/ImpData/set_1240_chr" + chr + "_imputed.support.txt";
                TextFile out   = new TextFile("/Volumes/iSnackHD/ImpData/rscripts/r" + chr + "-1240.R", TextFile.W);

                String script  = "library(snpStats);\n"
                        + "load(\"" + fileIn + "\");\n"
                        + "write.table(SNP.support, file=\""+fileOut+"\", quote=FALSE, sep=\"\t\");\n"
                        + "q()\n";
                out.writeln(script);
                out.close();
                
                // set_229_chr11_imputed
                fileIn  = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.Rdata";
                fileOut = "/Volumes/iSnackHD/ImpData/set_229_chr" + chr + "_imputed.support.txt";
                out   = new TextFile("/Volumes/iSnackHD/ImpData/rscripts/r" + chr + "-229.R", TextFile.W);

                script  = "library(snpStats);\n"
                        + "load(\"" + fileIn + "\");\n"
                        + "write.table(SNP.support, file=\""+fileOut+"\", quote=FALSE, sep=\"\t\");\n"
                        + "q()\n";
                
                out.writeln(script);
                out.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
