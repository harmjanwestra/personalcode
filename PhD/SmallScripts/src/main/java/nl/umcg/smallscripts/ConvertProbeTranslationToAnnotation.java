/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ConvertProbeTranslationToAnnotation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String annot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";

        String outfilePrefix = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-";
        boolean ht12v4 = true;
        try {


            TextFile in = new TextFile(annot, TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab); // skip the header

            // Platform        HT12v4-ArrayAddress     Symbol  Chr     ChrStart        ChrEnd  Probe   Seq

            elems = in.readLineElems(TextFile.tab);
            String platform = "HT12v3";
            int col = 5;
            if (ht12v4) {
                platform = "HT12v4";
                col = 7;
            }

            TextFile out = new TextFile(outfilePrefix + platform + ".txt", TextFile.W);
            out.writeln("Platform\t" + platform + "-ArrayAddress\tSymbol\tChr\tChrStart\tChrEnd\tProbe\tSeq");
            while (elems != null) {

                String platformId = elems[col];
                if (!platformId.equals("-")) {
                    String chrpos = elems[3];
                    String chrstart = "";
                    String chrend = "";
                    String[] chrposElems = chrpos.split(":");
                    if (chrposElems.length > 1) {
                        // exon boundary crossing probe
                        ArrayList<String> pos = new ArrayList<String>();
                        String firstElem = chrposElems[0];
                        String lastElem = chrposElems[chrposElems.length - 1];
                        String[] firstElemParts = firstElem.split("-");
                        String[] lastElemParts = lastElem.split("-");

                        chrstart = firstElemParts[0];
                        chrend = lastElemParts[1];

                    } else {
                        chrposElems = chrpos.split("-");
                        if (chrposElems.length > 1) {
                            chrstart = chrposElems[0];
                            chrend = chrposElems[1];
                        } else {
                            chrstart = "-";
                            chrend = "-";
                        }
                    }
//                    if (!elems[2].equals("-1")) {
                        out.writeln(platform + "\t" + platformId + "\t" + elems[4] + "\t" + elems[2] + "\t" + chrstart + "\t" + chrend + "\t" + elems[0] + "\t" + elems[1]);
//                    }
                }


                elems = in.readLineElems(TextFile.tab);
            }


            in.close();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
