/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.text.DecimalFormat;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harm-jan
 */
public class eQTLFileDecimalFormatter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String infile = "d:\\dropbox\\eQTLMeta\\Manuscript\\PackageFinal\\Supplementary Data\\Supplementary Table 3 - Trans-eQTLs with multiple probes per gene (FDR0.5).txt";
            String outfile = "d:\\dropbox\\eQTLMeta\\Manuscript\\PackageFinal\\Supplementary Data\\Supplementary Table 3 - Trans-eQTLs with multiple probes per gene (FDR0.5)-DecimalFormatted.txt";

            int metaZCol = -1;
            int dsZCol = -1;

            TextFile in = new TextFile(infile, TextFile.R);
            TextFile out = new TextFile(outfile, TextFile.W);

            String[] headerElems = in.readLineElems(TextFile.tab);
            for (int i = 0; i < headerElems.length; i++) {
                if (headerElems[i].equals("Meta-analysis z-score")) {
                    metaZCol = i;
                }
                if (headerElems[i].equals("Individual dataset z-scores")) {
                    dsZCol = i;
                }
            }
            out.writeln(Strings.concat(headerElems, Strings.tab));



            DecimalFormat df = new DecimalFormat("#.##");
            String[] elems = in.readLineElems(TextFile.tab);
            while (elems != null) {

                Double metaZ = Double.parseDouble(elems[metaZCol]);
                String dsZ = elems[dsZCol];

                String[] dsZElems = dsZ.split(";");
                if (dsZElems.length == 1) {
                    dsZElems = dsZ.split(",");
                }
                String[] output = new String[dsZElems.length];
                for (int i = 0; i < dsZElems.length; i++) {
                    try {
                        Double dsZVal = Double.parseDouble(dsZElems[i]);
                        if (i == 0) {
                            output[i] = df.format(dsZVal);
                        } else {
                            output[i] = " " + df.format(dsZVal);
                        }

                    } catch (NumberFormatException e) {
                        if (i == 0) {
                            output[i] = "N/A";
                        } else {
                            output[i] = " N/A";
                        }

                    }
                }
                elems[dsZCol] = Strings.concat(output, Strings.semicolon);
                elems[metaZCol] = df.format(metaZ);
                out.writeln(Strings.concat(elems, Strings.tab));
                elems = in.readLineElems(TextFile.tab);
            }

            out.close();
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
