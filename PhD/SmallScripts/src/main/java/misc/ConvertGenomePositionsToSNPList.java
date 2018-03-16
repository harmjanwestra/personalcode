/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ConvertGenomePositionsToSNPList {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String infile = "/Volumes/iSnackHD/Data/Projects/JavierGutierrez/2013-07-25-ListHarmJan.txt";
        String outfile = "/Volumes/iSnackHD/Data/Projects/JavierGutierrez/2013-07-25-ListHarmJan-SNPs.txt";
        String snpmapfile = "/Library/WebServer/Documents/Exchange/2013-08-05-SNPMappings.txt.gz";
        try {
            TextFile in = new TextFile(infile, TextFile.R);
            String[] elems = in.readLineElems(TextFile.tab);

            HashSet<Pair<Integer, Integer>> pairs = new HashSet<Pair<Integer, Integer>>();
            while (elems != null) {
                Integer chr = Integer.parseInt(elems[0]);
                Integer position = Integer.parseInt(elems[3]);
                pairs.add(new Pair<Integer, Integer>(chr, position));
                elems = in.readLineElems(TextFile.tab);
            }
            in.close();

            TextFile out = new TextFile(outfile, TextFile.W);
            TextFile tf2 = new TextFile(snpmapfile, TextFile.R);
            elems = tf2.readLineElems(TextFile.tab);
            while (elems != null) {
                Integer chr = Integer.parseInt(elems[0]);
                Integer pos = Integer.parseInt(elems[1]);
                String snp = elems[2];
                
                Pair<Integer, Integer> p = new Pair<Integer, Integer>(chr, pos);
                if (pairs.contains(p)) {
                    out.writeln(snp);
                }

                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
