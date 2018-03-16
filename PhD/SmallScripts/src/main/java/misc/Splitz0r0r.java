/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class Splitz0r0r {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String pbtfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        String infile = "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Meta/MetaAnalysisZScoreMatrix.txt";
        String outfile = "/Volumes/iSnackHD/Data/Projects/LudeFranke/2014-02-12-TransInteractionTerms/Meta/TopInteractionEffects.txt";
        String gwasCatalogLoc = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog.txt";

        double zthreshold = 6;

        HashMap<String, String> probeToGene = new HashMap<String, String>();
        ProbeTranslation pbt = new ProbeTranslation();

        try {

            GWASCatalog cataloig = new GWASCatalog(gwasCatalogLoc);
            probeToGene = pbt.getProbeTranslation(pbtfile, "HT12v3.txt", "Gene");
            DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(infile);

            TextFile out = new TextFile(outfile, TextFile.W);
            out.writeln("SNP\tProbe\tTransGene\tInteractionProbe\tInteractionProbeGene\tInteractionZScore\tAssociatedTraits");
            for (int col = 0; col < ds.nrCols; col++) {

                String[] elems = ds.colObjects.get(col).split("-");
                String snp = elems[0];
                String probe = elems[1];
                String gene = probeToGene.get(probe);

                GWASSNP snpObj = cataloig.getSnpToObj().get(snp);
                GWASTrait[] traits = snpObj.getAssociatedTraitsArray();
                ArrayList<String> traitsAssoc = new ArrayList<String>();
                for (GWASTrait t : traits) {
                    if (snpObj.getPValueAssociatedWithTrait(t) < 5E-8) {
                        traitsAssoc.add(t.getName());
                    }
                }

                for (int row = 0; row < ds.nrRows; row++) {
                    double z = ds.rawData[row][col];

                    String interactionbProbe = ds.rowObjects.get(row);
                    String interactionGene = probeToGene.get(interactionbProbe);

                    if (Math.abs(z) > zthreshold) {
                        out.writeln(snp + "\t" + probe + "\t" + gene + "\t" + interactionbProbe + "\t" + interactionGene + "\t" + z + "\t" + Strings.concat(traitsAssoc.toArray(new String[0]), Strings.semicolon));
                    }
                }
            }

            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
