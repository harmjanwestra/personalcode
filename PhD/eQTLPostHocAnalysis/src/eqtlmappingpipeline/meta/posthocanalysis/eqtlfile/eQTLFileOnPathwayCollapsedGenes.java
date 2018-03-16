/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFileOnPathwayCollapsedGenes {

    private static GWASCatalog catalog;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            catalog = new GWASCatalog();

            catalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");
            TextFile eqtlfile = new TextFile("/Volumes/iSnackHD/Data/Projects/LudeFranke/eqtlresults/2012-06-14-CISTRANS-Reactome_BPMean-DiseaseSNPs-CisRegressedOut-META/eQTLProbesFDR0.05.txt", TextFile.R);
            String[] elems = eqtlfile.readLineElems(TextFile.tab);
            elems = eqtlfile.readLineElems(TextFile.tab);
            TextFile out = new TextFile("/Volumes/iSnackHD/SkyDrive/PathwaysInfluencedByTraitAssocSNPs.txt" , TextFile.W);
            while (elems != null) {
                String snp = elems[1];
                String pw = elems[4];

                GWASSNP snpObj = catalog.getSnpToObj().get(snp);
                GWASTrait[] traits = snpObj.getAssociatedTraits().toArray(new GWASTrait[0]);
                ArrayList<String> traitnames = new ArrayList<String>();

                for (GWASTrait trait : traits) {
                    traitnames.add(trait.getName());
                }

                String[] list = traitnames.toArray(new String[0]);

                System.out.println(snp + "\t" + elems[2] + "\t" + elems[3] + "\t" + Strings.concat(list, Strings.comma) + "\t" + pw);
                out.writeln(snp + "\t" + elems[2] + "\t" + elems[3] + "\t" + Strings.concat(list, Strings.comma) + "\t" + pw);
                

                elems = eqtlfile.readLineElems(TextFile.tab);
            }
            out.close();
            eqtlfile.close();
        } catch (IOException e) {
        }
    }
}
