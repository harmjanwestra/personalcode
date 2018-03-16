/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.util.HashMap;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harm-jan
 */
public class CreateEQTLTable {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            GWASCatalog c = new GWASCatalog();
            c.read("d:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2011-09-24-gwascatalog.txt");
            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeToHT12 = pb.getProbeTranslation("d:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-04-23-Annot.txt", "Probe", "HumanHT-12_V3_0_R2_11283641_A.txt");
            String inFileName = "D:\\Skydrive\\latesteQTLs\\cisFDR0.5.txt.gz";
            TextFile inFile = new TextFile(inFileName, TextFile.R);
            TextFile outFile = new TextFile(inFileName+"-SuppTableFormat.txt.gz", TextFile.W);
            String header = "P-value	SNP Name	SNP Chromosome	SNP Chromosome Position	Probe Name	Probe Chromosome	Probe Center Chromosome Position	SNP Type	Z-score effect allele	Meta-analysis z-score	Datasets where SNP-probe combination passes QC	DatasetsZScores	Number of samples per dataset	Gene symbol (HUGO)	FDR	Associated GWAS traits";
            outFile.writeln(header);
            inFile.readLine();

            double threshold = 0.00352058461224396;
            String[] elems = inFile.readLineElems(TextFile.tab);

            while (elems != null) {

                double pval = Double.parseDouble(elems[0]);
                if (pval <= threshold) {
                    String[] output = new String[16];
                    output[0] = elems[0];
                    output[1] = elems[1];
                    output[2] = elems[2];
                    output[3] = elems[3];
                    output[4] = probeToHT12.get(elems[4]);
                    output[5] = elems[5];
                    output[6] = elems[6];
                    output[7] = elems[eQTLTextFile.SNPTYPE];
                    output[8] = elems[eQTLTextFile.ASESSEDALLELE];
                    output[9] = elems[eQTLTextFile.METAZ];
                    output[10] = elems[eQTLTextFile.DATASETNAMES];
                    output[11] = elems[eQTLTextFile.DATASETZSCORE];
                    output[12] = elems[eQTLTextFile.DATASETSIZE];
                    output[13] = elems[eQTLTextFile.HUGO];
                    output[14] = elems[elems.length - 1];
                    GWASSNP snpObj = c.getSnpToObj().get(elems[1]);
                    if (snpObj != null) {
                        GWASTrait[] traits = snpObj.getAssociatedTraitsArray();
                        output[15] = Strings.concat(traits, Strings.semicolon);
                    } else {
                        output[15] = "-";
                    }
                    outFile.writeln(Strings.concat(output, Strings.tab));
                }





                elems = inFile.readLineElems(TextFile.tab);
            }
            inFile.close();
            outFile.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
