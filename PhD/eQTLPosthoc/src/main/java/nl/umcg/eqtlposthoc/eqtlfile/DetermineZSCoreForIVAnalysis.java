/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author harmjan
 */
public class DetermineZSCoreForIVAnalysis {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String probemapfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String src = "Probe";
            String dst = "HumanHT-12_V3_0_R2_11283641_A.txt";
            String cis = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt";
            String trans = "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt";
            String IV = "/Volumes/iSnackHD/eQTLMetaAnalysis/DataForHanieh/cisAndTransFX.txt";
            String outfilename = "/Volumes/iSnackHD/eQTLMetaAnalysis/DataForHanieh/cisAndTransFX-WithEffectSizes-NotFiltered.txt";


            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeMap = pb.getProbeTranslation(probemapfile, src, dst);

            HashSet<Pair<String, String>> query = new HashSet<Pair<String, String>>();
            TextFile in = new TextFile(IV, TextFile.R);
            in.readLineElems(TextFile.tab);
            String[] elems = in.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[0];
                String cisprobe = elems[1];
                String transprobe = elems[2];

                query.add(new Pair<String, String>(snp, cisprobe));
                query.add(new Pair<String, String>(snp, transprobe));

                elems = in.readLineElems(TextFile.tab);
            }
            in.close();

            // load trans-eQTL file
            HashMap<Pair<String, String>, String> effects = new HashMap<Pair<String, String>, String>();

            TextFile in2 = new TextFile(cis, TextFile.R);
            in2.readLineElems(TextFile.tab);
            String[] elems2 = in2.readLineElems(TextFile.tab);
            while (elems2 != null) {
                String probe = probeMap.get(elems2[4]);
                String snp = elems2[1];
                Pair<String, String> p = new Pair<String, String>(snp, probe);
                if (query.contains(p)) {
                    String alleles = elems2[eQTLTextFile.ASESSEDALLELE - 1] + "\t" + elems2[eQTLTextFile.ASESSEDALLELE] + "\t" + elems2[eQTLTextFile.METAZ];
                    effects.put(new Pair<String, String>(snp, probe), alleles);
                }

                elems2 = in2.readLineElems(TextFile.tab);
            }
            in2.close();

            // load cis-eQTL file
            TextFile in3 = new TextFile(trans, TextFile.R);
            in3.readLineElems(TextFile.tab);
            String[] elems3 = in3.readLineElems(TextFile.tab);
            while (elems3 != null) {
                String probe = probeMap.get(elems3[4]);
                String snp = elems3[1];
                Pair<String, String> p = new Pair<String, String>(snp, probe);
                if (query.contains(p)) {
                    String alleles = elems3[eQTLTextFile.ASESSEDALLELE - 1] + "\t" + elems3[eQTLTextFile.ASESSEDALLELE] + "\t" + elems3[eQTLTextFile.METAZ];
                    effects.put(new Pair<String, String>(snp, probe), alleles);
                }
                elems3 = in3.readLineElems(TextFile.tab);
            }
            in3.close();

            TextFile out = new TextFile(outfilename, TextFile.W);
            out.writeln("snp\tcisprobe\talleles\tassessed\tz\ttransprobe\talleles\tassessed\tz");
            in = new TextFile(IV, TextFile.R);
            in.readLineElems(TextFile.tab);
            elems = in.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[0];
                String cisprobe = elems[1];
                String transprobe = elems[2];

                String cisEffect = effects.get(new Pair<String, String>(snp, cisprobe));
                String transEffect = effects.get(new Pair<String, String>(snp, transprobe));
                out.writeln(snp + "\t" + cisprobe + "\t" + cisEffect + "\t" + transprobe + "\t" + transEffect);
                elems = in.readLineElems(TextFile.tab);
            }
            in.close();
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
