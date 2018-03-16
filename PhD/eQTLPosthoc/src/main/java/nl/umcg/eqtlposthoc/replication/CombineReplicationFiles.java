/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.replication;

import java.util.ArrayList;
import java.util.HashMap;
import nl.umcg.eqtlposthoc.visualisation.CreateForrestPlots;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;

/**
 *
 * @author harm-jan
 */
public class CombineReplicationFiles {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        CombineReplicationFiles r = new CombineReplicationFiles();
        r.run();

    }

    class EQTL {

        public double pval;
        public double zscore;
        public String alleles;
        public String assessed;
        public String FDR;
    }

    public void run() {
        try {
            String originalFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/eQTLsFDR0.05.txt";
            String outFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/Replication/CombinedFile";

            ProbeTranslation pb = new ProbeTranslation();
            HashMap<String, String> probeToHT12 = pb.getProbeTranslation("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", "Probe", "HumanHT-12_V3_0_R2_11283641_A.txt");

            String[] replicationCohortNames = new String[]{"Monocytes", "B-Cells", "BSGS", "KORA", "BSGS - KORA Meta-analysis", "LCL"};
            String[] replicationFiles = new String[]{
                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BenFairFax/2012-11-25-UsedForPaper/FairFaxAllProbeNorm/FairFaxMonoFiltered/40_PC_all_flipped/eQTLsFDR.txt.gz",
                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BenFairFax/2012-11-25-UsedForPaper/FairFaxAllProbeNorm/FairFaxBCellFiltered/40_PC_all_flipped/eQTLsFDR.txt.gz",
                "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/bsgsrepl/Parsed/CisFXNotRemoved/FilteredForFDR0.05/eQTLsFDR.txt.gz",
                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2FilteredForFDR0.05SNPProbes/eQTLsFDR.txt.gz",
                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/BSGSFixedKoraMeta/MetaOut/eQTLsFDR-Meta.txt",
                "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/StrangerLCL/eQTLResults/2013-03-11-LCL-TRANS-PCACorr2/eQTLsFDR-Meta.txt.gz"
            };

            ArrayList<HashMap<String, EQTL>> eqtls = new ArrayList<HashMap<String, EQTL>>();
            for (String s : replicationFiles) {
                TextFile tf = new TextFile(s, TextFile.R);
                String[] elems = tf.readLineElems(TextFile.tab);
                elems = tf.readLineElems(TextFile.tab);
                HashMap<String, EQTL> eqtlsForDs = new HashMap<String, EQTL>();
                while (elems != null) {
                    String snpProbe = elems[1] + "-" + elems[4];
                    double pval = Double.parseDouble(elems[0]);
                    double z = Double.parseDouble(elems[eQTLTextFile.METAZ]);
                    String allele = elems[eQTLTextFile.ASESSEDALLELE];
                    String alleles = elems[eQTLTextFile.SNPTYPE];
                    EQTL e = new EQTL();
                    e.alleles = alleles;
                    e.assessed = allele;
                    e.zscore = z;
                    e.pval = pval;
                    e.FDR = elems[elems.length - 1];
                    eqtlsForDs.put(snpProbe, e);
                    elems = tf.readLineElems(TextFile.tab);
                }
                tf.close();
                eqtls.add(eqtlsForDs);
            }

            TextFile tf = new TextFile(originalFile, TextFile.R);
            TextFile tfOut = new TextFile(outFile, TextFile.W);

            String header = "SNP\tProbe\tGene\tSNP Type\tEffect Allele\tP-value\tZ-score\tFDR";
            for (int i = 0; i < replicationCohortNames.length; i++) {
                header += "\t" + replicationCohortNames[i] + " P-value\t" + replicationCohortNames[i] + " Z-score\t" + replicationCohortNames[i] + " FDR";
            }
            tfOut.writeln(header);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String snpProbe = elems[1] + "-" + elems[4];
                double pval = Double.parseDouble(elems[0]);
                double z = Double.parseDouble(elems[eQTLTextFile.METAZ]);
                String allele = elems[eQTLTextFile.ASESSEDALLELE];
                String alleles = elems[eQTLTextFile.SNPTYPE];

                String output = elems[1] + "\t" + probeToHT12.get(elems[4]) + "\t" + elems[eQTLTextFile.HUGO] + "\t" + alleles + "\t" + allele + "\t" + elems[0] + "\t" + elems[eQTLTextFile.METAZ] + "\t" + elems[elems.length - 1];


                for (int i = 0; i < replicationCohortNames.length; i++) {
                    HashMap<String, EQTL> replicationData = eqtls.get(i);
                    EQTL e = replicationData.get(snpProbe);
                    if (e == null) {
                        output += "\t-\t-\t-";
                    } else {

                        boolean flipZ = CreateForrestPlots.detmermineAlleleFlips(alleles, allele, e.alleles, e.assessed);
                        if (flipZ) {
                            output += "\t" + e.pval + "\t" + (-e.zscore) + "\t" + e.FDR;
                        } else {
                            output += "\t" + e.pval + "\t" + e.zscore + "\t" + e.FDR;
                        }


                    }

                }
                tfOut.writeln(output);


                elems = tf.readLineElems(TextFile.tab);
            }
            tfOut.close();
            tf.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
