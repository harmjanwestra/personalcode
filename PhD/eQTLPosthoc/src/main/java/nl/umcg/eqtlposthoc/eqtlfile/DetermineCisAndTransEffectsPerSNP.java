/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DetermineCisAndTransEffectsPerSNP {

    /**
     * @param args the command line arguments
     * 
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String cis = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-05-MetaAnalysisFinal/cis/2012-05-09-EGCUT+SHIP+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/PostQC/eQTLsFDR0.05.txt";
            String trans = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-05-MetaAnalysisFinal/trans/2012-05-02-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt";
            String probetrans = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            //String useMapping = "HumanHT-12_V3_0_R2_11283641_A.txt"; 
//            String useMapping = "HumanHT-12_V4_0_R1_15002873_B.txt";
             String useMapping = "H8v2ConvToHT12";
            String outfile = "/Volumes/iSnackHD/cisAndTransFX-H8v2.txt";
            run(cis, trans, probetrans, useMapping, outfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void run(String cis, String trans, String probeTranslation, String useMapping, String outfile) throws IOException {

// load trans effects
        ArrayList<String> uniqueTransSNPs = new ArrayList<String>();
        HashMap<String, ArrayList<String>> transFX = new HashMap<String, ArrayList<String>>();
        TextFile tf = new TextFile(trans, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab); // header
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            ArrayList<String> transprobes = transFX.get(snp);
            if (transprobes == null) {
                transprobes = new ArrayList<String>();
                uniqueTransSNPs.add(snp);
            }
            transprobes.add(probe);
            transFX.put(snp, transprobes);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

// load cis effects
        tf = new TextFile(cis, TextFile.R);
        elems = tf.readLineElems(TextFile.tab); // header
        HashMap<String, ArrayList<String>> cisFX = new HashMap<String, ArrayList<String>>();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String probe = elems[4];
            if (transFX.containsKey(snp)) {
                ArrayList<String> cisprobes = cisFX.get(snp);
                if (cisprobes == null) {
                    cisprobes = new ArrayList<String>();
                }
                cisprobes.add(probe);
                cisFX.put(snp, cisprobes);
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

// load probe translation table
        HashMap<String, String> probeToProbe = new HashMap<String, String>();

        tf = new TextFile(probeTranslation, TextFile.R);
        elems = tf.readLineElems(TextFile.tab); // header
        int mappingCol = -1;
        for (int i = 0; i < elems.length; i++) {
            if (elems[i].equals(useMapping)) {
                mappingCol = i;
            }
        }
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String metaprobe = elems[0];
            if (elems.length >= mappingCol) {
                String ht12v3probe = elems[mappingCol];
                probeToProbe.put(metaprobe, ht12v3probe);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

// output results
        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeln("SNPName\tcis_probe\ttrans_probe");
        for (String snp : uniqueTransSNPs) {
            ArrayList<String> cisprobes = cisFX.get(snp);
            ArrayList<String> transprobes = transFX.get(snp);
            if (cisprobes != null) {
                for (String cisf : cisprobes) {
                    String cisProbeArrayAddress = probeToProbe.get(cisf);
                    for (String transf : transprobes) {
                        String transProbeArrayAddress = probeToProbe.get(transf);

//                    System.out.println(snp + "\t" + probeToProbe.get(cisf) + "\t" + probeToProbe.get(transf));
                        if (cisProbeArrayAddress.equals("-") || transProbeArrayAddress.equals("-")) {
                        } else {
                            out.writeln(snp + "\t" + cisProbeArrayAddress + "\t" + transProbeArrayAddress);
                        }
                    }
                }
            }
        }
        out.close();
    }
}
