/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.tf;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.FisherExactTest;

/**
 *
 * @author harmjan
 */
public class CountEQTLsNearTranscriptionFactors {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String allCisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
        String ciseQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";
        String transeQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05.txt";
        String ensemblAnnot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt";
        String ensemblTF = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EnsemblTF/EnsemblGenesWithGO0001071Annotation.txt";
        String gwasSNPFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-GWAS-SNPs.txt";
        
        double pval = 1.3141052405861965E-4;
        try {
            CountEQTLsNearTranscriptionFactors.run(allCisEQTLFile, ciseQTLFile, pval, transeQTLFile, gwasSNPFile, ensemblAnnot, ensemblTF);
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }

    public static void run(String allCisEQTLFile, String ciseQTLFile, double probeLevelPvalueThreshold, String transeQTLFile, String gwasSNPFile, String ensemblAnnotation, String ensemblTranscriptionFactorFile) throws IOException {

        // load TFs
        HashSet<String> ensemblTranscriptionFactors = new HashSet<String>();
        TextFile in = new TextFile(ensemblTranscriptionFactorFile, TextFile.R);
        ensemblTranscriptionFactors.addAll(in.readAsArrayList());
        in.close();
        System.out.println("Loaded: " + ensemblTranscriptionFactors.size() + " transcription factors.");


        // determine all probes that have been tested
        HashSet<String> probesThatAreTested = new HashSet<String>();
        TextFile tf = new TextFile(allCisEQTLFile, TextFile.R);
        tf.readLine();
        String[] data = tf.readLineElems(TextFile.tab);
        while (data != null) {
            probesThatAreTested.add(data[4]);
            data = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println("Loaded: " + probesThatAreTested.size() + " probes that have been tested.");


        // read probes that are tested and are TFs
        HashSet<String> otherProbesthatHaveBeenTested = new HashSet<String>();
        HashSet<String> probesThatAreEnsemblTFsAndTested = new HashSet<String>();
        HashSet<String> probesThatAreEnsemblGenesAndTested = new HashSet<String>();
        HashMap<String, String> metaIdToEnsemblId = new HashMap<String, String>();

        TextFile eFile = new TextFile(ensemblAnnotation, TextFile.R);
        eFile.readLine();
        data = eFile.readLineElems(TextFile.tab);
        while (data != null) {
            String metaId = data[0];
            String ensembl = data[5];
            String[] ensemblElems = ensembl.split(",");

            metaIdToEnsemblId.put(metaId, ensembl);
            if (probesThatAreTested.contains(metaId)) {
                if (ensemblTranscriptionFactors.contains(ensembl)) {
                    probesThatAreEnsemblTFsAndTested.add(metaId);
                } else if (ensembl.trim().length() > 0) {
                    probesThatAreEnsemblGenesAndTested.add(metaId);
                } else {
                    otherProbesthatHaveBeenTested.add(metaId);
                }
            }

            data = eFile.readLineElems(TextFile.tab);
        }
        eFile.close();
        System.out.println("Loaded: " + probesThatAreEnsemblTFsAndTested.size() + " probes that have both been tested and are transcription factors");
        System.out.println("Loaded: " + probesThatAreEnsemblGenesAndTested.size() + " probes that have both been tested and are ensembl genes");

        // load significant cis-eQTLs (Probe Level FDR threshold)
        HashSet<String> probesThatShowCisEQTLEffect = new HashSet<String>();
        HashSet<String> probesThatShowCisEQTLEffectAndAreEnsemblGene = new HashSet<String>();
        HashSet<String> probesThatShowCisEQTLEffectAndAreEnsemblTFs = new HashSet<String>();

        TextFile ciseFile = new TextFile(ciseQTLFile, TextFile.R);
        ciseFile.readLine();
        data = ciseFile.readLineElems(TextFile.tab);
        while (data != null) {
            double pval = Double.parseDouble(data[0]);
            if (pval < probeLevelPvalueThreshold) {
                if (probesThatAreEnsemblTFsAndTested.contains(data[4])) {
                    probesThatShowCisEQTLEffectAndAreEnsemblTFs.add(data[4]);
                } else if (probesThatAreEnsemblGenesAndTested.contains(data[4])) {
                    probesThatShowCisEQTLEffectAndAreEnsemblGene.add(data[4]);
                } else {
                    probesThatShowCisEQTLEffect.add(data[4]);
                }
            }
            data = ciseFile.readLineElems(TextFile.tab);
        }
        ciseFile.close();

        System.out.println("Probes w/o gene name having cis-effect: " + probesThatShowCisEQTLEffect.size());
        System.out.println("TFs probes having cis-effect: " + probesThatShowCisEQTLEffectAndAreEnsemblTFs.size());
        System.out.println("Gene probes having cis-effect: " + probesThatShowCisEQTLEffectAndAreEnsemblGene.size());

        FisherExactTest fet = new FisherExactTest();
        int nrTFsTestedButNoEQTL = probesThatAreEnsemblTFsAndTested.size() - probesThatShowCisEQTLEffectAndAreEnsemblTFs.size();
        int nrTFsTestedAndEQTL = probesThatShowCisEQTLEffectAndAreEnsemblTFs.size();
        int nrOtherTestedButNoEQTL = probesThatAreEnsemblGenesAndTested.size() - probesThatShowCisEQTLEffectAndAreEnsemblGene.size();
        int nrOtherTestedAndEQTL = probesThatShowCisEQTLEffectAndAreEnsemblGene.size();
        System.out.println("2x2:");
        System.out.println("\tTF\tNoTF\tTot");
        System.out.println("EQTL\t" + nrTFsTestedAndEQTL + "\t" + nrOtherTestedAndEQTL + "\t" + (nrTFsTestedAndEQTL + nrOtherTestedAndEQTL));
        System.out.println("nEQTL\t" + nrTFsTestedButNoEQTL + "\t" + nrOtherTestedButNoEQTL + "\t" + (nrTFsTestedButNoEQTL + nrOtherTestedButNoEQTL));
        System.out.println("Tot:\t" + (nrTFsTestedAndEQTL + nrTFsTestedButNoEQTL) + "\t" + (nrOtherTestedAndEQTL + nrOtherTestedButNoEQTL) + "\t" + ((nrTFsTestedAndEQTL + nrOtherTestedAndEQTL) + (nrTFsTestedButNoEQTL + nrOtherTestedButNoEQTL)));

        double p = fet.getFisherPValue(nrTFsTestedAndEQTL, nrTFsTestedButNoEQTL, nrOtherTestedAndEQTL, nrOtherTestedButNoEQTL);
        System.out.println("p: " + p);
    }
}
