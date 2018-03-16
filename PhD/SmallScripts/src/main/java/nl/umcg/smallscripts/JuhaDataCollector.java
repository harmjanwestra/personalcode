/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class JuhaDataCollector {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String efile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";
            String anfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt";
            String qfile = "/Volumes/iSnackHD/Data/Projects/JuhaKarjalainen/2013-01-09-cisEQTLLookupForGWASSNPs/SNPsNewMeta.txt";
            String ofile = "/Volumes/iSnackHD/Data/Projects/JuhaKarjalainen/2013-01-09-cisEQTLLookupForGWASSNPs/SNPsNewMeta-eQTLs.txt";
            String ttdata = "/Data/GeneticalGenomicsDatasets/BloodHT12ImputeTriTyper/";

            TextFile tf = new TextFile(qfile, TextFile.R);
            HashSet<String> query = new HashSet<String>();
            query.addAll(tf.readAsArrayList());
            tf.close();

            TextFile tf3 = new TextFile(anfile, TextFile.R);
            HashMap<String, String> probeToEnsembl = (HashMap<String, String>) tf3.readAsHashMap(0, 4);
            tf3.close();

            TextFile tf2 = new TextFile(efile, TextFile.R);
            TextFile tf4 = new TextFile(ofile, TextFile.W);
            String[] elems = tf2.readLineElems(TextFile.tab);

            HashMap<String, String> topSNPs = new HashMap<String, String>();
            HashMap<String, String> topSNPPvals = new HashMap<String, String>();
            HashSet<String> visitedProbes = new HashSet<String>();

            TriTyperGenotypeData ds = new TriTyperGenotypeData();
            ds.load(ttdata);
            SNPLoader loader = ds.createSNPLoader();
            DetermineLD ldcalc = new DetermineLD();

            tf4.writeln("QuerySNPeQTLPVal\tQuerySNP\tTopPVal\tTopeSNP\tLD(R2In1240BloodSamples)\tProbe\tEnsembl");
            while (elems != null) {
                String snp = elems[1];
                if (!visitedProbes.contains(elems[4])) {
                    topSNPs.put(elems[4], snp);
                    topSNPPvals.put(elems[4], elems[0]);
                    visitedProbes.add(elems[4]);
                }
                if (query.contains(snp)) {
                    Integer snpId = ds.getSnpToSNPId().get(snp);
                    Integer snpId2 = ds.getSnpToSNPId().get(topSNPs.get(elems[4]));
                    Double ld = null;
                    if (snpId != null && snpId2 != null) {

                        SNP snpObj1 = ds.getSNPObject(snpId);
                        SNP snpObj2 = ds.getSNPObject(snpId2);
                        
                        
                        loader.loadGenotypes(snpObj1);
                        loader.loadGenotypes(snpObj2);
                        
                        ld = ldcalc.getRSquared(snpObj1, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                        
                        snpObj1.clearGenotypes();
                        snpObj2.clearGenotypes();
                    
                    }
                    String ann = probeToEnsembl.get(elems[4]);
                    tf4.writeln(elems[0] + "\t" + snp + "\t" + topSNPPvals.get(elems[4]) + "\t" + topSNPs.get(elems[4]) + "\t" + ld + "\t" + elems[4] + "\t" + ann);
                }
                elems = tf2.readLineElems(TextFile.tab);
            }
            tf2.close();
            tf4.close();
            loader.close();


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
