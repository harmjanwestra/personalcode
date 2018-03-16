/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class TableForLudeReplication {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here


        // /Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/MetaAnalysisAllSNPProbePairs/eQTLsFilteredForMetaCisFXRemovedFDR0.05.txtht12v3.txt
        // String transeqtlfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/PhD/2013-2012-2011-eQTLMeta/Pack-Rebuttal1/Supp/TransFXFDR0.5.txt";
        String transeqtlfile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/DataForRebuttal/2013-03-Trans40PCsCisFXNotRegressedOut/FDR0.5SNPProbes/eQTLs.txt-WSampleSize.txt.gzFilteredForMetaCisFXRemovedFDR0.05.txtht12v3.txt";
        String outfile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/PhD/2013-2012-2011-eQTLMeta/Pack-Rebuttal2/KoraReplication/MetaAnalysisTranseQTLsFDR0.05CisFXNotRegressedThatHaveBeenTestedInKORA.txt";
        String korafile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2FilteredForFDR0.05SNPProbes/eQTLsFDR.txt-HT12v3.gz";
        String blooddata = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";

        try {
            TextFile tf = new TextFile(transeqtlfile, TextFile.R);


            HashMap<Pair<String, String>, String> transeqtlsr2 = new HashMap<Pair<String, String>, String>();
            HashMap<Pair<String, String>, Double> transeqtlsz = new HashMap<Pair<String, String>, Double>();
            HashMap<Pair<String, String>, Pair<String, String>> transalleles = new HashMap<Pair<String, String>, Pair<String, String>>();

            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                String r2 = elems[elems.length - 1];
                String snp = elems[1];
                String probe = elems[4];

                String alleles = elems[8];
                String alleleAssessed = elems[9];
                Double z = Double.parseDouble(elems[10]);

                Pair<String, String> e = new Pair<String, String>(snp, probe);

                transalleles.put(e, new Pair<String, String>(alleles, alleleAssessed));
                transeqtlsz.put(e, z);
                System.out.println("Adding: " + e.toString());
                transeqtlsr2.put(e, r2);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(transeqtlsr2.size() + " transeqtls loaded");

            TriTyperGenotypeData ds = new TriTyperGenotypeData(blooddata);

            SNPLoader loader = ds.createSNPLoader();


            TextFile tfout = new TextFile(outfile, TextFile.W);
            TextFile kora = new TextFile(korafile, TextFile.R);
            kora.readLine();
            elems = kora.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];

                String samplesize = elems[13];
                double rkora = Double.parseDouble(elems[17]);
                double r2kora = rkora * rkora;

                String fdr = elems[elems.length - 1];

                Pair<String, String> e = new Pair<String, String>(snp, probe);

                String r2 = transeqtlsr2.get(e);
                if (r2 != null) {

                    String alleles = elems[8];
                    String alleleAssessed = elems[9];
                    Double z = Double.parseDouble(elems[10]);
                    Pair<String, String> transalleleforeqtl = transalleles.get(e);
                    Double metaz = transeqtlsz.get(e);
                    Boolean flip = BaseAnnot.flipalleles(transalleleforeqtl.getLeft(), transalleleforeqtl.getRight(), alleles, alleleAssessed);

                    if (flip == null) {
                        System.out.println("Incompatible alleles: " + transalleleforeqtl.getLeft() + "\t" + alleles);
                    } else if (flip) {
                        z *= -1;
                    }

                    boolean identicalDir = true;
                    if (z >= 0 && metaz < 0) {
                        identicalDir = false;
                    } else if (z < 0 && metaz >= 0) {
                        identicalDir = false;
                    }

                    Integer snpid = ds.getSnpToSNPId().get(snp);
                    if (snpid != null) {
                        SNP snpobj = ds.getSNPObject(snpid);
                        loader.loadGenotypes(snpobj);
                        double maf = snpobj.getMAF();
                        tfout.writeln(snp + "\t" + probe + "\t" + maf + "\t" + r2 + "\t" + samplesize + "\t" + r2kora + "\t" + fdr + "\t" + identicalDir);
                        snpobj.clearGenotypes();
                    } else {
                        System.out.println("SNP ID NULL: " + snp);
                    }
                } else {
                    System.out.println("R2 is null for eqtl: " + e.toString());
                }

                elems = kora.readLineElems(TextFile.tab);
            }

            loader.close();
            kora.close();
            tfout.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
