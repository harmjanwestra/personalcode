/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harm-jan
 */
public class IterativeConditionalAnalysisTabulator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String indir = "/Volumes/iSnackHD/Data/Projects/TonuEsko/eQTLResults/2013-01-16-HeightSNPsProbesWithin1MB-Cis-Run2/";
            int iteration = 7;
            double r2threshold = 0d;
            String queryfile = "/Volumes/iSnackHD/Data/Projects/TonuEsko/2013-02-11-snps.txt";
            String outdir = "/Volumes/iSnackHD/Data/Projects/TonuEsko/output/";
            String reference = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
            IterativeConditionalAnalysisTabulator t = new IterativeConditionalAnalysisTabulator();
            double fdr = 0.05;
            t.run(indir, iteration, queryfile, outdir, reference, r2threshold, fdr, IterativeConditionalAnalysisTabulator.type.ALL);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public enum type {

        PROBE, ALL
    };

    public void run(String indir, int nrIterations, String queryFile, String outdir, String reference, double r2threshold, double fdr, type t) throws IOException {
        // read query SNPs
        TextFile tf = new TextFile(queryFile, TextFile.R);
        HashSet<String> querySNPs = new HashSet<String>();
        ArrayList<String> querySNPsOrdered = tf.readAsArrayList();
        querySNPs.addAll(querySNPsOrdered);
        tf.close();

        HashMap<String, ArrayList<EQTL>> eqtls = new HashMap<String, ArrayList<EQTL>>();

        TriTyperGenotypeData ds = new TriTyperGenotypeData();
        ds.load(reference);
        SNPLoader loader = ds.createSNPLoader();
        // gather results over all iterations

        DetermineLD ldcalc = new DetermineLD();

        for (int it = 1; it < nrIterations + 1; it++) {
            String file = indir + "Iteration" + it + "/eQTLsProbesFDR" + fdr + ".txt";
            if (t == type.ALL) {
                file = indir + "Iteration" + it + "/eQTLsFDR" + fdr + ".txt";
            }

            TextFile tf2 = new TextFile(file, TextFile.R);
            String header = tf2.readLine();
            String[] elems = tf2.readLineElemsReturnObjects(TextFile.tab);
            while (elems != null) {
                // process the data..
                Double p = Double.parseDouble(elems[0]);
                String eSNP = elems[1];
                if (querySNPs.contains(eSNP)) {
                    // found a hit :-D
                    ArrayList<EQTL> eqtlsForSNP = eqtls.get(eSNP);
                    if (eqtlsForSNP == null) {
                        eqtlsForSNP = new ArrayList<EQTL>();
                    }
                    EQTL e = new EQTL();
                    e.eSNP = eSNP;
                    e.gene = elems[eQTLTextFile.HUGO];
                    e.iteration = it;
                    e.ld = 1;
                    e.probe = elems[4];
                    e.pval = p;
                    eqtlsForSNP.add(e);
                    eqtls.put(eSNP, eqtlsForSNP);
                } else {

                    // iterate over all query SNPs, check whether on the same chromosome
                    // if so, check whether in LD
                    Integer snpID = ds.getSnpToSNPId().get(eSNP);

                    if (snpID == null) {
                        System.out.println(eSNP + " is not in HapMap2?");

                    } else {
                        byte chr = -1;
                        SNP snpObj = null;

                        chr = ds.getChr(snpID);
                        snpObj = ds.getSNPObject(snpID);
                        loader.loadGenotypes(snpObj);

                        for (String snp : querySNPs) {

                            Integer snpID2 = ds.getSnpToSNPId().get(snp);
                            byte chr2 = ds.getChr(snpID2);

                            if (chr == chr2) {

                                SNP snpObj2 = ds.getSNPObject(snpID2);
                                loader.loadGenotypes(snpObj2);

                                double r2 = ldcalc.getRSquared(snpObj, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                if (r2 >= r2threshold) {
                                    ArrayList<EQTL> eqtlsForSNP = eqtls.get(snp);
                                    if (eqtlsForSNP == null) {
                                        eqtlsForSNP = new ArrayList<EQTL>();
                                    }
                                    EQTL e = new EQTL();
                                    e.eSNP = eSNP;
                                    e.gene = elems[eQTLTextFile.HUGO];
                                    e.iteration = it;
                                    e.ld = r2;
                                    e.probe = elems[4];
                                    e.pval = p;
                                    eqtlsForSNP.add(e);
                                    eqtls.put(eSNP, eqtlsForSNP);
                                }
                                snpObj2.clearGenotypes();
                            }
                        }
                        snpObj.clearGenotypes();
                    }
                }
                elems = tf2.readLineElemsReturnObjects(TextFile.tab);
            }
            tf2.close();
        }
        loader.close();

        int snpctr = 0;
        outdir = Gpio.formatAsDirectory(outdir);
        Gpio.createDir(outdir);
        String outfilename = outdir + "IterationSummary-R" + r2threshold + "-FDR" + fdr + ".txt";
        if (t == type.ALL) {
            outfilename = outdir + "IterationSummary-R" + r2threshold + "-FDR" + fdr + "-AllFX.txt";
        }
        TextFile outEQTLs = new TextFile(outfilename, TextFile.W);
        TextFile outNoEQTLs = new TextFile(outdir + "NoEQTLs.txt", TextFile.W);

        String header = "Nr\tHeightSNP";
        for (int it = 1; it < nrIterations + 1; it++) {
            header += "\tIt" + it + "Pval\tIt" + it + "eSNP\tIt" + it + "LD\tIt" + it + "Probe\tIt" + it + "Gene";
        }
        outEQTLs.writeln(header);

        for (String querySNP : querySNPsOrdered) {
            ArrayList<EQTL> eqtlsforSNP = eqtls.get(querySNP);
            if (eqtlsforSNP == null) {
                // no eQTLs
                String ln = snpctr + "\t" + querySNP;
                // outEQTLs.writeln(ln);
            } else {
                Collections.sort(eqtlsforSNP);
                EQTL[][] sortedEQTLs = new EQTL[nrIterations][eqtlsforSNP.size()];
                for (int it = 1; it < nrIterations + 1; it++) {
                    int ctr = 0;
                    for (EQTL e : eqtlsforSNP) {
                        if (e.iteration == it) {
                            sortedEQTLs[it - 1][ctr] = e;
                            ctr++;
                        }
                    }
                }

                // now we've finally sorted that shit out, output everything in one big momma of a table
                for (int e = 0; e < sortedEQTLs[0].length; e++) {
                    String ln = snpctr + "\t" + querySNP;
                    int nrEQTLsForRow = 0;
                    for (int it = 0; it < nrIterations; it++) {
                        if (sortedEQTLs[it][e] != null) {
                            EQTL eq = sortedEQTLs[it][e];
                            ln += "\t" + eq.pval + "\t" + eq.eSNP + "\t" + eq.ld + "\t" + eq.probe + "\t" + eq.gene;
                            nrEQTLsForRow++;
                        } else {
                            ln += "\t-\t-\t-\t-\t-";
                        }
                    }
                    if (nrEQTLsForRow == 0) {
                        break;
                    } else {
                        // print row
                        outEQTLs.writeln(ln);
                    }
                }

            }
            snpctr++;
        }
        outEQTLs.close();
        outNoEQTLs.close();
    }
    // Inner class Test2

    class EQTL implements Comparable<EQTL> {

        public int iteration;
        public String probe;
        public String gene;
        public String eSNP;
        public double pval;
        public double ld;
        public double fdr;

        @Override
        public int hashCode() {
            int hash = 5;
            hash = 29 * hash + this.iteration;
            hash = 29 * hash + (this.probe != null ? this.probe.hashCode() : 0);
            hash = 29 * hash + (this.gene != null ? this.gene.hashCode() : 0);
            hash = 29 * hash + (this.eSNP != null ? this.eSNP.hashCode() : 0);
            hash = 29 * hash + (int) (Double.doubleToLongBits(this.pval) ^ (Double.doubleToLongBits(this.pval) >>> 32));
            return hash;
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final EQTL other = (EQTL) obj;
            if ((this.probe == null) ? (other.probe != null) : !this.probe.equals(other.probe)) {
                return false;
            }
            if ((this.eSNP == null) ? (other.eSNP != null) : !this.eSNP.equals(other.eSNP)) {
                return false;
            }
            if (Double.doubleToLongBits(this.pval) != Double.doubleToLongBits(other.pval)) {
                return false;
            }
            return true;
        }

        @Override
        public int compareTo(EQTL o) {
            if (o.equals(this)) {
                return 0;
            } else if (o.pval > pval) {
                return -1;
            } else {
                return 1;
            }
        }
    }
}
