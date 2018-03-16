/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.snps;

import eqtlmappingpipeline.metaqtl3.FDR;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class SelectSimilarSNPsFromeQTLFiles {

    TriTyperGenotypeData ds = null;
    HashSet<String> allStartSNPs = null;
    double mafGranularity = 0.05;
    double distanceGranularity = 10000;
    private Features ensembl;
    private HashMap<String, Integer> querySNPmafBins;
    private HashMap<String, Integer> querySNPdistBins;
    private SNPLoader loader;
    private HashMap<Integer, ArrayList<String>> mafBins;
    private HashMap<Integer, ArrayList<String>> distBins;
    private HashSet<String> querySNPs;
    private HashSet<String> eQTLSNPs;

    public static void main(String[] args) {
        // TODO code application logic here
        SelectSimilarSNPsFromeQTLFiles selector = new SelectSimilarSNPsFromeQTLFiles();
        SelectSimilarSNPs q = new SelectSimilarSNPs();
        String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
        reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
//            String reference = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        String querySNPFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/Repruning/GwasSignificant-Clumped.txt";
        String cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
        // cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt";
//        String outdir = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/EQTLMetaRebuttal/2013-02-GWASSNPTransEQTLEnrichment-Run2/";
        String geneAnnotation = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl54_HG18/structures/structures_b54.txt";
        int nrOfSets = 100;
        double r2 = 0.2;
//        q.run(reference, querySNPFile, cisEQTLFile, outdir, nrOfSets, r2, geneAnnotation);

        try {
            for (int i = 0; i < 101; i++) {
                String indir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-Rebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/Set-" + (i - 1) + "/";
                String outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/eQTLs/Set-" + (i - 1) + "/";
                if (i == 0) {
                    indir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-Rebuttal/2013-02-GWASSNPTransEQTLEnrichment/eQTLs/QuerySet/";
                    outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/eQTLs/QuerySet/";
                }
                Gpio.createDir(outdir);

                selector.run(reference, geneAnnotation, cisEQTLFile, querySNPFile, indir, outdir, 1);


            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    private HashMap<String, String> replaceQueryByTheseProxies;

    public void run(String reference, String ensemblFeatures, String cisEQTLFile, String startSNPs, String indir, String outdir, int nrpermutations) throws IOException {
        if (ds == null) {
            ds = new TriTyperGenotypeData(reference);

            ensembl = new Features();
            ensembl.loadAnnotation(ensemblFeatures);

            // bins for query SNPs
            querySNPmafBins = new HashMap<String, Integer>();
            querySNPdistBins = new HashMap<String, Integer>();
            mafBins = new HashMap<Integer, ArrayList<String>>();
            distBins = new HashMap<Integer, ArrayList<String>>();
            loader = ds.createSNPLoader();

            TextFile queryFile = new TextFile(startSNPs, TextFile.R);
            querySNPs = new HashSet<String>();
            HashMap<String, ArrayList<String>> queryProxies = new HashMap<String, ArrayList<String>>();
            queryFile.readLine(); // skip that ugly header
            String ln = queryFile.readLine();
            while (ln != null) {
                while (ln.contains("  ")) {
                    ln = ln.replaceAll("  ", " ");
                }
                String[] qelems = ln.split(" ");

                if (qelems.length > 12) {
                    String lead = qelems[3];
                    String[] proxies = qelems[12].split(",");
                    querySNPs.add(lead);
                    System.out.println(lead);
                    ArrayList<String> proxiesForSNP = new ArrayList<String>();
                    for (String s : proxies) {
                        if (!s.equals("NONE")) {
                            s = s.replace("(1)", "");
                            proxiesForSNP.add(s);
                        }
                    }
                    queryProxies.put(lead, proxiesForSNP);
                }
                ln = queryFile.readLine();
            }
            queryFile.close();

            System.out.println(querySNPs.size() + " SNPs in query...");
            if (querySNPs.isEmpty()) {
                System.out.println("0 SNPs in query... ");
                System.exit(0);
            }

            // use the cis-eQTLs to determine the distance to the gene for each SNP.
            eQTLSNPs = new HashSet<String>();
            HashMap<String, Chromosome> chromosomes = ensembl.getChromosomeHash();

            if (cisEQTLFile != null) {
                TextFile efile = new TextFile(cisEQTLFile, TextFile.R);
                efile.readLine();
                String[] elems = efile.readLineElemsReturnObjects(TextFile.tab);
                int lnctr = 0;
                while (elems != null) {
                    String snp = elems[1];
                    eQTLSNPs.add(snp);
                    lnctr++;
                    if (lnctr % 1000000 == 0) {
                        System.out.println(lnctr + " eQTLs parsed");
                    }
                    elems = efile.readLineElemsReturnObjects(TextFile.tab);
                }
                efile.close();
                System.out.println("Loaded eQTL SNPs: " + eQTLSNPs.size());
            } else {
                eQTLSNPs.addAll(Arrays.asList(ds.getSNPs()));
            }
            HashSet<String> finalQuerySNPs = new HashSet<String>();
            for (String snp : querySNPs) {
                if(eQTLSNPs.contains(snp)){
                    finalQuerySNPs.add(snp);
                } else {
                    System.out.println(snp+" not in eQTL list");
                }
            }
            
            querySNPs = finalQuerySNPs;
            System.out.println(querySNPs.size()+" query snps remain.");

            int snpctr = 0;
            System.out.println("Binning your snps..");
            replaceQueryByTheseProxies = new HashMap<String, String>();
            for (String snp : eQTLSNPs) {
                Integer mafbin = SelectSimilarSNPs.getMAFBin(ds, loader, snp, mafGranularity);

                if (snp.equals("rs10440833")) {
                    System.out.println("rs10440833 " + mafbin);
                }

//                if (mafbin == null && querySNPs.contains(snp)) {
//                    ArrayList<String> proxies = queryProxies.get(snp);
//                    int ctr = 0;
//                    String selected = null;
//                    System.out.println("Query: " + snp + " not in ref. Selecting from: " + proxies.size() + " proxies.");
//                    while (ctr < proxies.size() && mafbin == null) {
//                        String proxySNP = proxies.get(ctr);
//                        mafbin = SelectSimilarSNPs.getMAFBin(ds, loader, proxySNP, mafGranularity);
//                        selected = proxySNP;
//                        ctr++;
//                    }
//                    if (mafbin != null) {
//                        System.out.println("replacing " + snp + " with " + selected);
//                        replaceQueryByTheseProxies.put(snp, selected);
//                        snp = selected;
//                    }
//                }

                if (mafbin != null) {
                    ArrayList<String> selectedSNPs = mafBins.get(mafbin);
                    if (selectedSNPs == null) {
                        selectedSNPs = new ArrayList<String>();
                    }

                    if (querySNPs.contains(snp)) {
                        querySNPmafBins.put(snp, mafbin);
                    } else {
                        selectedSNPs.add(snp);
                        mafBins.put(mafbin, selectedSNPs);
                    }

                    Integer distanceBin = SelectSimilarSNPs.getDistanceBin(ds, snp, chromosomes, distanceGranularity);
                    if (distanceBin != null) {
                        selectedSNPs = distBins.get(distanceBin);
                        if (selectedSNPs == null) {
                            selectedSNPs = new ArrayList<String>();
                        }

                        if (querySNPs.contains(snp)) {
                            querySNPdistBins.put(snp, distanceBin);
                        } else {
                            selectedSNPs.add(snp);
                            distBins.put(distanceBin, selectedSNPs);
                        }
                    }

                }

                snpctr++;
                if (snpctr % 100000 == 0) {
                    System.out.println(snpctr + " snps");
                }
            }


            loader.close();
        }

        // read eQTL file, select matching SNPs
        int minNrEQTLs = Integer.MAX_VALUE;
        for (int perm = 0; perm < nrpermutations + 1; perm++) {

            String infile = "";
            if (perm == 0) {
                infile = "eQTLs.txt";
                if (!Gpio.exists(infile)) {
                    infile = "eQTLs.txt.gz";
                }
            } else {
                infile = "PermutedEQTLsPermutationRound" + perm + ".txt.gz";
            }

            System.out.println("Parsing data from: " + indir);
            TextFile tf = new TextFile(indir + infile, TextFile.R);
            tf.readLine();
            HashSet<String> snpsInFile = new HashSet<String>();
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                snpsInFile.add(elems[1]);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println(snpsInFile.size() + " SNPs to select from");

            HashSet<String> selectedSubset = new HashSet<String>();
            int nrQuerySNPsNotFound = 0;
            for (String snp : querySNPs) {
                Integer mafbin = querySNPmafBins.get(snp);
                Integer distbin = querySNPdistBins.get(snp);

//                if (mafbin == null) {
//                    snp = replaceQueryByTheseProxies.get(snp);
//                    if (snp != null) {
//                        mafbin = querySNPmafBins.get(snp);
//                        distbin = querySNPdistBins.get(snp);
//                    }
//                }

                if (mafbin != null && distbin != null) {
                    ArrayList<String> mafbins = mafBins.get(mafbin);
                    ArrayList<String> distbins = distBins.get(distbin);
                    HashSet<String> mafbinhash = new HashSet<String>();
                    HashSet<String> disbinhash = new HashSet<String>();

                    mafbinhash.addAll(mafbins);
                    disbinhash.addAll(distbins);

                    for (String snp2 : snpsInFile) {
                        if (!selectedSubset.contains(snp2) && mafbinhash.contains(snp2) && disbinhash.contains(snp2)) {
                            selectedSubset.add(snp2);
                            break;
                        }
                    }
                } else {
                    System.out.println(snp + " query snp got thrown out... " + distbin + "\t" + mafbin);
                    nrQuerySNPsNotFound++;
                }


            }


            System.out.println(selectedSubset.size() + " snps matched to input: " + querySNPs.size() + "\tNr Not in Reference: " + nrQuerySNPsNotFound);

            if (selectedSubset.isEmpty() || ((double) selectedSubset.size() / querySNPs.size()) < 0.99) {
                System.out.println("ARGH: not enough SNPs to choose from... ");

                System.exit(0);
            }

            // now filter the eQTL files for the new subset.. subsequently run FDR.
            TextFile tfout = new TextFile(outdir + infile, TextFile.W);
            TextFile tfin = new TextFile(indir + infile, TextFile.R);
            tfout.writeln(tfin.readLine());
            elems = tfin.readLineElems(TextFile.tab);
            int ctr = 0;
            while (elems != null) {

                if (selectedSubset.contains(elems[1])) {
                    tfout.writeln(Strings.concat(elems, Strings.tab));
                    ctr++;
                }
                elems = tfin.readLineElems(TextFile.tab);

            }

            tfin.close();
            tfout.close();

            if (ctr < minNrEQTLs) {
                minNrEQTLs = ctr;
            }

        }

        FDR.calculateFDR(outdir, nrpermutations, minNrEQTLs, 0.05, false, null, null);


    }
}
