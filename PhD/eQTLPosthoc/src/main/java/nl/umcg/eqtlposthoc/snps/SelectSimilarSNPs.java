/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.snps;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Chromosome;
import umcg.genetica.containers.Gene;
import umcg.genetica.ensembl.Features;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class SelectSimilarSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {


            // use LD pruned input (!1!)
            SelectSimilarSNPs q = new SelectSimilarSNPs();
            String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
//            String reference = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
            String querySNPFile = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/Repruning/GwasSignificant-Clumped.txt";
            String cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
            // cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt";
            String outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/2013-04-Rebuttal2/2013-04-GWASSNPTransEQTLEnrichment/";
            Gpio.createDir(outdir);
            String geneAnnotation = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl54_HG18/structures/structures_b54.txt";
            int nrOfSets = 100;

            double r2 = 0.2;
            q.run(reference, querySNPFile, cisEQTLFile, outdir, nrOfSets, r2, geneAnnotation);


//            querySNPFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-KLF14TransEQTLs/KLF14SNP.txt";
//            String confineToTheseSNPs = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TestedSNPs.txt";
//            String outfile = "";
//            q.selectSNPsWithSimilarMaf(querySNPFile, confineToTheseSNPs, reference, outfile);


//            // use LD pruned input (!1!)
//            SelectSimilarSNPs q = new SelectSimilarSNPs();
//            String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
////            String reference = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
//            String querySNPFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/SNPSelectionsForFuncAnalysis/RealData/RealData.txt";
//            String cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz";
//            // cisEQTLFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/QC/TruePositives.txt";
//            String outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/RandomlySNPSelectionsForFuncAnalysisSet2/";
//            String geneAnnotation = "/Volumes/iSnackHD/Data/GenomeSequences/Ensembl54_HG18/structures/structures_b54.txt";
//            int nrOfSets = 10;
//            double r2 = 1; // omit LD threshold
////            q.run(reference, querySNPFile, cisEQTLFile, outdir, nrOfSets, r2, geneAnnotation);
//
//
////            for (int i = 1; i < 11; i++) {
//                String confineToTheseSNPs = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PreQC/TestedSNPs.txt";
//
//                outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/RandomlySNPSelectionsForFuncAnalysis/";
////                Gpio.createDir(dir);
////                String outfile = dir + "/SNPs.txt";
//                q.run(reference, querySNPFile, cisEQTLFile, outdir, nrOfSets, r2, geneAnnotation);
////            }
////            querySNPFile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/EQTLMetaRebuttal/2013-02-KLF14TransEQTLs/KLF14SNP.txt";

        } catch (IOException e) {

            e.printStackTrace();

        }
    }

    public void run(String referenceDataset, String plinkClumpedFile, String cisEQTLFile, String outdir, int nrOfSetsToPick, double r2threshold, String geneAnnotation) throws IOException {
//         Load the SNPs we want to query - classic kind of approach
//        TextFile queryFile = new TextFile(plinkClumpedFile, TextFile.R);
//        HashSet<String> querySNPs = new HashSet<String>();
//        querySNPs.addAll(queryFile.readAsArrayList());
//        queryFile.close();
//        System.out.println(querySNPs.size() + " SNPs loaded as query");

        // rather than looking at all SNPs, we actually want to look in the clumped file.
        TextFile queryFile = new TextFile(plinkClumpedFile, TextFile.R);
        HashSet<String> querySNPs = new HashSet<String>();
        HashMap<String, ArrayList<String>> queryProxies = new HashMap<String, ArrayList<String>>();
        queryFile.readLine(); // skip that ugly header
        String ln = queryFile.readLine();
        while (ln != null) {

            while (ln.contains("  ")) {
                ln = ln.replaceAll("  ", " ");
            }


            String[] qelems = ln.split(" ");



//            System.out.println(qelems.length);
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


        // bin on distance to snp
        // bin on maf

        Features ensembl = new Features();
        ensembl.loadAnnotation(geneAnnotation);

        System.out.println(querySNPs.size() + "\tSNPs loaded from: " + plinkClumpedFile);

        double mafGranularity = 0.05;
        double distanceGranularity = 10000;

        HashMap<Integer, ArrayList<String>> mafBins = new HashMap<Integer, ArrayList<String>>();
        HashMap<Integer, ArrayList<String>> distBins = new HashMap<Integer, ArrayList<String>>();

        // bins for query SNPs
        HashMap<String, Integer> querySNPmafBins = new HashMap<String, Integer>();
        HashMap<String, Integer> querySNPdistBins = new HashMap<String, Integer>();

        // load reference dataset
        TriTyperGenotypeData ds = new TriTyperGenotypeData(referenceDataset);
        SNPLoader loader = ds.createSNPLoader();

        // use the cis-eQTLs to determine the distance to the gene for each SNP.
        HashSet<String> eQTLSNPs = new HashSet<String>();
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
            if (eQTLSNPs.contains(snp)) {
                finalQuerySNPs.add(snp);
            } else {
                System.out.println(snp + " not in eQTL list");
            }
        }

        querySNPs = finalQuerySNPs;
        System.out.println(querySNPs.size() + " query snps remain.");

        int snpctr = 0;
        System.out.println("Binning your snps..");
        finalQuerySNPs = new HashSet<String>();
        for (String snp : eQTLSNPs) {
            Integer mafbin = getMAFBin(ds, loader, snp, mafGranularity);

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

                Integer distanceBin = getDistanceBin(ds, snp, chromosomes, distanceGranularity);
                if (distanceBin != null) {
                    selectedSNPs = distBins.get(distanceBin);
                    if (selectedSNPs == null) {
                        selectedSNPs = new ArrayList<String>();
                    }

                    if (querySNPs.contains(snp)) {
                        querySNPdistBins.put(snp, distanceBin);
                        finalQuerySNPs.add(snp);
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

        System.out.println(finalQuerySNPs.size() + " query SNPs are also in eQTL file");
        querySNPs = finalQuerySNPs;




        System.out.println("Loaded details for " + distBins.size() + " bins, containing " + querySNPdistBins.size() + " SNPs.");

        System.out.println("Query distribution");
        int[] queryDistbins = new int[distBins.size()];
        int[] queryMafbins = new int[mafBins.size()];
        int nrNotPresent = 0;
        HashSet<String> snpsNotInReference = new HashSet<String>();
        for (String s : querySNPs) {
            Integer distBinNo = querySNPdistBins.get(s);
            Integer mafBinNo = querySNPmafBins.get(s);
            if (distBinNo != null && mafBinNo == null) {
                System.out.println("There is an issue with SNP: " + s + " NO MAF BUT DOES HAVE DISTANCE");
            } else if (distBinNo != null && mafBinNo != null) {

                queryDistbins[distBinNo]++;
                queryMafbins[mafBinNo]++;
            } else {
                snpsNotInReference.add(s);
                nrNotPresent++;
            }
        }

        System.out.println(nrNotPresent + " of the query SNPs were not present in the data");

        TextFile outnotinref = new TextFile(outdir + "NotInReference.txt", TextFile.W);
        for (String s : snpsNotInReference) {
            outnotinref.writeln(s);
        }
        outnotinref.close();

        int sum = 0;
        for (int i = 0; i < distBins.size(); i++) {
            ArrayList<String> distbins = distBins.get(i);
            if (distbins != null) {
                System.out.println(i + "\t" + (i * distanceGranularity) + "\t" + distbins.size() + "\t" + queryDistbins[i]);
                sum += distbins.size();
            } else {
                System.out.println(i + "\t" + (i * distanceGranularity) + "\t" + 0 + "\t" + queryDistbins[i]);
            }
        }

        System.out.println("Total: " + sum);
        System.out.println("");

        sum = 0;
        for (int i = 0; i < mafBins.size(); i++) {
            ArrayList<String> distbins = mafBins.get(i);
            if (distbins != null) {
                System.out.println(i + "\t" + (i * mafGranularity) + "\t" + distbins.size() + "\t" + queryMafbins[i]);
            } else {
                System.out.println(i + "\t" + (i * mafGranularity) + "\t" + 0 + "\t" + queryMafbins[i]);
            }
        }

        System.out.println("Total: " + sum);
        System.out.println("");

        DetermineLD ldcalc = new DetermineLD();

        outdir = Gpio.formatAsDirectory(outdir);
        Gpio.createDir(outdir);

        TextFile q = new TextFile(outdir + "querySNPs.txt", TextFile.W);
        for (String snp : querySNPs) {
            q.writeln(snp);
        }
        q.close();

        // now picking randomsets //
        for (int set = 0; set < nrOfSetsToPick; set++) {

            System.out.println();
            HashSet<String> finalPickedSNPs = new HashSet<String>();
            // for each query SNP, find a similar one
            int sctr = 0;
            ProgressBar pb = new ProgressBar(querySNPs.size(), "Picking set " + set);
            for (String snp : querySNPs) {
                Integer qSNPMafBin = querySNPmafBins.get(snp);
                Integer qSNPDistanceBin = querySNPdistBins.get(snp);

                if (qSNPDistanceBin != null && qSNPMafBin != null) {

                    // select similar MAF:
                    ArrayList<String> snpsWithSimilarMAF = mafBins.get(qSNPMafBin);

                    // and then select the ones that are also in 
                    ArrayList<String> snpsWithSimilarDistance = distBins.get(qSNPDistanceBin);
                    HashSet<String> tmpDistSNPs = new HashSet<String>();
                    tmpDistSNPs.addAll(snpsWithSimilarDistance);
                    // intersect
                    ArrayList<String> intersect = new ArrayList<String>();
                    for (String s : snpsWithSimilarMAF) {
                        if (tmpDistSNPs.contains(s)) {
                            intersect.add(s);
                        }
                    }
                    if (intersect.size() > 0) {
                        int iterations = 0;
                        boolean stoprunning = false;
                        boolean snpPicked = false;

                        HashSet<Integer> indexesPicked = new HashSet<Integer>();
                        String finalSelection = null;
                        while (!stoprunning) {
                            if (indexesPicked.size() == intersect.size()) {
                                stoprunning = true;
                            }

                            if (!stoprunning) {
                                int randomIndex = (int) Math.floor(Math.random() * intersect.size());
                                String selectedSNP = intersect.get(randomIndex);

                                indexesPicked.add(randomIndex);

                                if (finalPickedSNPs.isEmpty()) {
                                    // when this is the first SNP to be picked
                                    finalPickedSNPs.add(selectedSNP);
                                    finalSelection = selectedSNP;
                                    snpPicked = true;
                                    stoprunning = true;
                                } else {

                                    // check if the SNP is not already in the list of picked SNPs. If so, this is not a proper pick
                                    if (!finalPickedSNPs.contains(selectedSNP)) {

                                        // check whether the selectedSNP is not in LD with any of the other added snps
                                        Integer id = ds.getSnpToSNPId().get(selectedSNP);
                                        byte chr = ds.getChr(id);
                                        SNP snpObj = ds.getSNPObject(id);
                                        loader.loadGenotypes(snpObj);
                                        boolean isInLD = false;

                                        // check all other already picked SNPs
                                        for (String s : finalPickedSNPs) {
                                            Integer id2 = ds.getSnpToSNPId().get(s);
                                            byte chr2 = ds.getChr(id);

                                            // check whether they are on the same chromosome...
                                            if (chr == chr2) {
                                                SNP snpObj2 = ds.getSNPObject(id2);
                                                loader.loadGenotypes(snpObj2);
                                                double r2 = ldcalc.getRSquared(snpObj, snpObj2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);

                                                // check whether there is linkage
                                                if (r2 > r2threshold) {
                                                    isInLD = true;
                                                }
                                                snpObj2.clearGenotypes();
                                            }
                                        }

                                        snpObj.clearGenotypes();
                                        if (!isInLD) {
                                            snpPicked = true;
                                            finalPickedSNPs.add(selectedSNP); // LD check 0.2
                                            finalSelection = selectedSNP;
                                            stoprunning = true;
                                        }
                                    }
                                }
                            }
                        }
                        if (!snpPicked) {
                            System.out.println(set + "\t" + sctr + "/" + querySNPs.size() + "\t" + intersect.size() + "\t" + snpPicked + "\t" + snp + "\t" + qSNPMafBin + "\t" + qSNPDistanceBin + "\t" + finalSelection);
                        }
                    }
                }

                pb.iterate();
                sctr++;

            }
            pb.close();
            System.out.println("SET: " + set + "\tPicked: " + finalPickedSNPs.size());


            TextFile out = new TextFile(outdir + "Set-" + set + ".txt", TextFile.W);

            int[] pickedDistbins = new int[distBins.size()];
            int[] pickedMafbins = new int[mafBins.size()];
            for (String s : finalPickedSNPs) {
                int mafb = getMAFBin(ds, loader, s, mafGranularity);
                int distb = getDistanceBin(ds, s, chromosomes, distanceGranularity);
                pickedDistbins[distb]++;
                pickedMafbins[mafb]++;
                out.writeln(s);
            }
            out.close();

            System.out.println("Distributions: ");
            System.out.println("");
            TextFile outDist = new TextFile(outdir + "Set-" + set + "-Dist.txt", TextFile.W);
            sum = 0;
            for (int i = 0; i < distBins.size(); i++) {
                ArrayList<String> distbins = distBins.get(i);
                if (distbins != null) {
                    System.out.println(i + "\t" + (i * distanceGranularity) + "\t" + distbins.size() + "\t" + queryDistbins[i] + "\t" + pickedDistbins[i]);
                    outDist.writeln(i + "\t" + (i * distanceGranularity) + "\t" + distbins.size() + "\t" + queryDistbins[i] + "\t" + pickedDistbins[i]);
                    sum += distbins.size();
                } else {
                    System.out.println(i + "\t" + (i * distanceGranularity) + "\t" + 0 + "\t" + queryDistbins[i] + "\t" + pickedDistbins[i]);
                    outDist.writeln(i + "\t" + (i * distanceGranularity) + "\t" + 0 + "\t" + queryDistbins[i] + "\t" + pickedDistbins[i]);
                }
            }

            outDist.writeln();
            System.out.println("Total: " + finalPickedSNPs.size());
            System.out.println("");

            sum = 0;
            for (int i = 0; i < mafBins.size(); i++) {
                ArrayList<String> distbins = mafBins.get(i);
                if (distbins != null) {
                    System.out.println(i + "\t" + (i * mafGranularity) + "\t" + distbins.size() + "\t" + queryMafbins[i] + "\t" + pickedMafbins[i]);
                    outDist.writeln(i + "\t" + (i * mafGranularity) + "\t" + distbins.size() + "\t" + queryMafbins[i] + "\t" + pickedMafbins[i]);
                } else {
                    System.out.println(i + "\t" + (i * mafGranularity) + "\t" + 0 + "\t" + queryMafbins[i] + "\t" + pickedMafbins[i]);
                    outDist.writeln(i + "\t" + (i * mafGranularity) + "\t" + 0 + "\t" + queryMafbins[i] + "\t" + pickedMafbins[i]);
                }
            }
            System.out.println("");
            System.out.println("---- ");
            System.out.println("");
            outDist.close();
        }

        loader.close();
    }

    public static Integer getMAFBin(TriTyperGenotypeData ds, SNPLoader loader, String snp, double mafGranularity) throws IOException {
        Integer snpId = ds.getSnpToSNPId().get(snp);
        if (snpId != null) {
            SNP snpobj = ds.getSNPObject(snpId);
            loader.loadGenotypes(snpobj);
            double maf = snpobj.getMAF();

            int bin = (int) Math.round(maf / mafGranularity); // 0.05 -> 1, 0.1 --> 2


            snpobj.clearGenotypes();
            return bin;

        }
        return null;
    }

    public static Integer getDistanceBin(TriTyperGenotypeData ds, String snp, HashMap<String, Chromosome> chromosomes, double distanceGranularity) throws IOException {
        Integer snpId = ds.getSnpToSNPId().get(snp);
        if (snpId != null) {
            String chr = "" + ds.getChr(snpId);
            int snppos = ds.getChrPos(snpId);

            Chromosome c = chromosomes.get(chr);
            if (c == null) {
                System.out.println(snp + "\t" + chr + "\t is NULL!");
            }
            ArrayList<Gene> genes = c.getGenes();
            int minDist = Integer.MAX_VALUE;
            for (Gene g : genes) {
                int geneStart = g.getStart();
                int distance = Math.abs(geneStart - snppos);
                if (distance < minDist) {
                    minDist = distance;
                }

            }
            int bin = (int) Math.round(minDist / distanceGranularity); // 0.05 -> 1, 0.1 --> 2
            return bin;
        }
        return null;

    }

    public static void selectSNPsWithSimilarMaf(String querySNPFile, String confineToTheseSNPs, String reference, String outfile) throws IOException {
        TextFile queryFile = new TextFile(querySNPFile, TextFile.R);
        HashSet<String> querySNPs = new HashSet<String>();
        querySNPs.addAll(queryFile.readAsArrayList());
        queryFile.close();
        System.out.println(querySNPs.size() + " SNPs loaded as query");

        TextFile confinefile = new TextFile(confineToTheseSNPs, TextFile.R);
        HashSet<String> confineSNPs = new HashSet<String>();
        confineSNPs.addAll(confinefile.readAsArrayList());
        confinefile.close();
        System.out.println(confineSNPs.size() + " SNPs loaded as confinement");

        String[] querySNPArr = querySNPs.toArray(new String[0]);

        // load reference dataset
        TriTyperGenotypeData ds = new TriTyperGenotypeData(reference);
        SNPLoader loader = ds.createSNPLoader();
        String[] snps = ds.getSNPs();

        TextFile out = new TextFile(outfile, TextFile.W);
        int[] mafbins = new int[querySNPs.size()];
        HashSet<String> selectedSNPs = new HashSet<String>();
        for (int i = 0; i < mafbins.length; i++) {
            System.out.println(querySNPArr[i]);
            int mafbin = getMAFBin(ds, loader, querySNPArr[i], 0.001);
            System.out.println(mafbin);
            boolean snpselected = false;
            HashSet<String> snpsToSelectFrom = new HashSet<String>();
            for (int s = 0; s < snps.length; s++) {
                String snp = snps[s];

                if (!querySNPs.contains(snp) && !selectedSNPs.contains(snp)) {
                    if (confineSNPs.contains(snp)) {
                        int mafbin2 = getMAFBin(ds, loader, snp, 0.001);
//                    System.out.println(mafbin + "\t" + snp + "\t"+mafbin2);
                        if (mafbin2 == mafbin) {
                            selectedSNPs.add(snp);
                            System.out.println(querySNPArr[i] + "\t" + snp);
                            out.writeln(snp);
                            snpselected = true;
                        }
                    }
                }
            }
        }

        out.close();
        loader.close();
    }

    static Integer getMAFBin(TriTyperGenotypeData ds, SNPLoader loader, String snp, double mafGranularity, boolean b) throws IOException {
        Integer snpId = ds.getSnpToSNPId().get(snp);
        if (snpId != null) {
            SNP snpobj = ds.getSNPObject(snpId);
            loader.loadGenotypes(snpobj);
            double maf = snpobj.getMAF();

            int bin = (int) Math.round(maf / mafGranularity); // 0.05 -> 1, 0.1 --> 2


            snpobj.clearGenotypes();
            return bin;

        } else {
            System.out.println("SNP not present in reference... ");
        }
        return null;
    }
}
