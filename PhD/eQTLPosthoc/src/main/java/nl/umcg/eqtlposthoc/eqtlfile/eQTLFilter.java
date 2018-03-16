/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import eqtlmappingpipeline.metaqtl3.FDR;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLFilter {

    public static void main(String[] args) {
        try {

//            int[] arrr = new int[]{25,30,35,40};
//            for(int q: arrr){
//                eQTLFilter.filterEQTLFileForSNPProbeCombos(
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\bcellmeta\\"+q+"PC\\",
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\bcellmetaFiltered\\"+q+"PC\\",
//                    10,
//                    "d:\\SkyDrive\\latesteQTLs\\transFDR0.05.txt.gz"
//                    ,1,4);  
//            }
//            
//            arrr = new int[]{20,30,35,40,50};
//            for(int q: arrr){
//                eQTLFilter.filterEQTLFileForSNPProbeCombos(
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometa\\"+q+"PC\\",
//                    "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometaFiltered\\"+q+"PC\\",
//                    10,
//                    "d:\\SkyDrive\\latesteQTLs\\transFDR0.05.txt.gz"
//                    ,1,4);  
//            }
//            
////            eQTLFilter.filterEQTLFileForProbesMappingInGenes(
////                    "/Volumes/iSnackHD/Skydrive/MetaAnalysisAnnotationFiles/2012-08-08-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation.txt",
////                    "/Volumes/iSnackHD/SkyDrive/latesteQTLs/cisProbeLevelFDR0.05.txt",
////                    "/Volumes/iSnackHD/SkyDrive/latesteQTLs/cisProbeLevelFDR0.05.txt-WithDistanceToGeneEnd.txt");

//            String probeList = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt";
//            String fileIn = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/MetaOut/eQTLs.txt";
//            String fileOut = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-02-KLF14TransEQTLs/MetaOutFiltered/eQTLs.txt";

//            eQTLFilter.filterForProbes(probeList, fileIn, fileOut);


//            String fileIn = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR.txt.gz";
//            String fileOut = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs.txt.gz";
//            double threshold = 1.3141052405861965E-4;
//            eQTLFilter.filterForPValueThreshold(threshold, fileIn, fileOut);
//            
            String fileIn = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs.txt.gz";
            String fileOut = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLsFDR-AllProbeLevelQTLs-OnlyGWASSignigicantSNPs.txt.gz";
            String snplist = "/Volumes/iSnackHD/SkyDrive2/SkyDrive/Rebuttal2/transeqtls/GenomeWideAccordingToLude.txt";
            eQTLFilter.filterForSNPs(snplist, fileIn, fileOut);

//            String probeList = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-03-Trans40PCsCisFXNotRegressedOut/ComparisonTo40PCCisFxRegressedOut/NonOverlappingProbes.txt";
//            String fileIn = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/Results/2012-01-11-MetaAnalysis/cis/2012-02-22-EGCUT+SHIP_TREND+Groningen+Rotterdam+DILGOM+INCHIANTI+HVH-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/eQTLProbesFDR0.05.txt";
//            String fileOut = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/2013-03-Trans40PCsCisFXNotRegressedOut/ComparisonTo40PCCisFxRegressedOut/CisFXForNonOverlappingProbesRun2.txt";
//            eQTLFilter.filterForProbes(probeList, fileIn, fileOut);
            
//            String indir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2CisEffectsRemoved/";
//            String outdir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/Replication/KatharinaHeim/Kora2CisEffectsRemovedFilteredForFDR0.05EQTLs/";
//            int nrPerm = 10;
//            String filterfile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/2012-07-27-Groningen+EGCUT+Rotterdam+DILGOM+SHIP+InChianti+HVH-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved-CisEffectsRegressedOut/PostQC/eQTLsFDR0.05-HT12v3.txt.gz";
//            eQTLFilter.filterEQTLFileForSNPProbeCombos(indir, outdir, nrPerm, filterfile, 1, 4);
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void filterForProbes(String probeList, String fileIn, String fileOut) throws IOException {
        HashSet<String> allowedProbes = new HashSet<String>();

        TextFile tf = new TextFile(probeList, TextFile.R);
        allowedProbes.addAll(tf.readAsArrayList());
        tf.close();

        TextFile fout = new TextFile(fileOut, TextFile.W);
        TextFile fin = new TextFile(fileIn, TextFile.R);
        fout.writeln(fin.readLine());
        String[] elems = fin.readLineElems(TextFile.tab);
        while (elems != null) {
            if (allowedProbes.contains(elems[4])) {
                fout.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = fin.readLineElems(TextFile.tab);
        }
        fin.close();
        fout.close();
    }

    private static void filterForPValueThreshold(double threshold, String fileIn, String fileOut) throws IOException {
        TextFile tf = new TextFile(fileIn, TextFile.R);
        TextFile tfout = new TextFile(fileOut, TextFile.W);
        tfout.writeln(tf.readLine());
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            double p = Double.parseDouble(elems[0]);
            if (p < threshold) {
                tfout.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        tfout.close();

    }

    private static void filterForSNPs(String snplist, String fileIn, String fileOut) throws IOException {
        HashSet<String> allowedSNPs = new HashSet<String>();

        TextFile tf = new TextFile(snplist, TextFile.R);
        allowedSNPs.addAll(tf.readAsArrayList());
        tf.close();

        TextFile fout = new TextFile(fileOut, TextFile.W);
        TextFile fin = new TextFile(fileIn, TextFile.R);
        fout.writeln(fin.readLine());
        String[] elems = fin.readLineElems(TextFile.tab);
        while (elems != null) {
            if (allowedSNPs.contains(elems[1])) {
                fout.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = fin.readLineElems(TextFile.tab);
        }
        fin.close();
        fout.close();
    }

    public void run(String queryfile, String eqtls, String outfile) throws IOException {
        HashSet<String> querySet = new HashSet<String>();
        TextFile query = new TextFile(queryfile, TextFile.R);
        querySet.addAll(query.readAsArrayList());
        query.close();

        System.out.println("Query consists of " + querySet.size() + " elements");
        TextFile eqtlfile = new TextFile(eqtls, TextFile.R);
        TextFile out = new TextFile(outfile, TextFile.W);

        String header = eqtlfile.readLine();
        String[] elems = eqtlfile.readLineElems(TextFile.tab);

        out.writeln(header);

        int ctr = 0;
        while (elems != null) {
            for (int i = 0; i < elems.length; i++) {
                if (querySet.contains(elems[i])) {
                    out.writeln(Strings.concat(elems, Strings.tab));
                    ctr++;
                }
            }
            elems = eqtlfile.readLineElems(TextFile.tab);
        }

        System.out.println(ctr + " eQTLs filtered for query.");
        out.close();
        eqtlfile.close();
    }

    public static void filterEQTLFileForSNPProbeCombos(String indir, String outdir, int nrPerm, String snpprobefile, int snpcol, int probecol) throws IOException {
        outdir += "/";
        indir += "/";
        Gpio.createDir(outdir);
        TextFile tf = new TextFile(snpprobefile, TextFile.R);
        HashSet<Pair<String, String>> snpProbeCombos = tf.readAsPairs(snpcol, probecol);
        tf.close();

        for (int perm = 0; perm < nrPerm + 1; perm++) {
            String infile = "";
            String outfile = "";
            if (perm == 0) {
                infile = indir + "/eQTLs.txt.gz";
                outfile = outdir + "/eQTLs.txt.gz";
            } else {
                infile = indir + "/PermutedEQTLsPermutationRound" + perm + ".txt.gz";
                outfile = outdir + "/PermutedEQTLsPermutationRound" + perm + ".txt.gz";

            }

            TextFile fileIn = new TextFile(infile, TextFile.R);
            TextFile fileOut = new TextFile(outfile, TextFile.W);
            fileOut.writeln(fileIn.readLine());
            String[] elems = fileIn.readLineElems(TextFile.tab);
            while (elems != null) {
                String snp = elems[1];
                String probe = elems[4];
                Pair<String, String> eqtl = new Pair<String, String>(snp, probe);

                if (snpProbeCombos.contains(eqtl)) {
                    fileOut.writeln(Strings.concat(elems, Strings.tab));
                }
                elems = fileIn.readLineElems(TextFile.tab);
            }
            fileIn.close();
            fileOut.close();
        }

        // now perform FDR
        String eQTLFile = outdir + "/eQTLs.txt.gz";
        TextFile t = new TextFile(eQTLFile, TextFile.R);
        int nreQTLs = t.countLines() - 1;
        t.close();

        if (nrPerm > 0) {
            FDR.calculateFDR(outdir, nrPerm, nreQTLs, 0.05, true, null, null);
        }
    }

    public static void filterEQTLFileForProbesMappingInGenes(String ensemblExonAnnotationFile, String eqtlFile, String outputfile) throws IOException {

        System.out.println("Writing to: " + outputfile);
        TextFile tf = new TextFile(ensemblExonAnnotationFile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        elems = tf.readLineElems(TextFile.tab);
        HashMap<String, Integer> probesThatMapToAGene = new HashMap<String, Integer>();
        HashMap<String, Integer> probesThatMapToAGeneEnd = new HashMap<String, Integer>();
        HashMap<String, Boolean> geneIsOnPositiveStrand = new HashMap<String, Boolean>();
        while (elems != null) {

            if (elems.length > 7 && !elems[4].equals("-")) {
                String probe = elems[0];
                String pos = elems[6];
                String strand = elems[7];
                boolean positiveStrand = true;
                if (strand.equals("-1")) {
                    positiveStrand = false;
                }
                String[] posElems = pos.split("-");
                Integer chrPos = Integer.parseInt(posElems[0]);
                Integer chrPosEnd = Integer.parseInt(posElems[posElems.length - 1]);

                System.out.println(chrPos + " - " + chrPosEnd + "\t" + (chrPosEnd - chrPos));
                probesThatMapToAGene.put(probe, chrPos);
                probesThatMapToAGeneEnd.put(probe, chrPosEnd);
                geneIsOnPositiveStrand.put(probe, positiveStrand);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile eqtl = new TextFile(eqtlFile, TextFile.R);
        TextFile tfout = new TextFile(outputfile, TextFile.W);
        elems = eqtl.readLineElems(TextFile.tab);
        elems[11] = "AbsMetaZ";


        elems[12] = "DistanceToGeneStart";
        elems[13] = "DistanceToGeneEnd";

        elems[14] = "GeneLength";
        elems[15] = "AbsDistanceToGeneStart";
        elems[16] = "SNP-ProbeDistance";
        elems[17] = "AbsSNP-ProbeDistance";

        elems[18] = "GeneStart";
        elems[19] = "GeneEnd";


        tfout.writeln(Strings.concat(elems, Strings.tab));
        elems = eqtl.readLineElems(TextFile.tab);
        double[] first8k = new double[7000];
        ArrayList<Double> lasteQTLs = new ArrayList<Double>();
        int ctr = 0;
        while (elems != null) {

            String probe = elems[4];
            Integer probePos = Integer.parseInt(elems[6]);
            Integer snpPos = Integer.parseInt(elems[3]);
            Integer genestart = probesThatMapToAGene.get(probe);
            int distanceToGeneStart = 0;
            int distanceToGeneEnd = 0;
            int snpProbeDistance = 0;
            snpProbeDistance = snpPos - probePos;

            Integer geneEnd = probesThatMapToAGeneEnd.get(probe);
            if (genestart != null) {
                Boolean strand = geneIsOnPositiveStrand.get(probe);
                if (!strand) {
//                    distanceToGeneStart = snpPos - geneEnd;
//                    distanceToGeneEnd = snpPos - genestart;

                    distanceToGeneStart = Integer.MAX_VALUE;
                    distanceToGeneEnd = Integer.MAX_VALUE;

                } else {
                    distanceToGeneStart = snpPos - genestart;
                    distanceToGeneEnd = snpPos - geneEnd;
                }


                int difference = distanceToGeneStart - distanceToGeneEnd;
                System.out.println(difference);
                Double z = Double.parseDouble(elems[10]);

                elems[11] = "" + Math.abs(z);

                elems[12] = "" + distanceToGeneStart;
                elems[13] = "" + distanceToGeneEnd;

                elems[14] = "" + (geneEnd - genestart);
                elems[15] = "" + Math.abs(distanceToGeneStart);
                elems[18] = "" + genestart;
                elems[19] = "" + geneEnd;

            } else {
                elems[11] = "-";
                elems[12] = "-";
                elems[13] = "-";
                elems[14] = "-";
                elems[15] = "-";
                elems[18] = "-";
                elems[19] = "-";
            }

            elems[16] = "" + snpProbeDistance;
            elems[17] = "" + Math.abs(snpProbeDistance);
            tfout.writeln(Strings.concat(elems, Strings.tab));
            if (ctr < first8k.length) {
                first8k[ctr] = snpProbeDistance;
                ctr++;
            } else {
                lasteQTLs.add((double) snpProbeDistance);
            }
            ctr++;
            elems = eqtl.readLineElems(TextFile.tab);
        }

        Double[] dataLastEQTL = lasteQTLs.toArray(new Double[0]);
        double[] lasteQTLsData = new double[dataLastEQTL.length];

        for (int d = 0; d < lasteQTLsData.length; d++) {
            lasteQTLsData[d] = dataLastEQTL[d];
        }

        WilcoxonMannWhitney m = new WilcoxonMannWhitney();
        double p = m.returnWilcoxonMannWhitneyPValue(first8k, lasteQTLsData);

        System.out.println(p);
        System.out.println(m.getAUC());

        tfout.close();
        eqtl.close();

    }
}
