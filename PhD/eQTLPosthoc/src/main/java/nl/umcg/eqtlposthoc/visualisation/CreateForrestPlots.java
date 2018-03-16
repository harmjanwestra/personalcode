/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.visualisation;

import eqtlmappingpipeline.metaqtl3.containers.WorkPackage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.graphics.ForestPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harm-jan
 */
public class CreateForrestPlots {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
//        String[] snpsToPlot = new String[]{"rs174548"};
        String[] snpsToPlot = new String[]{"rs653178", "rs174550", "rs4917014", "rs4788084", "rs3184504"};
//        String[] snpsToPlot = new String[]{"rs174550",
//            "rs174546",
//            "rs174547",
//            "rs174535",
//            "rs174536",
//            "rs102275",
//            "rs174583",
//            "rs174548",
//            "rs1535",
//            "rs174574",
//            "rs4246215",
//            "rs174538",
//            "rs174448",
//            "rs1000778",
//            "rs174570"
//        };
//      double[] significanceThresholds = new double[]{5.021869225607994, 4.3212320038155205}; // trans
        String[] datasetnames = new String[]{"EGCUT", "SHIP-Trend", "Fehrmann HT12", "Fehrmann H8v2", "Rotterdam Study", "DILGOM", "InCHIANTI", "HVH HT12v3", "HVH HT12v4", "Meta-analysis", "B-Cells", "Monocytes", "Kora", "BSGS"};
        int[] datasetweights = new int[]{891, 963, 1240, 229, 762, 509, 611, 43, 63, 5311, 10, 10, 10, 10};
        String xAxisName = "ZScore";

        // D:\\SkyDrive\latesteQTLs\\Replication\\FairFaxAllProbeNorm\\FairFaxBCellFiltered\\40_PC_all\\
        // D:\\SkyDrive\latesteQTLs\\Replication\\FairFaxAllProbeNorm\\FairFaxMonoFiltered\\40_PC_all\\
        // D:\\SkyDrive\\latesteQTLs\\Replication\\Kora2Filtered\\
        String[] replicationFiles = new String[]{"D:\\SkyDrive\\latesteQTLs\\Replication\\FairFaxAllProbeNorm\\FairFaxBCellFiltered\\40_PC_all\\eQTLsFDR.txt.gz",
            "D:\\SkyDrive\\latesteQTLs\\Replication\\FairFaxAllProbeNorm\\FairFaxMonoFiltered\\40_PC_all\\eQTLsFDR.txt.gz",
            "D:\\SkyDrive\\latesteQTLs\\Replication\\Kora2Filtered\\eQTLsFDR.txt.gz",
            "D:\\SkyDrive\\latesteQTLs\\Replication\\BSGS\\eQTLsFDR.txt.gz"};
        int metaRow = 9;
        String output = "PDF";
        String outdir = "d:\\work\\original\\";
        try {
            Gpio.createDir(outdir);
        } catch (Exception e) {
        }
        boolean useFlippedFile = true;
        double fdr = 0.05;
        for (int r = 0; r < 2; r++) {
            String suffix = "trans";
            double[] significanceThresholds = new double[]{5.021869225607994, 4.3212320038155205};
            double pvalThreshold = 1;
            if (r > 0) {
                suffix = "cis";
                pvalThreshold = 1.3141052405861965E-4;
                if (fdr == 0.5) {
                    pvalThreshold = 0.00352058461224396;
                }
                // 
                significanceThresholds = new double[]{3.8238746857983563, 2.918359859509293};
            }

            try {
                for (String snp : snpsToPlot) {
                    HashMap<String, String> probeToGene = new HashMap<String, String>();
                    String inFile = "";
                    if (useFlippedFile) {
                        inFile = "D:\\Skydrive\\latesteQTLs\\" + suffix + "FDR" + fdr + "-flipped.txt.gz";
                    } else {
                        inFile = "D:\\Skydrive\\latesteQTLs\\" + suffix + "FDR" + fdr + ".txt.gz";
                    }
                    TextFile tf = new TextFile(inFile, TextFile.R);
                    String[] data = tf.readLineElems(TextFile.tab);
                    data = tf.readLineElems(TextFile.tab);
                    ArrayList<SortableSNP> selectedGenes = new ArrayList<SortableSNP>();
                    int nr = 0;
                    HashMap<String, Double[]> effectSizes = new HashMap<String, Double[]>();
                    String[] rowNames = null;

                    HashMap<String, Integer> geneToId = new HashMap<String, Integer>();
                    Double max = Double.MIN_VALUE;
                    Double min = Double.MAX_VALUE;

                    String alleles = null;
                    String assessed = null;

                    while (data != null) {
                        double pval = Double.parseDouble(data[0]);
                        if (pval < pvalThreshold) {
                            if (data[1].equals(snp)) {

                                String gene = data[eQTLTextFile.HUGO];
                                if (!gene.equals("-")) {
                                    if (alleles == null) {
                                        alleles = data[eQTLTextFile.ASESSEDALLELE - 1];
                                        assessed = data[eQTLTextFile.ASESSEDALLELE];
                                    }
                                    String probeChr = data[eQTLTextFile.PROBECHR];
                                    int chrPos = Integer.parseInt(data[eQTLTextFile.PROBELOC]);
                                    byte probeChrB = ChrAnnotation.parseChr(probeChr);
                                    int replicates = 2;
                                    String origGene = gene;
                                    while (effectSizes.containsKey(gene)) {
                                        gene = origGene + "-" + replicates;
                                        replicates++;
                                    }
                                    SortableSNP probe = new SortableSNP(gene, nr, probeChrB, chrPos, SortableSNP.SORTBY.EFFECT);
                                    String dsZScores = data[eQTLTextFile.DATASETZSCORE];
                                    String dsSizes = data[eQTLTextFile.DATASETSIZE];
                                    String dsNames = data[eQTLTextFile.DATASETNAMES];



                                    String[] dsZElems = dsZScores.split(",");
                                    if (dsZElems.length == 1) {
                                        dsZElems = dsZScores.split(";");
                                    }

                                    String[] dsSElems = dsSizes.split(",");
                                    if (dsSElems.length == 1) {
                                        dsSElems = dsSizes.split(";");
                                    }

                                    String[] dsNElems = dsNames.split(",");
                                    if (dsNElems.length == 1) {
                                        dsNElems = dsNames.split(";");
                                    }


                                    if (rowNames == null) {
                                        rowNames = new String[datasetnames.length];
                                    }

                                    Double[] zScores = new Double[datasetnames.length];
                                    Integer[] weights = new Integer[datasetnames.length];
                                    int totalWeight = 0;
                                    for (int d = 0; d < dsZElems.length; d++) {
                                        if (dsZElems[d].equals("-")) {
                                            zScores[d] = null;
                                        } else {
                                            if (rowNames[d] == null) {
                                                rowNames[d] = dsNElems[d];
                                            }
//                                    if(dsSElems[d].equals("-")){
//                                        System.out.println("ERROR IN FILE? "+dsZElems[d]);
//                                        System.out.println(Strings.concat(data, Strings.tab));
//                                        System.exit(0);
//                                    }
                                            if (datasetweights == null) {
                                                weights[d] = Integer.parseInt(dsSElems[d]);
                                                totalWeight += weights[d];
                                            }

                                            double z = Double.parseDouble(dsZElems[d]);
                                            if (z > max) {
                                                max = z;
                                            }
                                            if (z < min) {
                                                min = z;
                                            }
                                            zScores[d] = z;
                                            System.out.println(d + "\t" + datasetnames[d] + "\t" + datasetweights[d] + "\t" + z);
                                        }
                                    }


                                    probeToGene.put(data[4], probe.name);
                                    double metaZ = Double.parseDouble(data[eQTLTextFile.METAZ]);
                                    SortableSNP probe2 = new SortableSNP(probe.name, probe.id, probe.chr, probe.chrpos, metaZ, SortableSNP.SORTBY.EFFECT);
                                    selectedGenes.add(probe2);
                                    if (metaZ > max) {
                                        max = metaZ;
                                    }
                                    if (metaZ < min) {
                                        min = metaZ;
                                    }
                                    rowNames[metaRow] = "Meta-Analysis";
                                    if (datasetweights != null) {
                                        for (int w : datasetweights) {
                                            totalWeight += w;
                                        }
                                    } else {
                                        weights[metaRow] = totalWeight;
                                    }
                                    zScores[metaRow] = metaZ;
                                    System.out.println("MetaZ: " + metaZ);
                                    effectSizes.put(gene, zScores);
                                    nr++;
                                }
                            }
                        }
                        data = tf.readLineElems(TextFile.tab);
                    }

                    if (selectedGenes == null || selectedGenes.isEmpty()) {
                        System.out.println("ERROR: no genes found matching to SNP: " + snp + "\t in file: " + tf.getFileName());
                    } else {


                        // now load the replication data..
                        for (int rep = 0; rep < replicationFiles.length; rep++) {
                            String fileName = replicationFiles[rep];


                            // iterate through the replication file, and search for SNP-probe combo's that are in the meta-analysis

                            // from probeToGene get the gene name
                            // get the effect sizes, and add replication effects (what about allele directions?)
                            TextFile tfrep = new TextFile(fileName, TextFile.R);
                            String[] repdata = tfrep.readLineElems(TextFile.tab);
                            // repdata = tfrep.readLineElems(TextFile.tab);
                            while (repdata != null) {
                                if (snp.equals(repdata[1])) {
                                    String gene = probeToGene.get(repdata[4]);
                                    if (gene != null) {

                                        String alleles2 = repdata[eQTLTextFile.ASESSEDALLELE - 1];
                                        String assessed2 = repdata[eQTLTextFile.ASESSEDALLELE];

                                        Double[] effects = effectSizes.get(gene);
                                        Double z = Double.parseDouble(repdata[eQTLTextFile.METAZ]);
                                        int samplesize = 0;

                                        try {
                                            samplesize = Integer.parseInt(repdata[eQTLTextFile.DATASETSIZE]);
                                        } catch (NumberFormatException e) {
                                            samplesize = 892;
                                        }
                                        datasetweights[metaRow + rep + 1] = samplesize;
                                        if (z > max) {
                                            max = z;
                                        }
                                        if (z < min) {
                                            min = z;
                                        }
                                        boolean flip = detmermineAlleleFlips(alleles, assessed, alleles2, assessed2);
                                        // check for flipped alleles
                                        if (flip) {
                                            z *= -1;
                                            System.out.println("Flipping effect: " + alleles + "\t" + assessed + "\t" + alleles2 + "\t" + assessed2);
                                        }
                                        double fdrforcombo = Double.parseDouble(repdata[repdata.length - 1]);
                                        if (fdrforcombo <= fdr) {
                                            effects[metaRow + rep + 1] = z;
                                        }
                                    }
                                }
                                repdata = tfrep.readLineElems(TextFile.tab);
                            }
                            tfrep.close();



                        }


                        Collections.sort(selectedGenes);
                        Double[][] plotvals = new Double[selectedGenes.size()][rowNames.length];
                        System.out.println("Final matrix is: " + selectedGenes.size() + " x " + rowNames.length);
                        String[] geneNames = new String[selectedGenes.size()];
                        byte[] chrs = new byte[geneNames.length];
                        int[] chrspos = new int[geneNames.length];
                        for (int i = 0; i < plotvals.length; i++) {
                            SortableSNP gene = selectedGenes.get(i);
                            String name = gene.name;
                            plotvals[i] = effectSizes.get(name);
                            geneNames[i] = selectedGenes.get(i).name;
                            chrs[i] = selectedGenes.get(i).chr;
                            chrspos[i] = selectedGenes.get(i).chrpos;
                        }
                        String[] yAxisNames = rowNames;
                        String filename = outdir + snp + "-" + suffix + "-FDR" + fdr + ".pdf";
                        if (output.equals("PNG")) {
                            filename = outdir + snp + "-" + suffix + "-FDR" + fdr + ".png";
                        }
                        ForestPlot p = new ForestPlot();

                        System.out.println("Min: " + min + "\tmax: " + max);

                        // equalize min and max..
                        max = Math.ceil(Math.abs(max));
                        System.out.println("Min: " + min + "\tmax: " + max);
                        double absmin = Math.ceil(Math.abs(min));
                        System.out.println("AbsMin: " + absmin + "\tmax: " + max);
                        if (max < absmin) {
                            max = absmin;
                        }
                        // we've found the boundraries, make them modulo 5

                        double remainder = max % 5;
                        double toAdd = 5 - remainder;
                        max += toAdd;
                        min = -max;

                        System.out.println("Min: " + min + "\tmax: " + max);
                        p.setGeneNames(geneNames);


                        if (output.equals("PNG")) {
                            p.drawMultiForrestPlot(xAxisName, datasetnames, plotvals, filename, ForestPlot.Output.PNG, significanceThresholds, min, max, datasetweights, chrs, chrspos, metaRow);
                        } else {
                            p.drawMultiForrestPlot(xAxisName, datasetnames, plotvals, filename, ForestPlot.Output.PDF, significanceThresholds, min, max, datasetweights, chrs, chrspos, metaRow);
                        }

                        tf.close();
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }

    public static boolean detmermineAlleleFlips(String allelesDs1, String alleleAssessedDs1, String allelesDs2, String alleleAssessedDs2) throws IOException {

        byte[] firstDatasetAlleles = new byte[2];
        byte[] secondDatasetAlleles = new byte[2];

        String[] alleles1 = allelesDs1.split("/");
        String[] alleles2 = allelesDs1.split("/");
        for (int i = 0; i < 2; i++) {
            firstDatasetAlleles[i] = BaseAnnot.toByte(alleles1[i]);
            secondDatasetAlleles[i] = BaseAnnot.toByte(alleles2[i]);
        }
        byte minor1 = BaseAnnot.toByte(alleleAssessedDs1);
        byte minor2 = BaseAnnot.toByte(alleleAssessedDs2);
        int nrAllelesIdentical = 0;
        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 2; b++) {
                if (firstDatasetAlleles[a] == secondDatasetAlleles[b]) {
                    nrAllelesIdentical++;
                }
            }
        }

        if (nrAllelesIdentical != 2) {
            //Alleles are different, take complimentary:
            secondDatasetAlleles = BaseAnnot.convertToComplementaryAlleles(secondDatasetAlleles);
            minor2 = BaseAnnot.getComplement(minor2);
        }

        nrAllelesIdentical = 0;

        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 2; b++) {
                if (firstDatasetAlleles[a] == secondDatasetAlleles[b]) {
                    nrAllelesIdentical++;
                }
            }
        }

        boolean flipalleles = false;
        if (nrAllelesIdentical != 2) {
            System.out.println("Alleles not identical!");
        } else {
            if (minor1 != minor2) {
                // error or warning or whatever
                //                        System.out.println("WARNING: minor allele is different for identical SNP: "+dSNP.getName() + ", probably due to high MAF.\nWill conform to allelic direction of dataset: "+m_gg[firstDatasetToPassQC].getSettings().name);
                //                        double[] allelefreq = dSNP.getAlleleFreq();
                //                        byte[] origAlleles = snps[firstDatasetToPassQC].getAlleles();
                //                        double[] origAlleleFreq = snps[firstDatasetToPassQC].getAlleleFreq();
                //                        System.out.println("Reference MAF:"+snps[firstDatasetToPassQC].getMAF()+"\tAssessed MAF:"+dSNP.getMAF());
                //                        for(int i=0; i<2; i++){
                //                            System.out.println("ref ds: "+m_gg[firstDatasetToPassQC].getSettings().name+"\t"+BaseAnnot.toString(origAlleles[i])+"\t("+origAlleleFreq[i]+")\tAssessed: "+m_gg[d].getSettings().name+"\t"+BaseAnnot.toString(allelesToCompare[i])+"\t("+allelefreq[i]+")");
                //                        }
                //                        System.out.println("");
                // take the orientation of the first dataset..., which is dataset

                flipalleles = true;


            }
        }
        return flipalleles;
    }
}
