/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import eqtlmappingpipeline.binarymeta.meta.MetaAnalyze;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.zip.DataFormatException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class TraitSNPGetAlleleAssessed extends MetaAnalyze {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            TraitSNPGetAlleleAssessed t = new TraitSNPGetAlleleAssessed();
            // "/Volumes/Data2/MetaAnalysisSoftware/settings/2012-06-24-CIS-40PCs-4GWAS-GeneticVectorsNotRemoved.xml"
            t.init("/Volumes/Data2/MetaAnalysisSoftware/settings/2012-03-12-TRANS-40PCs-4GWAS-GeneticVectorsNotRemoved-CisEffectsRegressedOut.xml", null, null);
            t.analyze();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public void analyze() throws IOException, DataFormatException, Exception {

        System.out.println("");
        System.out.println("Starting analysis!");

        ds = new BinaryResultDataset[m_settings.getDatasetlocations().size()];

        String[] datasets = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < m_settings.getDatasetnames().size(); i++) {
            datasets[i] = m_settings.getDatasetnames().get(i);
        }

        if (!m_settings.getOutput().endsWith("/")) {
            m_settings.setOutput(m_settings.getOutput() + "/MetaAnalysis/");
        }

        if (!Gpio.exists(m_settings.getOutput())) {
            Gpio.createDir(m_settings.getOutput());
        }

        String[] locations = new String[m_settings.getDatasetnames().size()];
        for (int i = 0; i < locations.length; i++) {
            locations[i] = m_settings.getDatasetlocations().get(i);
        }

        pvaluedistribution = null;
        eQTLBuffer = null;
        finalEQTLBuffer = null;
        nrInFinalBuffer = 0;

        uniqueProbes = new HashSet<String>();
        uniqueSNPs = new HashSet<String>();

        int numDatasets = ds.length;
        probes = new ArrayList<String>();

        snps = new ArrayList<String>();
        snpChr = new ArrayList<Byte>();
        snpChrPos = new ArrayList<Integer>();

        nrTotalSamples = 0;
        String[] probeName = probeTranslation.getProbes();
        probes.addAll(Arrays.asList(probeName));

        initdatasets(locations, 0, -1);

        TextFile out = new TextFile("/Volumes/iSnackHD/AlleleCodingForLude.txt", TextFile.W);

        int zscorecounter = 0;
        int nullctr = 0;
        for (int s = 0; s < snps.size(); s++) {
            BinaryResultSNP selectedSNP = null;
            for (int d = 0; d < ds.length; d++) {
                // get zScores.
                Integer snpId = snpTranslation[d][s];

                if (snpId != null) {

                    BinaryResultSNP snpObject = ds[d].getSnps()[snpId];

                    if (snpObject != null && selectedSNP == null) {
                        selectedSNP = snpObject;
                    }
//
//                    Float[] zscores = ds[d].readSNPZScores(snpObject);
//
//                    if (zscores == null) {
//                        System.out.println(snps.get(s) + "\thas null results!?");
//                        nullctr++;
//                    } else {
//                        for (int p = 0; p < probes.size(); p++) {
//                            Integer probeId = probeTranslationLookupTable[d][p];
//
//                            Byte probeChr = probeTranslation.getProbeChr(p);
//
//                            if (probeId != null && probeChr != null && probeChr != 22) {
//                                if (zscores[probeId] != null) {
////                                System.out.println("Probe chr: " + probeChr + "\tSNP Chr: " + snpChr.get(s));
//                                    zscorecounter++;
//                                }
//                            }
//                        }
//
//                    }








                    // 

                }
            }
            BinaryResultSNP snpObject = selectedSNP;

            StringBuilder b = new StringBuilder();
            b.append(snpObject.getName()).append("\t");
            b.append(snpObject.getChr()).append("\t");
            b.append(snpObject.getChrpos()).append("\t");
            b.append(BaseAnnot.toString(snpObject.getAlleles()[0])).append("/").append(BaseAnnot.toString(snpObject.getAlleles()[1])).append("\t");
            b.append(BaseAnnot.toString(snpObject.getAssessedAllele())).append("\t");

            out.writeln(b.toString());
        }



        System.out.println("NULLS:\t" + nullctr);

        System.out.println("Zs:\t" + zscorecounter);

        out.close();







    }
}
