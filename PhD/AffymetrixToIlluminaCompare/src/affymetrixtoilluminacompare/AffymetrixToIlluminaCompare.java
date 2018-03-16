/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package affymetrixtoilluminacompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.math.stats.*;
import umcg.genetica.util.RankDoubleArray;

/**
 *
 * @author harmjan
 */
public class AffymetrixToIlluminaCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            TextFile outfile = new TextFile("/Volumes/iSnackHD/Data/Projects/AffyWithIlluminaIntegration/2012-06-19-ComparisonBetweenHapMap-StrangerAndAltschuler-QNorm.txt", TextFile.W);
            TextFile yorfile = new TextFile("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Stranger/IDS_CEU.txt", TextFile.R);

            String[] ids = yorfile.readAsArray();

            outfile.writeln(ids.length + "\tindividuals preselected from file: " + yorfile.getFileName());
            yorfile.close();

            HashSet<String> allYORSamples = new HashSet<String>();
            allYORSamples.addAll(Arrays.asList(ids));

            TriTyperExpressionData affy = new TriTyperExpressionData();
            affy.load("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/AltschulerNormalized/ExpressionData.txt", null, null, true);

	    Log2Transform.log2transform(affy.getMatrix());
	    QuantileNormalization.quantilenormalize(affy.getMatrix());




            TriTyperExpressionData illumina = new TriTyperExpressionData();
            illumina.load("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Stranger/ExpressionData.txt", null, null, true);

	    Log2Transform.log2transform(illumina.getMatrix());
	    QuantileNormalization.quantilenormalize(illumina.getMatrix());

            String[] affyInds = affy.getIndividuals();
            HashSet<String> samplesInAffyAndYOR = new HashSet<String>();
            for (String s : affyInds) {
                if (allYORSamples.contains(s)) {
                    samplesInAffyAndYOR.add(s);
                }
            }

            System.out.println(samplesInAffyAndYOR.size() + "\tsamples shared between preselected set and Affy data");
            outfile.writeln(samplesInAffyAndYOR.size() + "\tsamples shared between preselected set and Affy data");

            String[] samplesInAffyAndYORArr = samplesInAffyAndYOR.toArray(new String[0]);

            ArrayList<String> samplesInBothAffyAndIllumina = new ArrayList<String>();
            for (String s : samplesInAffyAndYORArr) {
                if (illumina.getIndividualToId().get(s) != null) {
                    samplesInBothAffyAndIllumina.add(s);
                }
            }

            System.out.println(samplesInBothAffyAndIllumina.size() + "\tsamples in Affy, preselected set AND Illumina data.");



            int sharedsamples = samplesInBothAffyAndIllumina.size();
            outfile.writeln(sharedsamples + "\tsamples in Affy, preselected set AND Illumina data.");
            int[] samplemapaffy = new int[sharedsamples];
            int[] samplemapillumina = new int[sharedsamples];
            for (int i = 0; i < samplesInBothAffyAndIllumina.size(); i++) {
                String smpl = samplesInBothAffyAndIllumina.get(i);
                samplemapaffy[i] = affy.getIndividualId(smpl);
                samplemapillumina[i] = illumina.getIndividualId(smpl);
            }

//	    System.out.println(sharedsamples + " samples shared");

            // count the number of times unique genes are tagged by affy probes.
            String[] affyProbeGeneAnnotation = affy.getAnnotation();
            HashMap<String, Integer> probeAnnotationCounter = new HashMap<String, Integer>();
            for (String annot : affyProbeGeneAnnotation) {
                Integer ctr = probeAnnotationCounter.get(annot);
                if (ctr == null) {
                    ctr = 1;
                } else {
                    ctr++;
                }
                probeAnnotationCounter.put(annot, ctr);
            }

            ArrayList<Integer> genesWithSingleProbe = new ArrayList<Integer>();
            HashSet<String> uniqueGenesOnAffy = new HashSet<String>();
            uniqueGenesOnAffy.addAll(Arrays.asList(affyProbeGeneAnnotation));
            HashMap<String, ArrayList<Integer>> geneToProbeMapOnAffy = new HashMap<String, ArrayList<Integer>>();
            for (int i = 0; i < affyProbeGeneAnnotation.length; i++) {
                ArrayList<Integer> probeIds = geneToProbeMapOnAffy.get(affyProbeGeneAnnotation[i]);
                if (probeIds == null) {
                    probeIds = new ArrayList<Integer>();
                }
                probeIds.add(i);
                geneToProbeMapOnAffy.put(affyProbeGeneAnnotation[i], probeIds);
            }

            System.out.println(genesWithSingleProbe.size() + "\tAffy genes with single probe.");
            outfile.writeln(genesWithSingleProbe.size() + "\tAffy genes with single probe.");

            HashMap<String, ArrayList<Integer>> geneToProbeOnIllumina = new HashMap<String, ArrayList<Integer>>();
            String[] illuminaProbeGeneAnnotation = illumina.getAnnotation();
            for (int i = 0; i < illuminaProbeGeneAnnotation.length; i++) {
                ArrayList<Integer> probeIds = geneToProbeOnIllumina.get(illuminaProbeGeneAnnotation[i]);
                if (probeIds == null) {
                    probeIds = new ArrayList<Integer>();
                }
                probeIds.add(i);
                geneToProbeOnIllumina.put(illuminaProbeGeneAnnotation[i], probeIds);
            }

            int probesIllumina = illuminaProbeGeneAnnotation.length;
            int probeAffy = affy.getProbes().length;

            System.out.println(probesIllumina + " illumina probes. " + probeAffy + " affy probes");

            int ctr = 0;

            outfile.writeln("Now testing.");

            outfile.writeln("AffyProbe\tNrProbesOnAffyForGene\tIlluminaProbe\tNrProbesOnIlluminaForGene\tgene\tWilcoxonP\tWilcoxonAUC\tPearsonR\tPearsonR2\tSpearmanR\tSpearmanR2");
            String[] uniqueGenesOnAffyArray = uniqueGenesOnAffy.toArray(new String[0]);
            for (int p = 0; p < uniqueGenesOnAffyArray.length; p++) {

                ArrayList<Integer> probesForGeneOnAffy = geneToProbeMapOnAffy.get(uniqueGenesOnAffyArray[p]);
                ArrayList<Integer> probesForGeneOnIllumina = geneToProbeOnIllumina.get(uniqueGenesOnAffyArray[p]);
                if (probesForGeneOnAffy == null || probesForGeneOnIllumina == null || probesForGeneOnAffy.isEmpty() || probesForGeneOnIllumina.isEmpty()) {

//                    outfile.writeln(uniqueGenesOnAffyArray[p] + "\thas no probes on either Affy of Illumina");



                } else {
                    for (Integer affyp : probesForGeneOnAffy) {
                        for (Integer illup : probesForGeneOnIllumina) {
                            String gene = uniqueGenesOnAffyArray[p];
                            double[] affyvals = new double[sharedsamples];
                            double[] illuvals = new double[sharedsamples];

                            double[] affyOrigVals = affy.getMatrix()[affyp];
                            double[] illuOrigVals = illumina.getMatrix()[illup];


                            for (int i = 0; i < sharedsamples; i++) {
                                affyvals[i] = affyOrigVals[samplemapaffy[i]];
                                illuvals[i] = illuOrigVals[samplemapillumina[i]];
//			System.out.println(samplemapaffy[i] + "\t" + affyvals[i] + "\t" + samplemapillumina[i] + "\t" + illuvals[i]);
                            }

                            WilcoxonMannWhitney mwm = new WilcoxonMannWhitney();
                            double pval = mwm.returnWilcoxonMannWhitneyPValue(illuvals, affyvals);

                            double auc = mwm.getAUC();
                            RankDoubleArray rda = new RankDoubleArray();
                            double[] illuminaRanks = rda.rank(illuvals);
                            double[] affyranks = rda.rank(affyvals);
                            double correlation = Correlation.correlate(illuvals, affyvals);
                            double spearman = Correlation.correlate(illuminaRanks, affyranks);

                            String ln = affy.getProbes()[affyp] + "\t"
                                    + probesForGeneOnAffy.size()+ "\t"
                                    + illumina.getProbes()[illup] + "\t"
                                    + probesForGeneOnIllumina.size()+ "\t"
                                    + gene + "\t"
                                    + pval + "\t"
                                    + auc + "\t"
                                    + correlation + "\t"
                                    + (correlation * correlation) + "\t"
                                    + spearman + "\t"
                                    + (spearman * spearman);
                            System.out.println(ln);
                            outfile.writeln(ln);
                        }
                    }
                }

//                Integer affyProbeId = genesWithSingleProbe.get(p);
//
//                String hugoAffy = affyProbeGeneAnnotation[affyProbeId];
//                Integer id = geneToProbeOnIllumina.get(hugoAffy);
////		System.out.println(probe + "\t" + id);
////		System.out.println(affy.getAnnotation()[probe]);
////		System.out.println(affy.getAnnotation()[id]);
//
//                if (id != null) {
//
//                    
//
//
//                }
            }

            outfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
