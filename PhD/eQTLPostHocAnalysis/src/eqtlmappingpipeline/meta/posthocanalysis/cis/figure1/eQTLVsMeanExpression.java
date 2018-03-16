/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.cis.figure1;

import eqtlmappingpipeline.util.StringArrayByDoubleArraySort;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.eQTLTextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;

/**
 *
 * @author harm-jan
 */
public class eQTLVsMeanExpression {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
//            eQTLVsMeanExpression.run(
//                    "D:\\Skydrive\\MetaAnalysisAnnotationFiles\\2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt-WithHT12Identifiers.txt",
//                    "D:\\Skydrive\\ExpressionData\\ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//                    "D:\\Work\\2012-08-08-cis-PostQC\\PostQC\\eQTLProbesFDR0.05.txt",
//                    "D:\\Work\\2012-08-08-cis-PostQC\\PostQC\\2012-08-15-MeanExpressionOfProbes.txt");

            // cis
//            eQTLVsMeanExpression.run(
//                    "D:\\Skydrive\\MetaAnalysisAnnotationFiles\\2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt-WithHT12Identifiers.txt",
//                    "D:\\Skydrive\\ExpressionData\\ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//                    "D:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.05.txt",
//                    "D:\\SkyDrive\\latesteQTLs\\cisProbeLevelFDR0.05-MeanExpressionOfProbes.txt",
//                    "D:\\SkyDrive\\latesteQTLs\\preQC\\eQTLs.txt.gz", // add set of tested probes here
//                    null,
//                    "D:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", // probe translation
//                    true);
            
            eQTLVsMeanExpression.run(
                    "D:\\Skydrive\\MetaAnalysisAnnotationFiles\\2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt-WithHT12Identifiers.txt",
                    "D:\\Skydrive\\ExpressionData\\ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
                    "D:\\SkyDrive\\latesteQTLs\\transFDR0.05.txt.gz",
                    "D:\\SkyDrive\\latesteQTLs\\transVsGeneExpression\\transFDR0.05-MeanExpressionOfProbes.txt",
                    null, // add set of tested probes here
                    "D:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt",
                    "D:\\SkyDrive\\MetaAnalysisAnnotationFiles\\2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt", // probe translation
                    true);
//            eQTLVsMeanExpression.run(
//                    "/Volumes/iSnackHD/SkyDrive/MetaAnalysisAnnotationFiles//2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt-WithHT12Identifiers.txt",
//                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12OriginalExpressionData4SNPPCRemoved/ExpressionData.txt.gz.QuantileNormalized.Log2Transformed.txt.gz",
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/PostQC/eQTLProbesFDR0.05.txt",
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQCFDRFilteredForProbes/PostQC/2012-09-21-MeanExpressionOfProbes-WithGeneNameAnnotation-0PCsRemoved.txt",
//                    "/Volumes/iSnackHD/eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/PreQC/eQTLs.txt.gz",
//                    "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt",
//                    true);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void run(String ht12,
            String expressionFile,
            String eqtlfilename,
            String output,
            String eQTLFileTestOnlyTheseProbes,
            String eQTLFileTestOnlyTheseProbesList,
            String probeTranslation, boolean abs) throws IOException {



        HashSet<String> probesToTestFromeQTLFile = null;
        if (eQTLFileTestOnlyTheseProbes != null) {
            probesToTestFromeQTLFile = new HashSet<String>();
            TextFile f = new TextFile(eQTLFileTestOnlyTheseProbes, TextFile.R);
            String[] felems = f.readLineElems(TextFile.tab);
            felems = f.readLineElems(TextFile.tab);
            while (felems != null) {
                String probe = felems[4];
                probesToTestFromeQTLFile.add(probe);
                felems = f.readLineElems(TextFile.tab);
            }
            f.close();
        } else if (eQTLFileTestOnlyTheseProbesList != null) {
            probesToTestFromeQTLFile = new HashSet<String>();
            TextFile f = new TextFile(eQTLFileTestOnlyTheseProbesList, TextFile.R);
            probesToTestFromeQTLFile.addAll(f.readAsArrayList());
            f.close();
        }

        HashMap<String, String> metaProbeIdToGeneName = new HashMap<String, String>();

        TextFile probetf = new TextFile(probeTranslation, TextFile.R);

        String[] probetfelems = probetf.readLineElems(TextFile.tab);
        while (probetfelems != null) {
            String probe = probetfelems[0];
            String gene = probetfelems[4];
            metaProbeIdToGeneName.put(probe, gene);
            probetfelems = probetf.readLineElems(TextFile.tab);
        }
        probetf.close();

        TextFile tf = new TextFile(ht12, TextFile.R);
        String[] ht12Elems = tf.readLineElems(TextFile.tab);
        HashMap<String, String> metaToProbeId = new HashMap<String, String>();
        HashMap<String, String> probeToMetaId = new HashMap<String, String>();

        HashSet<String> ht12ProbesThatHaveBeenTested = new HashSet<String>();
        HashMap<String, String> metaZ = new HashMap<String, String>();
        while (ht12Elems != null) {
            metaToProbeId.put(ht12Elems[0], ht12Elems[1]);
            probeToMetaId.put(ht12Elems[1], ht12Elems[0]);
            if (probesToTestFromeQTLFile == null || probesToTestFromeQTLFile.contains(ht12Elems[0])) {
                ht12ProbesThatHaveBeenTested.add(ht12Elems[1]);
            }
            ht12Elems = tf.readLineElems(TextFile.tab);
        }

        tf.close();

        HashSet<String> probesThatHaveEQTLEffect = new HashSet<String>();
        TextFile eqtlfile = new TextFile(eqtlfilename, TextFile.R);
        String[] elems = eqtlfile.readLineElems(TextFile.tab); // header
        elems = eqtlfile.readLineElems(TextFile.tab);
        while (elems != null) {

            String z = elems[eQTLTextFile.METAZ];
            String probe = elems[4];
            metaZ.put(metaToProbeId.get(probe), z);

            probesThatHaveEQTLEffect.add(probe);
            elems = eqtlfile.readLineElems(TextFile.tab);
        }
        eqtlfile.close();

        HashSet<String> ht12probesThatHaveEQTLEFfect = new HashSet<String>();
        // convert to ht12 numbers
        for (String s : probesThatHaveEQTLEffect) {
            ht12probesThatHaveEQTLEFfect.add(metaToProbeId.get(s));
        }
        DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(expressionFile);

        double sum = 0;

        for (int i = 0; i < ds.nrCols; i++) {
            double intermediate = 0;
            int ctr = 0;
            for (int j = 0; j < ds.nrRows; j++) {
                String probeName = ds.rowObjects.get(j);
                if (ht12ProbesThatHaveBeenTested.contains(probeName)) {
                    intermediate += ds.rawData[j][i];
                    ctr++;
                }
            }
            sum += (intermediate / ctr);
        }

        double meanOverallExpression = sum / ds.nrCols;
        System.out.println("Mean expression: " + meanOverallExpression);

        // determine variance
        double overallVariance = 0;
        int ctr = 0;
        for (int i = 0; i < ds.nrCols; i++) {
            for (int j = 0; j < ds.nrRows; j++) {
                String probeName = ds.rowObjects.get(j);
                if (ht12ProbesThatHaveBeenTested.contains(probeName)) {
                    double d = ds.rawData[j][i];
                    double diff = d - meanOverallExpression;
                    overallVariance += (diff) * (diff);
                    ctr++;
                }
            }
        }
        overallVariance /= (ctr - 1);

        double sd = Math.sqrt(overallVariance);
        System.out.println("SD: " + sd);

        TextFile tfout = new TextFile(output, TextFile.W);
        // now determine the Z-score from the mean expression for eQTL probes
        tfout.writeln("Mean overall expression: " + meanOverallExpression);
        tfout.writeln("SD: " + sd);

        tfout.writeln("nr\tmetaId\tht12Id\tgenename\tisEQTL\tmeanExp\tstdev\tZ\tCV\tmetaZ");
        int nrCisEQTLProbes = 0;
        String[] probeNamesThatHaveBeenTestedInOrder = new String[ht12ProbesThatHaveBeenTested.size()];
        double[] probeMeansThatHaveBeenTestedInOrder = new double[ht12ProbesThatHaveBeenTested.size()];
        int probeCTR = 0;
        HashMap<String, Double> probeToCV = new HashMap<String, Double>();
        HashMap<String, Double> probeToSD = new HashMap<String, Double>();
        HashMap<String, Double> probeToZ = new HashMap<String, Double>();

        for (int i = 0; i < ds.nrRows; i++) {
            String probeName = ds.rowObjects.get(i);
            if (ht12ProbesThatHaveBeenTested.contains(probeName)) {

                double meanExp = 0;
                for (int j = 0; j < ds.nrCols; j++) {
                    meanExp += ds.rawData[i][j];
                }

                double probeMeanExp = Descriptives.mean(ds.rawData[i]);
                double probeSD = Math.sqrt(Descriptives.variance(ds.rawData[i], probeMeanExp));

                double cv = probeSD / probeMeanExp;



                meanExp /= ds.nrCols;
                double z = (meanExp - meanOverallExpression) / sd;
                boolean eQTLProbe = false;
                if (ht12probesThatHaveEQTLEFfect.contains(probeName)) {
                    eQTLProbe = true;
                    nrCisEQTLProbes++;
                }
                probeMeansThatHaveBeenTestedInOrder[probeCTR] = meanExp;
                probeNamesThatHaveBeenTestedInOrder[probeCTR] = probeName;

                probeToCV.put(probeName, cv);
                probeToSD.put(probeName, sd);
                probeToZ.put(probeName, z);

                String meta = probeToMetaId.get(probeName);
                String genename = metaProbeIdToGeneName.get(meta);

//                System.out.println(probeToMetaId.get(probeName) + "\t" + probeName + "\t" + eQTLProbe + "\t" + meanExp + "\t" + probeSD + "\t" + z + "\t" + cv + "\t" + metaZ.get(probeName));
                tfout.writeln(probeCTR + "\t" + probeToMetaId.get(probeName) + "\t" + probeName + "\t" + genename + "\t" + eQTLProbe + "\t" + meanExp + "\t" + probeSD + "\t" + z + "\t" + cv + "\t" + metaZ.get(probeName));
                probeCTR++;
            }
        }


        StringArrayByDoubleArraySort.sort(probeMeansThatHaveBeenTestedInOrder, probeNamesThatHaveBeenTestedInOrder);

//        for (int i = 0; i < 100; i++) {
//            System.out.println(i + "\t" + probeMeansThatHaveBeenTestedInOrder[i]);
//        }

        int increment = 2000;
        System.out.println(probeMeansThatHaveBeenTestedInOrder.length + " probes tested");
        System.out.println(ht12probesThatHaveEQTLEFfect.size() + " probes with eQTL");
        int previncrement = 0;
        int eQTLCounterCumulative = 0;
        for (int upperlimit = probeMeansThatHaveBeenTestedInOrder.length - 1; upperlimit > 0; upperlimit -= increment) {
            if (upperlimit - increment < 0) {
                int remainder = Math.abs(upperlimit - increment);
                increment -= remainder;
            }
            int lowerlimit = upperlimit - increment;

            for (int i = upperlimit; i > lowerlimit; i--) {
                if (ht12probesThatHaveEQTLEFfect.contains(probeNamesThatHaveBeenTestedInOrder[i])) {
                    eQTLCounterCumulative++;
                }
            }

            TextFile probesInBin = new TextFile(output + "-" + upperlimit + ".txt", TextFile.W);
            probesInBin.writeln("ht12v3Id\tmetaId\tHUGO\tmean\tsd\tz\tcv\tisEQTL");
            TextFile probeseQTLInBin = new TextFile(output + "-" + upperlimit + "-eQTLs.txt", TextFile.W);
            probeseQTLInBin.writeln("ht12v3Id\tmetaId\tHUGO\tmean\tsd\tz\tcv\tisEQTL");
            int nrEQTLsInBin = 0;
            for (int i = upperlimit; i > lowerlimit; i--) {
                String probe = probeNamesThatHaveBeenTestedInOrder[i];
                String meta = probeToMetaId.get(probe);
                String genename = metaProbeIdToGeneName.get(meta);
                double mean = probeMeansThatHaveBeenTestedInOrder[i];
                sd = probeToSD.get(probe);
                double cv = probeToCV.get(probe);
                double z = probeToZ.get(probe);
                String outputStr = probe + "\t" + meta + "\t" + genename + "\t" + mean + "\t" + sd + "\t" + z + "\t" + cv;
                if (ht12probesThatHaveEQTLEFfect.contains(probe)) {
                    nrEQTLsInBin++;
                    probeseQTLInBin.writeln(outputStr + "\tTRUE");
                } else {
                    probesInBin.writeln(outputStr + "\tFALSE");
                }
            }
            probesInBin.close();
            probeseQTLInBin.close();
            System.out.println(previncrement + "\t" + increment + "\t" + eQTLCounterCumulative + "\t" + nrEQTLsInBin);
            previncrement += increment;
        }

//        int nrEQTLsInBin = 0;
//
//        for (int i = upperlimit - increment; i < upperlimit; i++) {
//            if (ht12probesThatHaveEQTLEFfect.contains(probeNamesThatHaveBeenTestedInOrder[i])) {
//                nrEQTLsInBin++;
//            }
//        }
//        System.out.println(upperlimit + "\t" + eQTLCounter + "\t" + nrEQTLsInBin);


        System.out.println(nrCisEQTLProbes);
        tfout.close();
    }
}
