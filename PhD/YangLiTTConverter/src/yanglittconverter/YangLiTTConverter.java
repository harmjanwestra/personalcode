/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package yanglittconverter;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperExpressionData;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDataset;
import umcg.genetica.io.trityper.TriTyperGeneticalGenomicsDatasetSettings;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class YangLiTTConverter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        // TODO code application logic here

        String eqtldata = "/Volumes/iSnackHD/tmp/Choy - CEU/eQTLMapping-Uncorrected/eQTLProbesFDR0.05.txt";
        String outdir = "/Volumes/iSnackHD/tmp/ChoyYangLi/";

        TriTyperGeneticalGenomicsDatasetSettings settings = new TriTyperGeneticalGenomicsDatasetSettings();
        settings.cisAnalysis = true;
        settings.transAnalysis = true;
        settings.expressionLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/HapMap-Altschuler/AffymetrixExpressionDataOriginal.txt.TriTyperFormat.txt";
        settings.expressionplatform = "";
        settings.genotypeLocation = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/BloodHT12/BloodHT12ImputeTriTyper/";
        settings.genotypeToExpressionCoupling = "/Volumes/iSnackHD/tmp/Choy/GenotypeToExpressionIdOriginal-CEU.txt";
        settings.name = "Choy";
        settings.probeannotation = null;
        settings.transAnalysis = true;


        try {
            YangLiTTConverter conv = new YangLiTTConverter();
            conv.run(eqtldata, settings, outdir);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void run(String eqtlfile, TriTyperGeneticalGenomicsDatasetSettings settings, String outdir) throws IOException, Exception {
        HashSet<String> snpsToExport = new HashSet<String>();
        HashSet<String> probesToExport = new HashSet<String>();

        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            snpsToExport.add(elems[1]);
            probesToExport.add(elems[4]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        outdir = Gpio.formatAsDirectory(outdir);
        Gpio.createDir(outdir);

        TextFile genout = new TextFile(outdir + "genotypes.txt", TextFile.W);
        TextFile genoutAlleles = new TextFile(outdir + "genotypesAlleles.txt", TextFile.W);

        TriTyperGeneticalGenomicsDataset ds = new TriTyperGeneticalGenomicsDataset(settings);
//        TriTyperGenotypeData gds = ds.getGenotypeData();
//        SNPLoader loader = gds.createSNPLoader();
//        String[] individuals = gds.getIndividuals();
//        String header = "-";
//        for (int i = 0; i < individuals.length; i++) {
//            String ind = individuals[i];
//            if (ds.getGenotypeToExpressionCouplings().get(ind) != null) {
//                if (gds.getIsIncluded()[i]) {
//                    header += "\t" + individuals[i];
//                }
//            }
//        }
//        genout.writeln(header);
//
//
//
//        for (String s : snpsToExport) {
//
//            Integer SNPid = gds.getSnpToSNPId().get(s);
//            if (SNPid != null) {
//                SNP snpObj = gds.getSNPObject(SNPid);
//                loader.loadGenotypes(snpObj);
//                String output = s;
//                String alleleout = s + "\t" + BaseAnnot.toString(snpObj.getAlleles()[0]) + "\t" + BaseAnnot.toString(snpObj.getAlleles()[1]);
//                genoutAlleles.writeln(alleleout);
//                short[] genotypes = snpObj.getGenotypes();
//
//                for (int i = 0; i < individuals.length; i++) {
//                    String ind = individuals[i];
//                    if (ds.getGenotypeToExpressionCouplings().get(ind) != null) {
//                        if (gds.getIsIncluded()[i]) {
//                            if (genotypes[i] == -1) {
//                                output += "\tNA";
//                            } else {
//                                output += "\t" + genotypes[i];
//                            }
//                        }
//                    }
//
//                }
//                genout.writeln(output);
//                snpObj.clearGenotypes();
//            }
//        }
//
//        loader.close();
//        genout.close();
//        genoutAlleles.close();
//
//
//
//        TextFile expout = new TextFile(outdir + "expression.txt", TextFile.W);
//        String expheader = "-";
//        TriTyperExpressionData exp = ds.getExpressionData();
//        String[] inds = exp.getIndividuals();
//
//        for (int i = 0; i < inds.length; i++) {
//            String ind = inds[i];
//
//            expheader += "\t" + ind;
//        }
//
//        expout.writeln(expheader);
//
//        double[][] matrix = exp.getMatrix();
//        String[] probes = exp.getProbes();
////        for (String probe : probesToExport) {
////            Integer probeId = exp.getProbeToId().get(probe);
////            if (probeId != null) {
////                String output = probe;
////                for (int i = 0; i < inds.length; i++) {
////                    
////                    if (ds.getExpressionToGenotypeIdHash().get(i) != null) {
////                        probe+="\t"+matrix[probeId][i];
////                    }
////                }
////                expout.writeln(output);
////            }
////        }
//        for (int p = 0; p < probes.length; p++) {
//            String probe = probes[p];
//            String output = probe;
//            for (int i = 0; i < inds.length; i++) {
////                if (ds.getExpressionToGenotypeIdHash().get(i) != null) {
//                    output += "\t" + matrix[p][i];
////                }
//            }
//            expout.writeln(output);
//        }
//        expout.close();
    }
}
