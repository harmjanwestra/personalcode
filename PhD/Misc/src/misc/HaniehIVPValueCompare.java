/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.graphics.ScatterPlot;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.ProbeAnnotation;

/**
 *
 * @author harmjan
 */
public class HaniehIVPValueCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            String ivp = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/Results/2013-01-09-IVPilot40SNPs/2013-01-21-output.txtSTUDY_IV_ivp_DDMMYYYY.txt";
            String transp = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/Results/2013-01-09-IVPilot40SNPs/2013-01-21-output.txtSTUDY_IV_transp_DDMMYYYY.txt";
            String outdir = "/Volumes/iSnackHD/Data/Projects/2013-IVAnalysis/Results/2013-01-09-IVPilot40SNPs/ComparisonToTransP/";
            String probeTrans = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-EnsemblAnnotation-ProbeTranslationTable.txt";
            String probeAnnot = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-HT12v3.txt";
            HaniehIVPValueCompare.compare(ivp, transp, outdir, probeAnnot, probeTrans);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void compare(String IVP, String TRANSP, String outDir, String pbAnnot, String pbTrans) throws IOException {

        System.out.println("Loading pb trans");
        ProbeTranslation pbt = new ProbeTranslation();
        HashMap<String, String> metaToHT12v3Id = pbt.getProbeTranslation(pbTrans, "Probe", "HumanHT-12_V3_0_R2_11283641_A.txt");
        System.out.println("Loading pb annot");
        ProbeAnnotation pb = new ProbeAnnotation(pbAnnot);

        Gpio.createDir(outDir);
        TextFile tf = new TextFile(IVP, TextFile.R);
        int nrLnsIVP = tf.countLines() - 1;
        String[] elems = tf.readLineElems(TextFile.tab); // skip header
        String[] transprobes = new String[elems.length - 2];
        for (int q = 2; q < elems.length; q++) {
            transprobes[q - 2] = elems[q];
        }
        double[][] IVpvalsLog = new double[nrLnsIVP][elems.length - 2];
        double[][] TRANspvalsLog = new double[nrLnsIVP][elems.length - 2];

        String[] cisProbes = new String[nrLnsIVP];
        String[] snps = new String[nrLnsIVP];

        double[][] IVpvals = new double[nrLnsIVP][elems.length - 2];
        double[][] TRANspvals = new double[nrLnsIVP][elems.length - 2];
        elems = tf.readLineElems(TextFile.tab);
        int lnctr = 0;

        while (elems != null) {
            snps[lnctr] = elems[0];
            cisProbes[lnctr] = elems[1];
            for (int d = 2; d < elems.length; d++) {
                IVpvalsLog[lnctr][d - 2] = -Math.log(Double.parseDouble(elems[d])); 
                IVpvals[lnctr][d - 2] = Double.parseDouble(elems[d]);
                if (IVpvalsLog[lnctr][d - 2] > 400) {
                    IVpvalsLog[lnctr][d - 2] = 400;
                }
            }
            lnctr++;
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(TRANSP, TextFile.R);
        String[] transpelems = tf2.readLineElems(TextFile.tab);  // this file has three headers
        transpelems = tf2.readLineElems(TextFile.tab);
        lnctr = 0;
        while (transpelems != null) {
            for (int d = 2; d < transpelems.length; d++) {
                TRANspvalsLog[lnctr][d - 2] = -Math.log(Double.parseDouble(transpelems[d]));
                TRANspvals[lnctr][d - 2] = Double.parseDouble(transpelems[d]);
                if (TRANspvalsLog[lnctr][d - 2] > 400) {
                    TRANspvalsLog[lnctr][d - 2] = 400;
                }
            }

            transpelems = tf2.readLineElems(TextFile.tab);

            lnctr++;
        }

        tf2.close();

        for (int d = 0; d < TRANspvals.length; d++) {

            double[] xaxis2 = new double[TRANspvals[d].length];
            double[] yaxis2 = new double[TRANspvals[d].length];
            for (int i = 0; i < TRANspvals[d].length; i++) {
                xaxis2[i] = i + 1;
                if (cisProbes[d].equals(transprobes[i])) {
                    yaxis2[i] = 0;
                    IVpvals[d][i] = 0;
                    TRANspvals[d][i] = 0;
                    IVpvalsLog[d][i] = 0;
                    TRANspvalsLog[d][i] = 0;
                } else {
                    yaxis2[i] = Math.abs(IVpvalsLog[d][i] - TRANspvalsLog[d][i]);
                }
//                if (snps[d].equals("rs496173")) {
                if (yaxis2[i] > 0.5) {

                    String ht12v3cis = metaToHT12v3Id.get(cisProbes[d]);
                    String ht12v3trans = metaToHT12v3Id.get(transprobes[i]);

                    if (ht12v3cis == null || ht12v3trans == null) {
                        System.out.println("Could not find: " + cisProbes[d] + "\tor " + transprobes[i]);
                    } else {
                        Integer cisId = pb.getProbeToProbeId().get(ht12v3cis);
                        Integer transId = pb.getProbeToProbeId().get(ht12v3trans);



                        if (cisId != null && transId != null && pb.getChr()[cisId] == pb.getChr()[transId]) {
                            Integer midCis = (pb.getChrEnd()[cisId] + pb.getChrStart()[cisId]) / 2;
                            Integer midTrans = (pb.getChrEnd()[transId] + pb.getChrStart()[transId]) / 2;
                            if (Math.abs(midCis - midTrans) < 250000) {
                                System.out.println("CIS EFFECT WTFBBQ!");
                                yaxis2[i] = 0;
                                IVpvals[d][i] = 0;
                                TRANspvals[d][i] = 0;
                                IVpvalsLog[d][i] = 0;
                                TRANspvalsLog[d][i] = 0;
                            } else if (Math.abs(midCis - midTrans) < 1000000) {
                                System.out.println("Possible CIS: " + Math.abs(midCis - midTrans));
                                yaxis2[i] = 0;
                                IVpvals[d][i] = 0;
                                TRANspvals[d][i] = 0;
                                IVpvalsLog[d][i] = 0;
                                TRANspvalsLog[d][i] = 0;
                            }
                        } else {
                            System.out.println("Different CHR: " + pb.getChr()[cisId] + "\t" + pb.getChr()[transId] + "\t" + ht12v3cis + "\t" + ht12v3trans);
                        }

                    }


                    System.out.println(snps[d] + "\t" + cisProbes[d] + "\t" + (i + 2) + "\t" + transprobes[i] + "\t" + IVpvalsLog[d][i] + "\t" + TRANspvalsLog[d][i] + "\t" + IVpvals[d][i] + "\t" + TRANspvals[d][i]+"\t"+(IVpvals[d][i]-TRANspvals[d][i]));
                    System.out.println("");
                }
//                }
            }

//            ScatterPlot q = new ScatterPlot(500, 500, xaxis2, yaxis2, outDir + snps[d] + "Log.png");
//            ScatterPlot q2 = new ScatterPlot(500, 500, IVpvalsLog[d], TRANspvalsLog[d], outDir + snps[d] + "IVvsTrans.png");

            double[] xaxis = new double[TRANspvals[d].length];
            double[] yaxis = new double[TRANspvals[d].length];
            for (int i = 0; i < TRANspvals[d].length; i++) {
                xaxis[i] = i + 1;
                if (cisProbes[d].equals(transprobes[i])) {
                    yaxis[i] = 0;
                } else {
                    yaxis[i] = IVpvals[d][i] - TRANspvals[d][i];
                }

            }

//            ScatterPlot p = new ScatterPlot(500, 500, xaxis, yaxis, outDir + snps[d] + ".png");



//            ScatterPlot p = new ScatterPlot(500, 500, TRANspvals[d], IVpvals[d], outDir + snps.get(d) + ".png");
//            ScatterPlot plot = new ScatterPlot(520,520, 0.01);
//            plot.plot(TRANspvals[d], IVpvals[d]);
//            plot.draw(outDir+snps.get(d)+".png");
//            }

        }

        System.out.println("Total combos: " + (nrLnsIVP * TRANspvals[0].length));
    }
}
