package nl.harmjanwestra.playground.cis;

import com.itextpdf.text.DocumentException;
import eqtlmappingpipeline.binarymeta.meta.ZScoreComparison;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.graphics.panels.SpacerPanel;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

public class CompareMultipleTissues {


    public static void main(String[] args) {

        CompareMultipleTissues c = new CompareMultipleTissues();
        String[] tissues = new String[]{
                "Cerebellum",
                "Cortex",
                "MotorCortex",
                "SpinalCord"
        };
        String basedir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\ciseqtls\\nopc";
        String compdir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\ciseqtls\\tissuerepl";
        String outputfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\ciseqtls\\tissuecomp.png";
        try {
//            c.run(tissues, basedir, compdir, outputfilename);
            String eqtlgenb38 = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";
            basedir = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\ciseqtls\\bloodcombos";
            outputfilename = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-TargetALS\\ciseqtls\\eqtlgencom-fdr.png";
            c.compareWithEQTLGen(basedir, eqtlgenb38, tissues, outputfilename);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }
    }

    private void compareWithEQTLGen(String basedir, String eqtlgenb38, String[] tissues, String outputfilename) throws IOException, DocumentException {
        Range r = new Range(-201, -25, 201, 25);
        String tissue1 = "eQTLGen";
        QTLTextFile tf = new QTLTextFile(eqtlgenb38, QTLTextFile.R);
        EQTL[] t1eqtls = tf.read();
        HashMap<String, EQTL> eqtlmap = new HashMap<>();
        for (EQTL e : t1eqtls) {
            eqtlmap.put(e.getRsName() + "-" + e.getProbe(), e);
        }
        tf.close();
        System.out.println(eqtlmap.size() + " eqtls loaded from " + eqtlgenb38);


        Grid g = new Grid(250, 250, 1, tissues.length, 100, 100);
        for (int i = 0; i < tissues.length; i++) {
            String tissue2 = tissues[i];
            QTLTextFile tf2 = new QTLTextFile(basedir + "/" + tissue2 + "/eQTLsFDR0.05.txt.gz", QTLTextFile.R);
            EQTL[] t2eqtls = tf2.read();
            tf2.close();
            System.out.println(t2eqtls.length + " eqtls loaded for " + tissue2);

            ArrayList<Double> x = new ArrayList<Double>();
            ArrayList<Double> y = new ArrayList<Double>();

            // compare directions
            int samedir = 0;
            for (EQTL e : t2eqtls) {
                EQTL t1eqtl = eqtlmap.get(e.getRsName() + "-" + e.getProbe());
                if (t1eqtl != null) {

                    Boolean flip = BaseAnnot.flipalleles(t1eqtl.getAlleles(), t1eqtl.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
                    if (flip != null) {
                        Double z2 = e.getZscore();
                        if (flip) {
                            z2 *= -1;
                        }

                        x.add(t1eqtl.getZscore());
                        y.add(z2);
                        if ((z2 > 0 && t1eqtl.getZscore() > 0) || (z2 < 0 && t1eqtl.getZscore() < 0)) {
                            samedir++;
                        }
                    }

                }
            }

            ScatterplotPanel p = new ScatterplotPanel(1, 1);

            p.setAlpha(0.5f);
            p.setPlotElems(true, false);
            p.setLabels(tissue1, tissue2);
            double perc = (double) samedir / y.size();
            perc *= 100;
            DecimalFormat df = new DecimalFormat("##.##");
            p.setTitle(y.size() + " shared, " + df.format(perc) + "% same direction");
            p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
            p.setDataRange(r);

            g.addPanel(p);
        }
        g.draw(outputfilename);

    }

    public void run(String[] tissues, String basedir, String compdir, String outputfilename) throws IOException, DocumentException {

        Grid g = new Grid(250, 250, tissues.length, tissues.length, 100, 100);
        Range r = new Range(-25, -25, 25, 25);
        for (int i = 0; i < tissues.length; i++) {
            String tissue1 = tissues[i];
            QTLTextFile tf = new QTLTextFile(basedir + "/" + tissue1 + "/eQTLs.txt.gz", QTLTextFile.R);
            EQTL[] t1eqtls = tf.read();
            HashMap<String, EQTL> eqtlmap = new HashMap<>();
            for (EQTL e : t1eqtls) {
                eqtlmap.put(e.getRsName() + "-" + e.getProbe(), e);

            }
            tf.close();

            for (int j = 0; j < tissues.length; j++) {
                System.out.println(i + "\t" + j);
                if (j == i) {
                    g.addPanel(new SpacerPanel(1, 1));
                } else {
                    String tissue2 = tissues[j];
                    QTLTextFile tf2 = new QTLTextFile(compdir + "/" + tissue2 + "-" + tissue1 + "/eQTLs.txt.gz", QTLTextFile.R);
                    EQTL[] t2eqtls = tf2.read();
                    tf2.close();

                    ArrayList<Double> x = new ArrayList<Double>();
                    ArrayList<Double> y = new ArrayList<Double>();

                    // compare directions
                    int samedir = 0;
                    for (EQTL e : t2eqtls) {
                        EQTL t1eqtl = eqtlmap.get(e.getRsName() + "-" + e.getProbe());
                        if (t1eqtl != null) {

                            Boolean flip = BaseAnnot.flipalleles(t1eqtl.getAlleles(), t1eqtl.getAlleleAssessed(), e.getAlleles(), e.getAlleleAssessed());
                            if (flip != null) {
                                Double z2 = e.getZscore();
                                if (flip) {
                                    z2 *= -1;
                                }

                                x.add(t1eqtl.getZscore());
                                y.add(z2);
                                if ((z2 > 0 && t1eqtl.getZscore() > 0) || (z2 < 0 && t1eqtl.getZscore() < 0)) {
                                    samedir++;
                                }
                            }

                        }
                    }

                    ScatterplotPanel p = new ScatterplotPanel(1, 1);

                    p.setAlpha(0.5f);
                    p.setPlotElems(true, false);
                    p.setLabels(tissue1, tissue2);
                    double perc = (double) samedir / y.size();
                    perc *= 100;
                    DecimalFormat df = new DecimalFormat("##.##");
                    p.setTitle(y.size() + " shared, " + df.format(perc) + "% same direction");
                    p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
                    p.setDataRange(r);

                    g.addPanel(p);

                }
            }


        }

        g.draw(outputfilename);

    }

}
