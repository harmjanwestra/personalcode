package nl.harmjanwestra.playground.biogen.eqtlforestplot;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.biogen.mrvisualisation.MRPanel;
import org.checkerframework.checker.units.qual.A;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class EqtlForestPlot {
    public static void main(String[] args) {
        String qtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\eQTLs-wDatasetZScores.txt.gz";

//        qtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-transzscores\\eQTLs.txt.gz";
//        qtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-ciszscores\\eQTLs.txt.gz";
//        qtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-ciszscores\\eQTLs.txt.gz";
        qtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\Cortex-EUR-AFR-noENA\\eQTLs.txt.gz";

        EqtlForestPlot e = new EqtlForestPlot();
        String snp = "7:12244161:rs1990622:A_G";
        String gene = null;//"ENSG00000172137.19"; // "CellMapNNLS_Neuron";
//        gene = "ENSG00000172137.19"; // "CellMapNNLS_Neuron";
        String snpFileName = snp.replaceAll(":", "_");
//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\" + snpFileName + "-" + gene + ".pdf";
//        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-transzscores\\" + snpFileName + "-" + gene + ".pdf";
//        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-ciszscores\\" + snpFileName + "-" + gene + ".pdf";
//        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\trans\\rs1990622-ciszscores\\eQTLs.pdf";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2021-02-04-QTLTablesForPaper\\cellcounts\\Cortex-EUR-AFR-noENA\\" + snpFileName + "-" + gene + ".pdf";
        boolean flip = false;
        boolean sidetoside = true;
        try {
            String[] datasetOrder = new String[]{
                    "AFR-AMPAD-MSBB-V2",
                    "EUR-AMPAD-MAYO-V2",
                    "EUR-AMPAD-MSBB-V2",
                    "EUR-AMPAD-ROSMAP-V2",
                    "AFR-CMC",
                    "EUR-CMC",
                    "AFR-CMC_HBCC_set1",
                    "AFR-CMC_HBCC_set2",
                    "EUR-CMC_HBCC_set2",
                    "EUR-CMC_HBCC_set3",
                    "EUR-GTEx",
                    "EUR-GVEX",
                    "EUR-BrainGVEX-V2",
                    "AFR-LIBD_1M",
                    "EUR-LIBD_1M",
                    "AFR-LIBD_h650",
                    "EUR-LIBD_h650",
                    "EUR-NABEC-H550",
                    "EUR-NABEC-H610",
                    "EUR-TargetALS",
                    "EUR-UCLA_ASD"
            };
            datasetOrder = new String[]{
                    "AFR-AMPAD-MSBB-V2",
                    "EUR-AMPAD-MAYO-V2",
                    "EUR-AMPAD-MSBB-V2",
                    "EUR-AMPAD-ROSMAP-V2",
                    "EUR-BrainGVEX-V2",
                    "AFR-CMC",
                    "EUR-CMC",
                    "AFR-CMC_HBCC_set1",
                    "AFR-CMC_HBCC_set2",
                    "EUR-CMC_HBCC_set2",
                    "EUR-CMC_HBCC_set3",
                    "EUR-GTEx",
                    "EUR-GVEX",
                    "AFR-LIBD_1M",
                    "EUR-LIBD_1M",
                    "AFR-LIBD_h650",
                    "EUR-LIBD_h650",
                    "EUR-NABEC-H550",
                    "EUR-NABEC-H610",
                    "EUR-UCLA_ASD"
            };


            e.plot(qtlfile, snp, gene, output, flip, sidetoside, datasetOrder);
        } catch (IOException | DocumentException ioException) {
            ioException.printStackTrace();
        }
    }

    public void plot(String qtlfile, String snpquery, String genequery, String output, boolean flip, boolean sideToSide, String[] datasetOrder) throws IOException, DocumentException {

        HashMap<String, Integer> datasetMap = null;
        if (datasetOrder != null) {
            datasetMap = new HashMap<>();
            for (int i = 0; i < datasetOrder.length; i++) {
                datasetMap.put(datasetOrder[i], i);
            }
        }

        TextFile tf = new TextFile(qtlfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashMap<String, ArrayList<ForestPlotEQTL>> qtls = new HashMap<String, ArrayList<ForestPlotEQTL>>();

        Grid grid = new Grid(Grid.SIZE.LETTER, 500, 500, 10, 10, 100, 100);
        while (elems != null) {
            Double p = Double.parseDouble(elems[0]);
            String snp = elems[1];
            String gene = elems[4];

            if ((snpquery == null || snp.equals(snpquery)) && (genequery == null || genequery.equals(gene))) {
                String[] snpElems = snp.split(":");

                if (sideToSide) {
                    qtls = new HashMap<String, ArrayList<ForestPlotEQTL>>();
                }

                gene = gene.replaceAll("CellMapNNLS_", "");

                ArrayList<ForestPlotEQTL> events = new ArrayList<>();

                String metabetaStr = elems[18].replaceAll("\\(", "");
                metabetaStr = metabetaStr.replaceAll("\\)", "");
                String[] metabetaElems = metabetaStr.split(" ");
                double metabeta = Double.parseDouble(metabetaElems[0]);
                double metabetase = Double.parseDouble(metabetaElems[1]);

                String[] betas = elems[19].split(";");
                String[] samplesizes = elems[13].split(";");
                String[] datasets = elems[11].split(";");
                int sumN = 0;

                String alleles = elems[8];
                String effectAllele = elems[9];
                if (flip) {
                    String[] alElems = alleles.split("/");
                    if (alElems[0].equals(effectAllele)) {
                        effectAllele = alElems[1];
                    } else {
                        effectAllele = alElems[0];
                    }
                }

                if (datasetMap != null) {
                    for (int i = 0; i < datasetOrder.length; i++) {
                        ForestPlotEQTL qtl = new ForestPlotEQTL();
                        qtl.dataset = datasetOrder[i];
                        events.add(qtl);
                    }
                }


                for (int i = 0; i < datasets.length; i++) {
                    if (!datasets[i].equals("-")) {
                        String betastr = betas[i].replaceAll("\\(", "");
                        betastr = betastr.replaceAll("\\)", "");
                        String[] betaelems = betastr.split(" ");
                        double beta = Double.parseDouble(betaelems[0]);
                        double se = Double.parseDouble(betaelems[1]);


                        String dataset = datasets[i];
                        if (datasetMap == null || datasetMap.containsKey(dataset)) {
                            ForestPlotEQTL qtl = null;
                            if (datasetMap == null) {
                                qtl = new ForestPlotEQTL();
                            } else {
                                Integer id = datasetMap.get(dataset);
                                qtl = events.get(id);
                            }

                            qtl.dataset = datasets[i];
                            int n = Integer.parseInt(samplesizes[i]);
                            qtl.dataset = qtl.dataset.replaceAll("-V2", "");
                            if (qtl.dataset.contains("EUR")) {
                                qtl.dataset = qtl.dataset.replaceAll("EUR-", "");
                                qtl.dataset = qtl.dataset + " (EUR; n=" + n + ")";
                            } else if (qtl.dataset.contains("AFR")) {
                                qtl.dataset = qtl.dataset.replaceAll("AFR-", "");
                                qtl.dataset = qtl.dataset + " (AFR; n=" + n + ")";
                            } else {
                                qtl.dataset = qtl.dataset + " (n=" + n + ")";
                            }

                            if (flip) {
                                qtl.beta = -beta;

                            } else {
                                qtl.beta = beta;
                            }

                            qtl.se = se;
                            qtl.EA = effectAllele;
                            qtl.NONEA = alleles;
                            qtl.gene = gene;
                            qtl.snp = snp;

                            sumN += n;
                            qtl.samplesize = n;
                            if (datasetMap == null) {
                                events.add(qtl);
                            }
                        }
                    }
                }


                ForestPlotEQTL qtl = new ForestPlotEQTL();
                qtl.dataset = "Meta-Analysis";
                qtl.isMeta = true;
                if (flip) {
                    qtl.beta = -metabeta;
                } else {
                    qtl.beta = metabeta;
                }

                qtl.se = metabetase;
                qtl.EA = effectAllele;
                qtl.NONEA = alleles;
                qtl.gene = gene;
                qtl.snp = snp;
                qtl.samplesize = sumN;
                events.add(qtl);
                String genehgnc = elems[16];
                snp = snpElems[2] + " - " + alleles + " - " + effectAllele;
                String qtlstr = snp + " - " + genehgnc;
                qtls.put(qtlstr, events);

                if (sideToSide) {
                    ForestplotPanel panel = new ForestplotPanel(1, 1);
                    // panel.determineRangeOverAllEvents();
                    panel.overlapBoxes();
                    grid.addPanel(panel);
                    panel.setData(qtls);
                    panel.setShading(false);
                }

            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();


        System.out.println("plotting");

        if (!sideToSide) {
            ForestplotPanel panel = new ForestplotPanel(1, 1);
            panel.determineRangeOverAllEvents();
            panel.overlapBoxes();
            panel.setRange(new Range(-1, -1, 1, 1));
            grid.addPanel(panel);
            panel.setData(qtls);
            panel.setShading(false);
        }

        grid.draw(output);
        System.out.println("Output: " + output);

    }


}
