package nl.harmjanwestra.playground.biogen.freeze2dot1.gtex;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class Replicate {


    public static void main(String[] args) {

        String eqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cortex-EUR\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
        String gtex = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_EUR_eQTL_all_associations-flat\\";
        String gtexsignificant = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_eQTL\\";
        String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cortex-EUR\\gtexV8Comparison\\";
        String outSignificant = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cortex-EUR\\gtexV8ComparisonSignificant\\";
        double pvalthreshold = 5e-6;


        Replicate r = new Replicate();

//            r.run(eqtlfile, gtex, out);
//            r.run(eqtlfile, gtex, out, false, pvalthreshold);
//            r.run(eqtlfile, gtexsignificant, outSignificant, true, pvalthreshold);
        eqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cerebellum-EUR\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
        out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cerebellum-EUR\\gtexV8Comparison\\";
        outSignificant = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-cerebellum-EUR\\gtexV8ComparisonSignificant\\";

//            r.run(eqtlfile, gtexsignificant, outSignificant, true, pvalthreshold);
//            r.run(eqtlfile, gtex, out, false, pvalthreshold);

        eqtlfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-basalganglia-EUR\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
        out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-basalganglia-EUR\\gtexV8Comparison\\";
        outSignificant = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-04-13-eqtls\\cis-basalganglia-EUR\\gtexV8ComparisonSignificant\\";
//            r.run(eqtlfile, gtex, out, false, pvalthreshold);
//            r.run(eqtlfile, gtexsignificant, outSignificant, true, pvalthreshold);


        // 2020-05-26 release

        String outall = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisReplicationGTEx\\";


        String[] eqtlfiles = new String[]{
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Cortex-EUR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Basalganglia-EUR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Cerebellum-EUR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Cortex-AFR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Cortex-EAS-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Hippocampus-EUR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz",
                "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cis\\2020-05-26-Spinalcord-EUR-Iteration1-eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        };
        String[] tissue = new String[]{
                "Cortex-EUR-Iteration1",
                "Basalganglia-EUR-Iteration1",
                "Cerebellum-EUR-Iteration1",
                "Cortex-AFR-Iteration1",
                "Cortex-EAS-Iteration1",
                "Hippocampus-EUR-Iteration1",
                "Spinalcord-EUR-Iteration1",
        };
        IntStream.range(0, eqtlfiles.length).parallel().forEach(v -> {
            String efile = eqtlfiles[v];
            String outallprefix = outall + tissue[v] + "-all";
            String outsigprefix = outall + tissue[v] + "-sig";
//            try {
//                r.run(efile, gtex, outallprefix, false, pvalthreshold, false);
//                r.run(efile, gtexsignificant, outsigprefix, true, pvalthreshold, false);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }

        });

        String mergedout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisReplicationGTEx\\2020-05-26-AllTissues";
        try {
            r.mergeFiles(outall, tissue, mergedout);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private void mergeFiles(String outall, String[] tissueList, String mergedout) throws IOException {

        HashSet<String> gtexTissues = new HashSet<>();
        HashMap<String, HashMap<String, Pair<Integer, Double>>> data = new HashMap<>();
        for (int t = 0; t < tissueList.length; t++) {
            String outsigprefix = outall + tissueList[t] + "-sig-summary.txt";
            TextFile tf = new TextFile(outsigprefix, TextFile.R);
            tf.readLine();//header
            String[] elems = tf.readLineElems(TextFile.tab);
            HashMap<String, Pair<Integer, Double>> tissuedata = new HashMap<>();
            while (elems != null) {
                String tissue = elems[0];
                Double perc = Double.parseDouble(elems[3]);
                Integer shared = Integer.parseInt(elems[1]);
                tissuedata.put(tissue, new Pair<Integer, Double>(shared, perc));
                gtexTissues.add(tissue);
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            data.put(tissueList[t], tissuedata);
        }

        TextFile sigmerged = new TextFile(mergedout + "-merged-sig.txt", TextFile.W);
        String header = "-";
        for (int t = 0; t < tissueList.length; t++) {
            header += "\t" + tissueList[t] + "-Shared\t" + tissueList[t] + "-Concordant";
        }

        sigmerged.writeln(header);
        for (String tissue : gtexTissues) {
            String ln = tissue;
            for (int t = 0; t < tissueList.length; t++) {
                HashMap<String, Pair<Integer, Double>> d = data.get(tissueList[t]);
                if (d != null) {
                    Pair<Integer, Double> pair = d.get(tissue);

                    if (pair != null) {
                        ln += "\t" + pair.getLeft() + "\t" + pair.getRight();
                    } else {
                        ln += "\t" + 0 + "\t" + 0;
                    }
                } else {
                    ln += "\t" + 0 + "\t" + 0;
                }
            }
            sigmerged.writeln(ln);
        }
        sigmerged.close();

    }


    class EQTL {
        String alleles;
        String assessed;
        double z;

        double otherZ = Double.NaN;

    }

    public void run(String eqtlfile, String gtexfolder, String output, boolean GTEXsignificantFolder,
                    double pvalthreshold, boolean writeTissueOutput) throws IOException {

        Gpio.createDir(output);
        HashMap<String, EQTL> eqtls = new HashMap<String, EQTL>();

        TextFile tf = new TextFile(eqtlfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String gene = elems[4].split("\\.")[0];
            String eqtl = elems[2] + ":" + elems[3] + ":" + gene;
            EQTL eqtlobj = new EQTL();
            eqtlobj.alleles = elems[8];
            eqtlobj.assessed = elems[9];
            eqtlobj.z = Double.parseDouble(elems[10]);
            eqtls.put(eqtl, eqtlobj);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        String[] gtexfiles = Gpio.getListOfFiles(gtexfolder);
        HashSet<String> tissues = new HashSet<String>();
        for (String s : gtexfiles) {
            // Adipose_Subcutaneous.v8.EUR.allpairs.chr1.parquet.txt.gz
            String tissue = s.split("\\.")[0];
            tissues.add(tissue);
            System.out.println("Tissue found: " + tissue);
        }

        System.out.println(tissues.size() + " total tissues.");

        TextFile summary = new TextFile(output + "-summary.txt", TextFile.W);
        if (GTEXsignificantFolder) {
            summary.writeln("Tissue\tShared\tSharedSameDir\tPercSharedSameDir\tSignificantGenes(<gtexfdrthreshold)");
        } else {
            summary.writeln("Tissue\tShared\tSharedSameDir\tPercSharedSameDir\tSignificantGenes(p<" + pvalthreshold + ")");
        }

        int tctr = 1;
        ProgressBar pb = new ProgressBar(tissues.size());
        for (String tissue : tissues) {
            System.out.println("Processing: " + tissue);
            ArrayList<String> files = new ArrayList<>();

            if (GTEXsignificantFolder) {
                files.add(tissue + ".v8.signif_variant_gene_pairs.txt.gz");
            } else {
                for (int c = 1; c < 23; c++) {
                    files.add(tissue + ".v8.EUR.allpairs.chr" + c + ".parquet.txt.gz");
                }
            }

            // reset other eqtl zscores
            for (String key : eqtls.keySet()) {
                eqtls.get(key).otherZ = Double.NaN;
            }

            AtomicInteger sharedF = new AtomicInteger();
            AtomicInteger sharedSameDirF = new AtomicInteger();
            AtomicInteger significantGenes = new AtomicInteger();
            IntStream.range(0, files.size()).parallel().forEach(f -> {
                String file = gtexfolder + files.get(f);
                try {
                    if (GTEXsignificantFolder) {
                        processFileSignificant(eqtls, file, significantGenes, sharedF, sharedSameDirF);
                    } else {
                        processFile(eqtls, file, significantGenes, sharedF, sharedSameDirF, pvalthreshold);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });

            // write results
            TextFile out = null;
            if (writeTissueOutput) {
                out = new TextFile(output + tissue + ".txt.gz", TextFile.W);
                out.writeln("EQTL\tZ\totherZ");
            }

            int shared = 0;
            int shareddirection = 0;
            for (String key : eqtls.keySet()) {
                EQTL e = eqtls.get(key);
                if (!Double.isNaN(e.otherZ)) {
                    if (writeTissueOutput) {
                        out.writeln(key + "\t" + e.z + "\t" + e.otherZ);
                    }
                    shared++;
                    if ((e.z >= 0 && e.otherZ >= 0) || (e.z < 0 && e.otherZ < 0)) {
                        shareddirection++;
                    }
                }
            }
            if (writeTissueOutput) {
                out.close();
            }

            double perc = (double) shareddirection / shared;
            System.out.println(tctr + "/" + tissues.size() + "\t" + tissue + " Done.\t" + shared + "\t" + shareddirection + "\t" + perc);
            summary.writeln(tissue + "\t" + shared + "\t" + shareddirection + "\t" + perc + "\t" + significantGenes.get());
            summary.flush();
            pb.set(tctr);
            tctr++;

            System.out.println();
        }
        pb.close();
        summary.close();

    }


    private void processFileSignificant(HashMap<String, EQTL> eqtls, String file, AtomicInteger
            significantGenes, AtomicInteger sharedF, AtomicInteger sharedSameDirF) throws IOException {

        TextFile tf = new TextFile(file, TextFile.R, 8 * 1048576);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int lnctr = 0;
        HashSet<String> uniqueGenes = new HashSet<String>();
        while (elems != null) {
            String gene = elems[1].split("\\.")[0];
            String[] snpelems = elems[0].split("_");
            uniqueGenes.add(gene);
            String snp = snpelems[0].replaceAll("chr", "") + ":" + snpelems[1];
            String key = snp + ":" + gene;
            EQTL e = eqtls.get(key);
            if (e != null) {
                String alleles = snpelems[2] + "/" + snpelems[3];
                String assessed = snpelems[3];
                Boolean b = BaseAnnot.flipalleles(e.alleles, e.assessed, alleles, assessed);
                if (b != null) {
                    double slope = Double.parseDouble(elems[7]);
                    double slopese = Double.parseDouble(elems[8]);
                    double z = slope / slopese;
                    if (b) {
                        z *= -1;
                    }
                    e.otherZ = z;
                    sharedF.getAndIncrement();
                    if ((e.z >= 0 && e.otherZ >= 0) || (e.z < 0 && e.otherZ < 0)) {
                        sharedSameDirF.getAndIncrement();
                    }
                }
            }
            lnctr++;
            if (lnctr % 5000000 == 0) {
                System.out.println(file + "\t" + lnctr + " lines parsed\t" + sharedF.get() + " eqtls shared\t" + sharedSameDirF.get() + " same direction.");
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        significantGenes.getAndAdd(uniqueGenes.size());
    }

    private void processFile(HashMap<String, EQTL> eqtls, String file, AtomicInteger
            significantGenes, AtomicInteger sharedF, AtomicInteger sharedSameDirF, double pvalthreshold) throws IOException {


        TextFile tf = new TextFile(file, TextFile.R, 8 * 1048576);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int lnctr = 0;
        HashSet<String> uniqueGenes = new HashSet<String>();
        while (elems != null) {
            String gene = elems[0].split("\\.")[0];
            String[] snpelems = elems[1].split("_");

            String snp = snpelems[0].replaceAll("chr", "") + ":" + snpelems[1];
            String key = snp + ":" + gene;
            EQTL e = eqtls.get(key);
            if (e != null) {
                String alleles = snpelems[2] + "/" + snpelems[3];
                String assessed = snpelems[3];
                Boolean b = BaseAnnot.flipalleles(e.alleles, e.assessed, alleles, assessed);
                if (b != null) {
                    double p = Double.parseDouble(elems[4]);
                    if (p < pvalthreshold) {
                        uniqueGenes.add(gene);
                    }
                    double z = Double.parseDouble(elems[3]);
                    if (b) {
                        z *= -1;
                    }
                    e.otherZ = z;
                    sharedF.getAndIncrement();
                    if ((e.z >= 0 && e.otherZ >= 0) || (e.z < 0 && e.otherZ < 0)) {
                        sharedSameDirF.getAndIncrement();
                    }
                }
            }
            lnctr++;
            if (lnctr % 5000000 == 0) {
                System.out.println(file + "\t" + lnctr + " lines parsed\t" + sharedF.get() + " eqtls shared\t" + sharedSameDirF.get() + " same direction.");
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        significantGenes.getAndAdd(uniqueGenes.size());
    }
}
