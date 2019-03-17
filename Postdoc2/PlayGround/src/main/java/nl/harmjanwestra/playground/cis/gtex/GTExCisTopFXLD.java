package nl.harmjanwestra.playground.cis.gtex;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.legacy.vcf.DetermineLD;
import nl.harmjanwestra.playground.legacy.vcf.VCFTabix;
import nl.harmjanwestra.playground.legacy.vcf.VCFVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.ViolinBoxPlot;
import umcg.genetica.graphics.panels.HistogramPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class GTExCisTopFXLD {


    public static void main(String[] args) {

        String eqtlgenfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
        String eqtlgenfiletop = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-topfx.txt.gz";
        String gtextoploc = "D:\\Sync\\SyncThing\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v7_eQTL\\";
        String tabixprefix = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\chrCHR.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz";
        String samplelimit = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-p35va-europeans.txt";
        String outfileloc = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\ld-eur\\ldwithtopfx.txt";

        String tissuelimit = "Brain_Nucleus_accumbens_basal_ganglia.";

        GTExCisTopFXLD ld = new GTExCisTopFXLD();

        try {
//            ld.gettopfx(eqtlgenfile, eqtlgenfiletop);
            ld.run(gtextoploc, eqtlgenfiletop, tabixprefix, samplelimit, tissuelimit, outfileloc);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }
    }

    private void gettopfx(String eqtlgenfile, String output) throws IOException {
        HashMap<String, EQTL> eqtlgentopfx = new HashMap<>();
        TextFile tf = new TextFile(eqtlgenfile, TextFile.R);
        TextFile tfo = new TextFile(output, TextFile.W);
        tfo.writeln(tf.readLine());
        String[] elems = tf.readLineElems(TextFile.tab);
        System.out.println("Parsing: " + eqtlgenfile);
        int ectr = 0;
        while (elems != null) {
            String gene = elems[4];
            if (!eqtlgentopfx.containsKey(gene)) {
                EQTL e = new EQTL();
                e.setRsName(elems[1]);
                e.setRsChr(Byte.parseByte(elems[2]));
                e.setRsChrPos(Integer.parseInt(elems[3]));
                e.setZscore(Double.parseDouble(elems[10]));
                e.setProbeChrPos(Integer.parseInt(elems[6]));
                eqtlgentopfx.put(gene, e);
                tfo.writeln(Strings.concat(elems, Strings.tab));
            }

            ectr++;
            if (ectr % 1000000 == 0) {
                System.out.println(ectr + " eqtls parsed. " + eqtlgentopfx.size() + " loaded.");
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tfo.close();
    }


    public void run(String gtextoploc,
                    String eqtlgenfile,
                    String tabixprefix,
                    String samplelimit,
                    String tissuelimit,
                    String outfileloc) throws IOException, DocumentException {

        HashMap<String, EQTL> eqtlgentopfx = new HashMap<>();

        TextFile tf = new TextFile(eqtlgenfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        System.out.println("Parsing: " + eqtlgenfile);
        int ectr = 0;
        while (elems != null) {
            String gene = elems[4];
            if (!eqtlgentopfx.containsKey(gene)) {
                EQTL e = new EQTL();
                e.setRsName(elems[1]);
                e.setRsChr(Byte.parseByte(elems[2]));
                e.setRsChrPos(Integer.parseInt(elems[3]));
                e.setZscore(Double.parseDouble(elems[10]));
                e.setProbeChrPos(Integer.parseInt(elems[6]));
                eqtlgentopfx.put(gene, e);
            }

            ectr++;
            if (ectr % 1000000 == 0) {
                System.out.println(ectr + " eqtls parsed. " + eqtlgentopfx.size() + " loaded.");
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(eqtlgentopfx.size() + " top fx loaded.");
        GTExCisRepl c = new GTExCisRepl();
        HashMap<String, HashMap<String, EQTL>> gtextopfx = c.getGTExTopFx(gtextoploc, "v7.egenes.txt.gz", null, true, true);


        ArrayList<String> keys = new ArrayList<>();
        keys.addAll(gtextopfx.keySet());
        Collections.sort(keys);

        DetermineLD ld = new DetermineLD();

        int nrbins = 10;

        TextFile out = new TextFile(outfileloc, TextFile.W);

        String header = "tissue\tnrAll\tnrSig\tmeanAll\tmedianall\tAboveThreshold\tmeanSig\tmedianSig\tAboveThresholdAll";
        for (int b = 0; b < nrbins; b++) {
            header += "\tBinAll" + b;
        }
        for (int b = 0; b < nrbins; b++) {
            header += "\tBinSig" + b;
        }
        out.writeln(header);


        HashMap<SNPFeature, VCFVariant> variantHash = new HashMap<>();

        int tctr = 0;


        int nrcols = 5;
        int remainder = gtextopfx.size() % nrcols;
        int nrrows = (int) Math.floor(gtextopfx.size() / nrcols);
        if (remainder > 0) {
            nrrows++;
        }

        Grid gridsig = new Grid(200, 200, nrrows, nrcols, 75, 100);
        Grid gridall = new Grid(200, 200, nrrows, nrcols, 75, 100);

        int colctr = 0;
        int rowctr = 0;
// write bottom list of eqtls


        for (String tissue : keys) {

            if (tissuelimit == null || tissue.equals(tissuelimit)) {
                String tissue2 = tissue.replaceAll("\\.", "");
                TextFile tfz = new TextFile(outfileloc + "-" + tissue2 + "-LDComparison.txt", TextFile.W);
                String headerld = "Bin\trsq\tGene\tGTExRS\tGTExChrpos\tGTExZ\tGTExDist\tEQTLRS\tEQTLChrpos\tEQTLZ\tEQTLDist";
                tfz.writeln(headerld);
                String nicetissue = tissue.replaceAll("_", " ");
                nicetissue = nicetissue.replaceAll("\\.", "");
                System.out.println("Tissue: " + tissue + " " + tctr + "/" + gtextopfx.size());
                int[] binssig = new int[nrbins];
                int[] binsall = new int[nrbins];
                ArrayList<Double> varsSig = new ArrayList<>();
                ArrayList<Double> varsAll = new ArrayList<>();


                int abovethresholdsig = 0;
                int abovethresholdall = 0;
                HashMap<String, EQTL> tissuetopfx = gtextopfx.get(tissue);
                int gctr = 0;
                int gtestedctr = 0;


                ArrayList<ArrayList<Double>> distancesGTEx = new ArrayList<>();
                ArrayList<ArrayList<Double>> distancesEQTLGen = new ArrayList<>();
                ArrayList<ArrayList<Double>> zsGTEx = new ArrayList<>();
                ArrayList<ArrayList<Double>> zsEQTLGen = new ArrayList<>();
                for (int i = 0; i < nrbins; i++) {
                    distancesGTEx.add(new ArrayList<>());
                    distancesEQTLGen.add(new ArrayList<>());
                    zsGTEx.add(new ArrayList<>());
                    zsEQTLGen.add(new ArrayList<>());
                }


                for (String gene : tissuetopfx.keySet()) {
                    EQTL gtexe = tissuetopfx.get(gene);
                    EQTL eqtlgene = eqtlgentopfx.get(gene);

                    if (gtexe != null && eqtlgene != null) {
                        boolean gtexsig = false;
                        if (gtexe.getPvalue() < gtexe.getPvalueAbs()) {
                            gtexsig = true;
                        }

                        // get both variants
                        String tabixfile = tabixprefix.replaceAll("CHR", "" + gtexe.getRsChr());
                        VCFTabix t = new VCFTabix(tabixfile);
                        boolean[] includesamples = t.getSampleFilter(samplelimit);


                        SNPFeature gtexf = tofeat(gtexe);
                        SNPFeature eqtlf = tofeat(eqtlgene);


                        VCFVariant gtexvar = variantHash.get(gtexf);
                        boolean updatehash = false;
                        if (gtexvar == null) {
                            gtexvar = t.getVariant(tofeat(gtexe), includesamples);
//                        variantHash.put(gtexf, gtexvar);
                            updatehash = true;
                        }

                        VCFVariant eqtlvar = variantHash.get(eqtlf);
                        if (eqtlvar == null) {
                            eqtlvar = t.getVariant(tofeat(eqtlgene), includesamples);
                            variantHash.put(eqtlf, eqtlvar);
                            updatehash = true;

                        }
//                    if (updatehash) {
//                        System.out.println(variantHash.size() + " vars in hash");
//                    }

                        if (gtexvar != null && eqtlvar != null) {
                            Pair<Double, Double> ldvars = ld.getLD(gtexvar, eqtlvar);
                            double rsq = ldvars.getRight();
                            if (!Double.isNaN(rsq)) {
                                gtestedctr++;
                                if (rsq > 0.8) {
                                    abovethresholdall++;
                                }
                                int binno = (int) Math.floor(rsq * nrbins);
                                if (binno > nrbins - 1) {
                                    binno = nrbins - 1;
                                }
//                            if (varsAll.size() > 1000) {
//                                break;
//                            }


                                double genedistEQTLGen = Math.abs(eqtlgene.getProbeChrPos() - eqtlgene.getRsChrPos());
                                double genedistGTEx = Math.abs(eqtlgene.getProbeChrPos() - gtexe.getRsChrPos());

                                if (genedistEQTLGen > 1000000) {
                                    genedistEQTLGen = 1000000;

                                }
                                if (genedistGTEx > 1000000) {
                                    genedistGTEx = 1000000;
                                }


                                distancesEQTLGen.get(binno).add(genedistEQTLGen / 1000000);
                                distancesGTEx.get(binno).add(genedistGTEx / 1000000);

                                if (Math.abs(eqtlgene.getZscore()) > 15) {
                                    zsEQTLGen.get(binno).add(15d);
                                } else {
                                    zsEQTLGen.get(binno).add(Math.abs(eqtlgene.getZscore()));
                                }

                                if (Math.abs(gtexe.getZscore()) > 15) {
                                    zsGTEx.get(binno).add(Math.abs(15d));
                                } else {
                                    zsGTEx.get(binno).add(Math.abs(gtexe.getZscore()));
                                }



                                String lnout = binno + "\t" + rsq + "\t" + gtexe.getProbe() + "\t" + gtexe.getRsName() + "\t" + gtexe.getRsChrPos() + "\t" + Math.abs(gtexe.getZscore()) + "\t" + genedistGTEx
                                        + "\t" + eqtlgene.getRsName() + "\t" + eqtlgene.getRsChrPos() + "\t" + Math.abs(eqtlgene.getZscore()) + "\t" + genedistEQTLGen;
                                tfz.writeln(lnout);

                                varsAll.add(rsq);
                                binsall[binno]++;
                                if (gtexsig) {
                                    varsSig.add(rsq);
                                    if (rsq > 0.8) {
                                        abovethresholdsig++;
                                    }
                                    binssig[binno]++;
                                }
                            }
                        }
                    }
                    gctr++;
                    if (gctr % 1000 == 0) {
                        System.out.print("\r" + gctr + " genes processed out of " + tissuetopfx.size() + "\t" + gtestedctr + " tested. " + variantHash.size() + " vars in hash.");
                    }

                }
                tfz.close();
                System.out.println();


                double[][][] distances = new double[2][nrbins][];
                double[][][] zs = new double[2][nrbins][];
                String[][] xlabels = new String[2][nrbins];
                for (int i = 0; i < nrbins; i++) {
                    distances[0][i] = Primitives.toPrimitiveArr(distancesEQTLGen.get(i));
                    distances[1][i] = Primitives.toPrimitiveArr(distancesGTEx.get(i));
                    zs[0][i] = Primitives.toPrimitiveArr(zsEQTLGen.get(i));
                    zs[1][i] = Primitives.toPrimitiveArr(zsGTEx.get(i));
                    xlabels[0][i] = "" + i;
                    xlabels[1][i] = "" + i;
                }
                String[] datasetnames = new String[]{"eQTLGen", "GTEx"};

                String ylabel = "Distance (bp)";
                ViolinBoxPlot vbp = new ViolinBoxPlot();
                String outfilename = outfileloc + "-" + tissue2 + "-LDComparison-Distance.pdf";
                vbp.draw(distances, datasetnames, xlabels, ylabel, ViolinBoxPlot.Output.PDF, outfilename);

                vbp = new ViolinBoxPlot();
                ylabel = "Absolute Z-score (capped)";
                outfilename = outfileloc + "-" + tissue2 + "-LDComparison-ZScore.pdf";
                vbp.draw(zs, datasetnames, xlabels, ylabel, ViolinBoxPlot.Output.PDF, outfilename);

                // vals[group][dataset][value]


                // TODO: mean LD for all variants that have same direction too.. ??

                double[] varsallarr = Primitives.toPrimitiveArr(varsAll);
                double[] varssigarr = Primitives.toPrimitiveArr(varsSig);

                int maxall = JSci.maths.ArrayMath.max(binsall);
                int maxsig = JSci.maths.ArrayMath.max(binssig);

                double meanall = JSci.maths.ArrayMath.mean(varsallarr);
                double meansig = JSci.maths.ArrayMath.mean(varssigarr);
                double medianall = JSci.maths.ArrayMath.median(varsallarr);
                double mediansig = JSci.maths.ArrayMath.median(varssigarr);


                // add a line or something
                String ln = nicetissue
                        + "\t" + varsAll.size()
                        + "\t" + varsSig.size()
                        + "\t" + meanall
                        + "\t" + medianall
                        + "\t" + abovethresholdall
                        + "\t" + meansig
                        + "\t" + mediansig
                        + "\t" + abovethresholdsig
                        + "\t" + Strings.concat(binsall, Strings.tab)
                        + "\t" + Strings.concat(binssig, Strings.tab);

                out.writeln(ln);
                System.out.println(varsAll.size() + "all\t" + meanall + " mean all\t" + varsSig.size() + " sig\t" + meansig + " mean sig\t" + medianall + " median all\t" + mediansig + " median sig");
                Arrays.sort(varsallarr);
                tctr++;


                DecimalFormat format = new DecimalFormat("#.#");
                DecimalFormat usdf = (DecimalFormat) DecimalFormat.getNumberInstance(Locale.US);

                HistogramPanel p = new HistogramPanel(1, 1);
                p.setAxisLabels("LD between top variants", "Frequency");
                p.setData(binssig);
                p.setTitle(nicetissue + " (n=" + usdf.format(varsSig.size()) + ")");
                p.setRange(new Range(0, 0, 1, maxsig));
                String percabovesig = format.format(((double) abovethresholdsig / varsSig.size()) * 100);
                p.setMarginBetweenBinClusters(0);

                p.setThreshold(0.8, "" + usdf.format(abovethresholdsig) + "\n(" + percabovesig + "%)");
                gridsig.addPanel(p);

                HistogramPanel p2 = new HistogramPanel(1, 1);
                p2.setAxisLabels("LD between top variants", "Frequency");
                Range rangeall = new Range(0, 0, 1, maxall);
                System.out.println("Input range: " + rangeall);
                p2.setRange(rangeall);


                p2.setData(binsall);
                p2.setTitle(nicetissue + " (n=" + usdf.format(varsAll.size()) + ")");

                String percaboveall = format.format(((double) abovethresholdall / varsAll.size()) * 100);
                p2.setThreshold(0.8, "" + usdf.format(abovethresholdall) + "\n(" + percaboveall + "%)");
                p2.setMarginBetweenBinClusters(0);
                gridall.addPanel(p2);


//                for (int i = 0; i < bottomLDEqtlgen.size(); i++) {
//
//                }
//                for (int i = 0; i < topLDEqtlgen.size(); i++) {
//                    EQTL egt = topLDGTEx.get(i);
//                    EQTL eqe = topLDEqtlgen.get(i);
//                    int deltaEQE = Math.abs(eqe.getProbeChrPos() - eqe.getRsChrPos());
//                    int deltaEGT = Math.abs(eqe.getProbeChrPos() - egt.getRsChrPos());
//                    String lnout = (nrbins - 1) + "\t" + egt.getProbe() + "\t" + egt.getRsName() + "\t" + egt.getRsChrPos() + "\t" + Math.abs(egt.getZscore()) + "\t" + deltaEGT
//                            + "\t" + eqe.getRsName() + "\t" + eqe.getRsChrPos() + "\t" + Math.abs(eqe.getZscore()) + "\t" + deltaEQE;
//                    tfz.writeln(lnout);
//                }
                tfz.close();

            }


        }
        out.close();

        gridall.draw(outfileloc + "gridAll.pdf");
        gridsig.draw(outfileloc + "gridSig.pdf");

    }

    private SNPFeature tofeat(EQTL gtexe) {

        SNPFeature f = new SNPFeature();

        f.setChromosome(Chromosome.parseChr("" + gtexe.getRsChr()));
        f.setStart(gtexe.getRsChrPos());
        f.setStop(gtexe.getRsChrPos());

        return f;

    }

}
