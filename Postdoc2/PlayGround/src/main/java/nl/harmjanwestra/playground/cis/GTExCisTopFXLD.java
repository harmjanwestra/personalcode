package nl.harmjanwestra.playground.cis;

import nl.harmjanwestra.utilities.enums.Chromosome;
import nl.harmjanwestra.utilities.legacy.genetica.containers.Pair;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.features.*;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;
import nl.harmjanwestra.utilities.math.DetermineLD;
import nl.harmjanwestra.utilities.vcf.VCFTabix;
import nl.harmjanwestra.utilities.vcf.VCFVariant;
import org.jfree.util.ArrayUtils;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class GTExCisTopFXLD {


    public static void main(String[] args) {

        String eqtlgenfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
        String gtextoploc = "D:\\Sync\\SyncThing\\Data\\eQTLs\\GTEx\\GTEx_Analysis_v7_eQTL\\";
        String tabixprefix = "D:\\Sync\\SyncThing\\Data\\Ref\\1kg-eqtlgen\\eur\\chrCHR.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.gz";
        String samplelimit = null;
        String outfileloc = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-06-05-cis-gtexrepl\\ld\\ldwithtopfx.txt";


        GTExCisTopFXLD ld = new GTExCisTopFXLD();

        try {
            ld.run(gtextoploc, eqtlgenfile, tabixprefix, samplelimit, outfileloc);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void run(String gtextoploc, String eqtlgenfile, String tabixprefix, String samplelimit, String outfileloc) throws IOException {

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

        DetermineLD ld = new DetermineLD();

        int nrbins = 10;

        TextFile out = new TextFile(outfileloc, TextFile.W);

        String header = "tissue\tnrAll\tnrSig\tmeanAll\tmedianall\tmeanSig\tmedianSig";
        for (int b = 0; b < nrbins; b++) {
            header += "\tBinAll" + b;
        }
        for (int b = 0; b < nrbins; b++) {
            header += "\tBinSig" + b;
        }
        out.writeln(header);


        HashMap<SNPFeature, VCFVariant> variantHash = new HashMap<>();

        int tctr = 0;
        for (String tissue : gtextopfx.keySet()) {
            System.out.println("Tissue: " + tissue + " " + tctr + "/" + gtextopfx.size());
            int[] binssig = new int[nrbins];
            int[] binsall = new int[nrbins];
            ArrayList<Double> varsSig = new ArrayList<>();
            ArrayList<Double> varsAll = new ArrayList<>();


            HashMap<String, EQTL> tissuetopfx = gtextopfx.get(tissue);
            int gctr = 0;
            int gtestedctr = 0;
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
                            int binno = (int) Math.floor(rsq * nrbins);
                            if (binno > nrbins - 1) {
                                binno = nrbins - 1;
                            }

                            varsAll.add(rsq);
                            binsall[binno]++;
                            if (gtexsig) {
                                varsSig.add(rsq);
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
            System.out.println();
            // TODO: mean LD for all variants that have same direction too.. ??
            double[] varsallarr = Primitives.toPrimitiveArr(varsAll);
            double[] varssigarr = Primitives.toPrimitiveArr(varsSig);


            double meanall = JSci.maths.ArrayMath.mean(varsallarr);
            double meansig = JSci.maths.ArrayMath.mean(varssigarr);
            double medianall = JSci.maths.ArrayMath.median(varsallarr);
            double mediansig = JSci.maths.ArrayMath.median(varssigarr);


            // add a line or something
            String ln = tissue
                    + "\t" + varsAll.size()
                    + "\t" + varsSig.size()
                    + "\t" + meanall
                    + "\t" + medianall
                    + "\t" + meansig
                    + "\t" + mediansig
                    + "\t" + Strings.concat(binsall, Strings.tab)
                    + "\t" + Strings.concat(binssig, Strings.tab);

            out.writeln(ln);
            System.out.println(varsAll.size() + "all\t" + meanall + " mean all\t" + varsSig.size() + " sig\t" + meansig + " mean sig\t" + medianall + " median all\t" + mediansig + " median sig");
            Arrays.sort(varsallarr);
            tctr++;
        }
        out.close();

    }

    private SNPFeature tofeat(EQTL gtexe) {

        SNPFeature f = new SNPFeature();

        f.setChromosome(Chromosome.parseChr("" + gtexe.getRsChr()));
        f.setStart(gtexe.getRsChrPos());
        f.setStop(gtexe.getRsChrPos());

        return f;

    }

}
