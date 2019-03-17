package nl.harmjanwestra.playground.cis.gtex;

import org.apache.commons.io.comparator.NameFileComparator;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class GTExCisRepl {

    public static void main(String[] args) {

        if (args.length < 4) {
            System.out.println("Usage: toploc allloc eqtlfile outfile");
            System.exit(-1);
        }
        GTExCisRepl c = new GTExCisRepl();
        try {
            c.run(args[0],
                    args[1],
                    args[2],
                    args[3]);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void run(String gtextoploc, String gtexloc, String eqtlfile, String outfile) throws IOException {
        // load thresholds
        HashMap<String, HashMap<String, EQTL>> topeffects = getGTExTopFx(gtextoploc, "v7.egenes.txt.gz", null, true, true);

        // load GTEx eqtls
        ArrayList<Pair<String, HashMap<String, EQTL>>> alleffects = getGTExFx(gtexloc, "allpairs.txt.gz", null, true, true);

        // now build the table
        String header = "Gene" +
                "\tGene-Chr" +
                "\tGene-Pos" +
                "\tRsId" +
                "\tSNP-Chr" +
                "\tSNP-Pos" +
                "\tAlleles" +
                "\tAlleleAssessed";
        String header0 = "-\t-\t-\t-\t-\t-\t-\t-";

        // Z       RSq     P       FDR
        for (int d = -1; d < alleffects.size(); d++) {
            if (d == -1) {
                header0 += "\teQTLGen\t\t\t";
                header += "\tZ\tRsq\tP\tSignificant";
            } else {
                header0 += "\t" + alleffects.get(d).getLeft() + "\t\t";
                header += "\tZ\tP\tSignificant";
            }

        }
        TextFile out = new TextFile(outfile, TextFile.W);
        out.writeln(header0);
        out.writeln(header);

        TextFile in = new TextFile(eqtlfile, TextFile.R);
        in.readLine();
        String ln = in.readLine();
        while (ln != null) {
            String[] elems = ln.split("\t");
            EQTL reference = EQTL.fromString(elems, "-", Strings.semicolon);

            int n = Descriptives.sum(Primitives.toPrimitiveArr(reference.getDatasetsSamples()));
            double r = ZScores.zToR(reference.getZscore(), n);

            String lnout = reference.getProbe() + "\t" + reference.getProbeChr() + "\t" + reference.getProbeChrPos()
                    + "\t" + reference.getRsName() + "\t" + reference.getRsChr() + "\t" + reference.getRsChrPos()
                    + "\t" + reference.getAlleles() + "\t" + reference.getAlleleAssessed()
                    + "\t" + reference.getZscore() + "\t" + (r * r) + "\t" + reference.getPvalue() + "\t" + (reference.getFDR() < 0.05);


            String query = reference.getProbe() + "_" + reference.getRsChr() + "_" + reference.getRsChrPos();
            for (int d = 0; d < alleffects.size(); d++) {
                Pair<String, HashMap<String, EQTL>> eqtlpair = alleffects.get(d);
                String tissue = eqtlpair.getLeft();
                HashMap<String, EQTL> eqtls = eqtlpair.getRight();
                EQTL test = eqtls.get(query);
                if (test == null) {
                    lnout += "\t-\t-\t-";
                } else {

                    Boolean flip = BaseAnnot.flipalleles(reference.getAlleles(), reference.getAlleleAssessed(), test.getAlleles(), test.getAlleleAssessed());
                    if (flip == null) {
                        lnout += "\t-\t-\t-";
                    } else {
                        double z = test.getZscore();
                        if (flip) {
                            z *= -1;
                        }

                        // get whether gene is significant
                        HashMap<String, EQTL> topfx = topeffects.get(tissue);
                        EQTL topeqtl = topfx.get(reference.getProbe());
                        if (topeqtl == null) {
                            System.err.println("ERROR could not find gene " + reference.getProbe() + " in tissue " + tissue);
                            System.exit(-1);
                        }
                        double threshold = topeqtl.getPvalueAbs();
                        boolean significant = false;
                        if (test.getPvalue() < threshold) {
                            significant = true;
                        }
                        lnout += "\t" + z + "\t" + test.getPvalue() + "\t" + significant;
                    }
                }
            }
            out.writeln(lnout);
            ln = in.readLine();
        }
        in.close();
        out.close();

    }

    private ArrayList<Pair<String, HashMap<String, EQTL>>> getGTExFx(String gtexloc, String filefilter, HashSet<String> genelimit, boolean includeblood, boolean includesexchr) throws IOException {
        ArrayList<Pair<String, HashMap<String, EQTL>>> eqtlspertissue = new ArrayList<>();
        File folder = new File(gtexloc);
        System.out.println(gtexloc + "\texists: " + folder.exists());
        if (!folder.exists()) {
            System.exit(-1);
        }
        File[] listOfFiles = folder.listFiles();

        // sort files
        Arrays.sort(listOfFiles, NameFileComparator.NAME_COMPARATOR);


        for (File file : listOfFiles) {
            if (file.getName().endsWith(filefilter)) {
                TextFile br = new TextFile(file, TextFile.R); // new BufferedReader(new InputStreamReader(tarInput)); // Read directly from tarInput
                br.readLine();
//				System.out.println("For File = " + currentEntry.getName());
                String tissue = file.getName();

                if (includeblood || (!includeblood && !tissue.contains("Whole_Blood.v7.egenes"))) {
                    String tissuename = tissue.replaceAll(filefilter, "");
                    HashMap<String, EQTL> topeqtlgencispergenegtextissue = new HashMap<>();
                    String line;
                    while ((line = br.readLine()) != null) {
                        String[] elems = line.split("\t");
                        String geneid = (elems[0].split("\\.")[0]);
                        if (genelimit == null || genelimit.contains(geneid)) {
                            // check whether this is the top fx
							/*
							0 gene_id
							1 variant_id
							2 tss_distance
							3 ma_samples
							4 ma_count
							5 maf
							6 pval_nominal
							7 slope
							8 slope_se
							 */
                            try {
                                double beta = Double.parseDouble(elems[7]);
                                double se = Double.parseDouble(elems[8]);
                                double zscore = beta / se;
                                String[] variantElems = elems[1].split("_"); // 1_754105_C_T_b37
                                Chromosome chr = Chromosome.parseChr(variantElems[0]);
                                if (chr.isAutosome() || (!chr.isAutosome() && includesexchr)) {
                                    EQTL e = new EQTL();
                                    e.setRsChr((byte) chr.getNumber());
                                    e.setRsChrPos(Integer.parseInt(variantElems[1]));
                                    e.setProbe(geneid);
                                    e.setAlleles((variantElems[2] + "/" + variantElems[3]));
                                    e.setAlleleAssessed((variantElems[3]));
                                    e.setZscore(zscore);

                                    Double pval = Double.parseDouble(elems[6]);
                                    e.setPvalue(pval);
                                    String eqtlid = geneid + "_" + variantElems[0] + "_" + variantElems[1];
                                    topeqtlgencispergenegtextissue.put(eqtlid, e);
                                }
                            } catch (NumberFormatException e) {
                                System.out.println("Issue with beta or SE: " + line);
                            }
                        }
                    }
                    eqtlspertissue.add(new Pair<>(tissuename, topeqtlgencispergenegtextissue));
                    System.out.println(tissuename + "\t" + topeqtlgencispergenegtextissue.size());
                }


            }
        }
        return eqtlspertissue;
    }

    public HashMap<String, HashMap<String, EQTL>> getGTExTopFx(String gtexloc,
                                                               String filefilter,
                                                               Set<String> genelimit,
                                                               boolean includesexchr,
                                                               boolean includeblood) throws IOException {

        HashMap<String, HashMap<String, EQTL>> eqtlspertissue = new HashMap<>();

        // now iterate the GTEx file
//		TarArchiveInputStream tarInput = new TarArchiveInputStream(new GzipCompressorInputStream(new FileInputStream(gtextarball)));
//		TarArchiveEntry currentEntry = tarInput.getNextTarEntry();
//		BufferedReader br = null;
//		StringBuilder sb = new StringBuilder();
        File folder = new File(gtexloc);
        System.out.println(gtexloc + "\texists: " + folder.exists());
        if (!folder.exists()) {
            System.exit(-1);
        }
        File[] listOfFiles = folder.listFiles();

        // sort files
        Arrays.sort(listOfFiles, NameFileComparator.NAME_COMPARATOR);

        for (File file : listOfFiles) {
            if (file.getName().endsWith(filefilter)) {
                TextFile br = new TextFile(file, TextFile.R); // new BufferedReader(new InputStreamReader(tarInput)); // Read directly from tarInput
                br.readLine();
//				System.out.println("For File = " + currentEntry.getName());
                String tissue = file.getName();
                if (includeblood || (!includeblood && !tissue.contains("Whole_Blood.v7.egenes"))) {
                    HashMap<String, EQTL> topeqtlgencispergenegtextissue = new HashMap<>();
                    String line;
                    String tissuename = tissue.replaceAll(filefilter, "");
                    while ((line = br.readLine()) != null) {
                        String[] elems = line.split("\t");
					/*
0	gene_id
2	gene_chr
3	gene_start
4	gene_end
10	variant_id
12	chr
13	snp_pos
14	ref
15	alt
16	rs_id_dbSNP142_GRCh37p13
22	pval_nominal
23	slope
24	slope_se
27	qval
28	pval_nominal_threshold
					 */
/* v7:
0	gene_id
1	gene_name
2	gene_chr
3	gene_start
4	gene_end
5	strand
6	num_var
7	beta_shape1
8	beta_shape2
9	true_df
10	pval_true_df
11	variant_id
12	tss_distance
13	chr
14	pos
15	ref
16	alt
17	num_alt_per_site
18	rs_id_dbSNP147_GRCh37p13
19	minor_allele_samples
20	minor_allele_count
21	maf
22	ref_factor
23	pval_nominal
24	slope
25	slope_se
26	pval_perm
27	pval_beta
28	qval
29	pval_nominal_threshold
30	log2_aFC
31	log2_aFC_lower
32	log2_aFC_upper

*/
                        String geneid = (elems[0].split("\\.")[0]);
                        if (genelimit == null || genelimit.contains(geneid)) {
                            // check whether this is the top fx

                            EQTL top = topeqtlgencispergenegtextissue.get(geneid);
                            double beta = Double.parseDouble(elems[24]);
                            double se = Double.parseDouble(elems[25]);
                            double zscore = beta / se;
                            if (top == null) {
                                // if not use this as the top

                                Chromosome chr = Chromosome.parseChr(elems[13]);
                                Chromosome chr2 = Chromosome.parseChr(elems[2]);
                                if (chr.isAutosome() || (!chr.isAutosome() && includesexchr)) {
                                    if (chr2.isAutosome() || (!chr2.isAutosome() && includesexchr)) {
                                        EQTL e = new EQTL();
                                        e.setRsName((elems[18]));
                                        e.setRsChr((byte) chr.getNumber());
                                        e.setRsChrPos(Integer.parseInt(elems[14]));
                                        e.setProbe(geneid);
                                        e.setFDR(Double.parseDouble(elems[28]));
                                        e.setProbeChr((byte) chr2.getNumber());
                                        e.setAlleles((elems[15] + "/" + elems[16]));
                                        e.setAlleleAssessed((elems[16]));
                                        e.setZscore(zscore);

                                        Double pval = Double.parseDouble(elems[23]);
                                        e.setPvalue(pval);
                                        Double pvalthreshold = Double.parseDouble(elems[29]);
                                        e.setPvalueAbs(pvalthreshold);
                                        topeqtlgencispergenegtextissue.put(geneid, e);
                                    }
                                }

                            } else {
                                // compare absolute z-scores
                                double zref = Math.abs(top.getZscore());
                                double ztest = Math.abs(zscore);
                                if (zref < ztest) {
                                    Chromosome chr = Chromosome.parseChr(elems[13]);
                                    Chromosome chr2 = Chromosome.parseChr(elems[2]);
                                    if (chr.isAutosome() || (!chr.isAutosome() && includesexchr)) {
                                        if (chr2.isAutosome() || (!chr2.isAutosome() && includesexchr)) {
                                            EQTL e = new EQTL();
                                            e.setRsName((elems[18]));
                                            e.setRsChr((byte) chr.getNumber());
                                            e.setRsChrPos(Integer.parseInt(elems[14]));
                                            e.setProbe(geneid);
                                            e.setFDR(Double.parseDouble(elems[28]));
                                            e.setProbeChr((byte) chr2.getNumber());
                                            e.setAlleles((elems[15] + "/" + elems[16]));
                                            e.setAlleleAssessed((elems[16]));
                                            e.setZscore(zscore);

                                            Double pval = Double.parseDouble(elems[23]);
                                            e.setPvalue(pval);
                                            Double pvalthreshold = Double.parseDouble(elems[29]);
                                            e.setPvalueAbs(pvalthreshold);
                                            topeqtlgencispergenegtextissue.put(geneid, e);

                                            topeqtlgencispergenegtextissue.put(geneid, e);
                                        }
                                    }
                                }
                            }
                        }
                    }
//				br.close();

                    eqtlspertissue.put(tissuename, topeqtlgencispergenegtextissue);

                    System.out.println(tissuename + "\t" + topeqtlgencispergenegtextissue.size());

//                    if (eqtlspertissue.size() == 1) {
//                        break;
//                    }
                }

            }
        }
        return eqtlspertissue;
    }
}
