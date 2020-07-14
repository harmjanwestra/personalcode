package nl.harmjanwestra.playground.biogen.freeze2dot1;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class ReplaceIDs {


    public static void main(String[] args) {
        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\eQTLGen-32k-Blood-significantSNPProbeCombos-trans-b38.txt.gz";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\eQTLGen-32k-Blood-significantSNPProbeCombos-trans-b38-MetaBrainFreeze2dot1IDs.txt.gz";
        int genecol = 1;
        int snpcol = 0;
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        String filelistfile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-06-02-SNPIds\\lsof.txt";

        input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene-combos.txt";
        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene-combos-MetaBrainFreeze2dot1IDs.txt";

        input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene.txt.gz";
        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-01-31-eqtlgen-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38-topfxpergene-MetaBrainFreeze2dot1IDs.txt.gz";
        genecol = 4;
        snpcol = 1;

        input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-b38-MetaBrainFreeze2dot1IDs.txt.gz";


        ReplaceIDs r = new ReplaceIDs();
        try {
            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Freeze2CisEQTLs-combos.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Freeze2CisEQTLs-combos-MetaBrain2dot1IDs.txt";
            snpcol = 0;
            genecol = 1;
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);
            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Freeze2TransEQTLs-combos.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Freeze2TransEQTLs-combos-MetaBrain2dot1IDs.txt";
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\Freeze2CisEQTLs-MetaBrain2dot1IDs.txt.gz";
            snpcol = 1;
            genecol = 4;
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);
            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Trans\\eQTLsFDR0.05.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\Freeze2TransEQTLs-MetaBrain2dot1IDs.txt.gz";
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);

            String genelist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-18-eqtls\\freeze2comp\\trans\\oppositeGenes.txt";
//            r.geneToGeneSymbol(genelist, gtf, null);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-PatchSequenceIssue\\MAPT\\MAPT-metabrain-Cis-Cortex-EUR.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-09-24-PatchSequenceIssue\\MAPT\\MAPT-metabrain-Cis-Cortex-EUR-metaBrainFreeze2dot1Ids.txt";
            snpcol = 1;
            genecol = 4;
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-MetaBrain2IDs.txt.gz";
            snpcol = 1;
            genecol = 4;
//            r.run(filelistfile, gtf, input, output, snpcol, genecol);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added-metabrain2dot1ids.txt.gz";
//            r.replaceSNPIds(input, output, filelistfile, 0);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\binary\\ZScoreMatrix.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\binary\\ZScoreMatrix-MetaBrain2dot1IDs.txt.gz";
//            r.replaceSNPIds(input, output, filelistfile, 0);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2019-12-11-eQTLProbesFDR0.05-ProbeLevel.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2019-12-11-cis-eQTLProbesFDR0.05-ProbeLevel-MetaBrain2dot1IDs.txt.gz";
//            r.run(filelistfile, gtf, input, output, 1, 4);
            String outputUpdatedPos = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\eqtlgen\\2019-12-11-cis-eQTLProbesFDR0.05-ProbeLevel-MetaBrain2dot1IDsandPos.txt.gz";
            String geneAnnotation = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz";
            String snpmap = "";
//            r.updatePositions(output, snpmap, geneAnnotation, outputUpdatedPos);

            // D:\Sync\SyncThing\Postdoc2\2019-BioGen\data\2020-01-Freeze2dot1\2020-05-26-assoc\cisBIOS\010517-BIOS-independent-eQTLs.txt
            // rewrite bios ids
            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-independent-eQTLs.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-independent-eQTLs-MetaBrain2dot1IDs.txt";
            // r.run(filelistfile, gtf, input, output, 0, 1);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-testedGenes.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\cisBIOS\\010517-BIOS-testedGenes-MetaBrain2dot1IDs.txt";
//            r.run(filelistfile, gtf, input, output, -1, 0);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2019-12-11-eQTLProbesFDR0.05-ProbeLevel.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\eqtlgen\\2019-12-11-eqtlgen-cis-eQTLProbesFDR0.05-ProbeLevel-MetaBrain2dot1IDs.txt.gz";
//            r.run(filelistfile, gtf, input, output, 1, 4);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\eqtlgen\\2019-12-11-eqtlgen-cis-eQTLsFDR-ProbeLevel-MetaBrain2dot1IDs.txt.gz";
//            r.run(filelistfile, gtf, input, output, 1, 4);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-eQTLMeta\\data\\2018-05-22-assoc\\trans\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt.gz";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-26-assoc\\eqtlgen\\2019-12-11-eqtlgen-trans-eQTLsFDR-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-MetaBrain2dot1IDs.txt.gz";
//            r.run(filelistfile, gtf, input, output, 1, 4);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2019-12-18-BenteSNPs\\Alzheimer.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2019-12-18-BenteSNPs\\Alzheimer-metabrainfreeze2ids.txt";
            r.run(filelistfile, gtf, input, output, 0, -1);

            input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2019-12-18-BenteSNPs\\Depression.txt";
            output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2019-12-18-BenteSNPs\\Depression-metabrainfreeze2ids.txt";
            r.run(filelistfile, gtf, input, output, 0, -1);

        } catch (IOException e) {


            e.printStackTrace();
        }
    }

    private void updatePositions(String output, String snpmap, String probeannotation, String outputUpdatedPos) throws IOException {


//        HashMap<String, String> ids = new HashMap<String, String>();
//        TextFile tf = new TextFile(output, TextFile.R);
//        tf.readLine();
//        String[] elems = tf.readLineElems(TextFile.tab);
//        while (elems != null) {
//            String id = elems[1];
//            ids.put(id, null);
//            elems = tf.readLineElems(TextFile.tab);
//        }
//        tf.close();

//        tf = new TextFile(snpmap, TextFile.R);
//        elems = tf.readLineElems(TextFile.tab);
//        while (elems != null) {
//            String id = elems[2];
//            if (ids.containsKey(id)) {
//                ids.put(id, elems[0] + ":" + elems[1]);
//            }
//            elems = tf.readLineElems(TextFile.tab);
//        }
//        tf.close();


        HashMap<String, String> genetopos = new HashMap<String, String>();
        TextFile p1 = new TextFile(probeannotation, TextFile.R);
        p1.readLine();
        String[] elems = p1.readLineElems(TextFile.tab);
        while (elems != null) {
            genetopos.put(elems[1], elems[3] + ":" + elems[4]);
            elems = p1.readLineElems(TextFile.tab);
        }
        p1.close();

        TextFile tf = new TextFile(output, TextFile.R);
        TextFile tfout = new TextFile(outputUpdatedPos, TextFile.W);
        tfout.writeln(tf.readLine());
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[1];
            String[] idelems = id.split(":");
            if (idelems.length > 1) {
                elems[2] = idelems[0];
                elems[3] = idelems[1];
            } else {
                elems[2] = "-1";
                elems[3] = "-1";
            }


            String gene = elems[4];
            String genepos = genetopos.get(gene);
            if (genepos == null) {
                elems[5] = "-1";
                elems[6] = "-1";
            } else {
                String[] poselems = genepos.split(":");
                elems[5] = poselems[0];
                elems[6] = poselems[1];
            }

            tfout.writeln(Strings.concat(elems, Strings.tab));

            elems = tf.readLineElems(TextFile.tab);
        }
        tfout.close();

    }

    private void replaceSNPIds(String input, String output, String filelistfile, int snpcol) throws IOException {
        HashMap<String, String> rsToId = new HashMap<String, String>();
        TextFile tf = new TextFile(input, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            String snpid = elems[snpcol];
            rsToId.put(snpid, "NotInDs-" + snpid);
            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 100000 == 0) {
                System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
            }
        }
        tf.close();
        System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
        System.out.println();
        TextFile tf1 = new TextFile(filelistfile, TextFile.R);
        ArrayList<String> list = tf1.readAsArrayList();
        tf1.close();

        for (String s : list) {
            TextFile tf2 = new TextFile(s, TextFile.R);

            String ln = tf2.readLine();

            while (ln != null) {
                elems = ln.split(":");
                if (rsToId.containsKey(ln)) {
                    rsToId.put(ln, ln);
                } else if (rsToId.containsKey(elems[2])) {
                    rsToId.put(elems[2], ln);
                }
                ln = tf2.readLine();
            }
            tf2.close();
            System.out.println(rsToId.size() + " ids after reading: " + s);
        }

        tf = new TextFile(input, TextFile.R);
        TextFile tf2 = new TextFile(output, TextFile.W);
        tf2.writeln(tf.readLine());
        elems = tf.readLineElems(TextFile.tab);
        ctr = 0;
        while (elems != null) {
            String snpid = elems[snpcol];
            elems[snpcol] = rsToId.get(snpid);
            tf2.writeln(Strings.concat(elems, Strings.tab));
            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 100000 == 0) {
                System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
            }
        }
        tf.close();
        tf2.close();
        System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
        System.out.println();
    }

    public void geneToGeneSymbol(String list, String gtf, String out) throws IOException {
        TextFile tf = new TextFile(list, TextFile.R);
        ArrayList<String> arrl = tf.readAsArrayList();
        tf.close();

        HashSet<String> set = new HashSet<String>();
        set.addAll(arrl);

        GTFAnnotation gtfAnnotation = new GTFAnnotation(gtf);
        Collection<Gene> genes = gtfAnnotation.getGenes();
        for (Gene g : genes) {
            String name = g.getName().split("\\.")[0];
            if (set.contains(name)) {
                System.out.println(name + "\t" + g.getGeneSymbol());
            }
        }


    }

    public void run(String filelistfile, String gtf, String input, String output, int snpcol, int genecol) throws IOException {

        HashMap<String, String> rsToId = new HashMap<String, String>();
        if (snpcol > -1) {
            TextFile tf = new TextFile(input, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);
            int ctr = 0;
            while (elems != null) {
                String snpid = elems[snpcol];
                rsToId.put(snpid, "NotInDs-" + snpid);
                elems = tf.readLineElems(TextFile.tab);
                ctr++;
                if (ctr % 100000 == 0) {
                    System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
                }
            }
            tf.close();
            System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
            System.out.println();
            TextFile tf1 = new TextFile(filelistfile, TextFile.R);
            ArrayList<String> list = tf1.readAsArrayList();
            tf1.close();


            for (String s : list) {
                TextFile tf2 = new TextFile(s, TextFile.R);

                String ln = tf2.readLine();

                while (ln != null) {
                    elems = ln.split(":");
                    if (rsToId.containsKey(ln)) {
                        rsToId.put(ln, ln);
                    } else if (rsToId.containsKey(elems[2])) {
                        rsToId.put(elems[2], ln);
                    }
                    ln = tf2.readLine();
                }
                tf2.close();
                System.out.println(rsToId.size() + " ids after reading: " + s);
            }
        }

        GTFAnnotation g = new GTFAnnotation(gtf);
        Collection<Gene> genes = g.getGenes();
        HashMap<String, String> genemap = new HashMap<String, String>();
        for (Gene gene : genes) {
            String id = gene.getName().split("\\.")[0];
            genemap.put(id, gene.getName());
        }

        TextFile tf = new TextFile(input, TextFile.R);
        TextFile out = new TextFile(output, TextFile.W);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            if (snpcol > -1) {
                String snpid = elems[snpcol];
                elems[snpcol] = rsToId.get(snpid);
            }
            if (genecol > -1) {
                String geneid = elems[genecol].split("\\.")[0];
                String replacements = genemap.get(geneid);
                if (replacements == null) {
                    replacements = "NotInDs-" + geneid;
                }
                elems[genecol] = replacements;
            }
            out.writeln(Strings.concat(elems, Strings.tab));
            elems = tf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 100000 == 0) {
                System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
            }
        }
        System.out.print(ctr + " lines processed. " + rsToId.size() + " rsids found sofar.\r");
        System.out.println();
        tf.close();
        out.close();

    }
}
