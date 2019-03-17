package nl.harmjanwestra.playground.cis;


import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class LiftOverEQTLFile {

    public static void main(String[] args) {

        String efile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.txt.gz";
        String bedout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved.bed";

        LiftOverEQTLFile t = new LiftOverEQTLFile();

        try {
//            t.toBed(efile, bedout);

            String liftfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38.bed";
            String efileout = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-b38.txt.gz";
            String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\annotation\\gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz";
            t.toEQTL(liftfile, efile, efileout, gtf);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void toBed(String efilein, String bedout) throws IOException {
        TextFile tf = new TextFile(efilein, TextFile.R);
        TextFile bedoutf = new TextFile(bedout, TextFile.W);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String rs = elems[1];
            String chr = elems[2];
            Integer pos = Integer.parseInt(elems[3]);
            bedoutf.writeln("chr" + chr + "\t" + pos + "\t" + (pos + 1) + "\t" + rs);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        bedoutf.close();
    }

    public void toEQTL(String liftfile, String efilein, String efileout, String gtfannot) throws IOException {
        HashMap<String, Pair<String, String>> rsToPos = new HashMap<String, Pair<String, String>>();
        TextFile tf = new TextFile(liftfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String rs = elems[3];
            String chr = "" + Chromosome.parseChr(elems[0]);
            String pos = elems[1];
            rsToPos.put(rs, new Pair<>(chr, pos));
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        GTFAnnotation gtf = new GTFAnnotation(gtfannot);
        TextFile ef = new TextFile(efilein, TextFile.R);
        TextFile efo = new TextFile(efileout, TextFile.W);
        HashMap<String, Gene> genemap = new HashMap<>();
        efo.writeln(ef.readLine());
        HashSet<String> genesnotfound = new HashSet<>();
        elems = ef.readLineElems(TextFile.tab);
        while (elems != null) {
            String rs = elems[1];
            String chr = elems[2];

            Pair<String, String> newpos = rsToPos.get(rs);
            if (newpos != null) {
                String gene = elems[4];
                Gene replacementGene = genemap.get(gene);
                if (replacementGene == null && !genesnotfound.contains(gene)) {

                    for (Gene g : gtf.getGenes()) {
                        String name = g.getName();
                        if (name.contains(gene)) {
                            genemap.put(gene, g);

                            replacementGene = g;
                            break;
                        }
                    }

                    if (replacementGene == null) {
                        genesnotfound.add(gene);
                        System.out.println(gene + " not found in GTF");
                    }
                }
                if (replacementGene != null) {
                    Chromosome chrnew = Chromosome.parseChr(newpos.getLeft());
                    if (chrnew.equals(replacementGene.getChromosome())) {
                        Integer pos = Integer.parseInt(newpos.getRight());
                        int gmidpos = (replacementGene.getStop() + replacementGene.getStart()) / 2;
                        if (Math.abs(pos - gmidpos) > 1000000) {
                            System.out.println("New build: " + rs + "\t" + elems[4] + "\t" + replacementGene.getName() + "\t" + pos + "\t" + gmidpos + "\t" + Math.abs(pos - gmidpos));
                        } else {
                            elems[2] = "" + Chromosome.parseChr(newpos.getLeft()).getNumber();
                            elems[3] = newpos.getRight();
                            elems[4] = replacementGene.getName();
                            efo.writeln(Strings.concat(elems, Strings.tab));
                        }
                    } else {
                        System.out.println("Gene " + elems[4] + "\t" + replacementGene.getName() + " are not on same chr: " + chr + "\t" + replacementGene.getChromosome());
                    }
                } else {

                }
            } else {
                System.out.println("SNP " + rs + " not lifted over.");
            }


            elems = ef.readLineElems(TextFile.tab);
        }
        ef.close();
        efo.close();
    }
}
