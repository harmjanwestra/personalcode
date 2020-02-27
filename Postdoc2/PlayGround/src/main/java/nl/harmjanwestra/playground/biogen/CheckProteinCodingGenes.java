package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class CheckProteinCodingGenes {

    public static void main(String[] args) {
        String geneset = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\genesBeforeFiltering.txt";
        String genesettypeout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\genesBeforeFiltering-type.txt";
        String oldbuild = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2019-06-05-protein_coding_genes_gencode24.txt";
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        String afterfilter = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\availablegenes.txt";
        CheckProteinCodingGenes c = new CheckProteinCodingGenes();
        try {
            c.getType(geneset, gtf, oldbuild, afterfilter, genesettypeout);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void getType(String geneset, String gtf, String oldbuild, String afterfilter, String genesettypeout) throws IOException {


        HashSet<String> genes = new HashSet<String>();
        HashSet<String> proteincodingv28genes = new HashSet<String>();
        HashSet<String> genesafterfilter = new HashSet<String>();


        TextFile tf = new TextFile(afterfilter, TextFile.R);
        tf.readLine();
        String line = tf.readLine();
        while (line != null) {
            String[] gene = Strings.dot.split(line);
            genesafterfilter.add(gene[0]);
            line = tf.readLine();
        }
        tf.close();

        tf = new TextFile(oldbuild, TextFile.R);
        tf.readLine();
        line = tf.readLine();
        while (line != null) {
            String[] gene = Strings.dot.split(line);
            proteincodingv28genes.add(gene[0]);
            line = tf.readLine();
        }
        tf.close();

        tf = new TextFile(geneset, TextFile.R);
        tf.readLine();
        line = tf.readLine();
        while (line != null) {
            String[] gene = Strings.dot.split(line);
            genes.add(gene[0]);
            line = tf.readLine();
        }
        tf.close();

        GTFAnnotation g = new GTFAnnotation(gtf);
        ArrayList<Gene> genelist = g.getGenesAsArrayList();
        TextFile outf = new TextFile(genesettypeout, TextFile.W);
        outf.writeln("Gene\tType\tIsProteinCodingInb28\tInFileAfterQC");
        for (Gene gene : genelist) {
            String name = gene.getName();
            String[] geneelems = Strings.dot.split(name);
            if (genes.contains(geneelems[0])) {
                outf.writeln(gene.getName() + "\t" + gene.getType() + "\t" + proteincodingv28genes.contains(geneelems[0]) + "\t" + genesafterfilter.contains(geneelems[0]));
            }
        }
        outf.close();
    }
}
