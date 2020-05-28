package nl.harmjanwestra.playground.biogen.freeze2dot1.ddg2p;

import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class RewriteGeneIDs {

    public static void main(String[] args) {
        String gtf = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz";
        String file = "U:\\2020-05-DDG2P\\DDG2P_11_5_2020.csv.gz";
        String fileout = "U:\\2020-05-DDG2P\\DDG2P_11_5_2020-ensg.txt";
        RewriteGeneIDs r = new RewriteGeneIDs();

        String synonyms = "U:\\2020-05-DDG2P\\Homo_sapiens.gene_info.gz";
        try {
            r.rewrite(gtf, synonyms, file, fileout);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void rewrite(String gtffile, String synonymfile, String file, String out) throws IOException {
        GTFAnnotation gtf = new GTFAnnotation(gtffile);

        HashMap<String, String> symbolToENSG = new HashMap<String, String>();
        Collection<Gene> genes = gtf.getGenes();
        for (Gene g : genes) {
            symbolToENSG.put(g.getGeneSymbol(), g.getName());
        }


        HashMap<String, ArrayList<String>> synonyms = new HashMap<String, ArrayList<String>>();
        TextFile tf2 = new TextFile(synonymfile, TextFile.R);
        String[] elems2 = tf2.readLineElems(TextFile.tab);
        while (elems2 != null) {
            String hygo = elems2[2];

            String syns = elems2[4];
            String[] synelems = Strings.pipe.split(syns);
            ArrayList<String> allSyns = new ArrayList<>();
            for (String s : synelems) {
                allSyns.add(s);
            }
            synonyms.put(hygo, allSyns);
            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();

        TextFile tf = new TextFile(file, TextFile.R);
        TextFile tfo = new TextFile(out, TextFile.W);

        int ctr = 0;
        int found = 0;
        String[] header = tf.readLineElems(TextFile.comma);
        tfo.writeln(Strings.concat(header, Strings.tab) + "\tENSG");
        String[] elems = tf.readLineElems(TextFile.comma);
        while (elems != null) {
            String id = elems[0];
            if(id.equals("KIFBP")){
                System.out.println("!");
            }
            String ensg = symbolToENSG.get(id);

            if (ensg == null) {
                // try synonym
                ArrayList<String> synonym = synonyms.get(id);
                if (synonym != null) {
                    for (String s : synonym) {
                        ensg = symbolToENSG.get(s);
                        if (ensg != null) {
                            break;
                        }
                    }
                }
            }

            if (ensg != null) {
                found++;

            } else {
                System.out.println("Can't find:\t" + id);
            }
            ctr++;
            tfo.writeln(Strings.concat(elems, Strings.tab).replaceAll(" ","_") + "\t" + ensg);
            elems = tf.readLineElems(TextFile.comma);
        }
        tf.close();
        tfo.close();
        System.out.println(found + "/" + ctr);
    }


}
