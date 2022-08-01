package nl.harmjanwestra.playground.qtltest;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class GeneAnnotation {

    HashMap<String, Integer> geneToId = new HashMap<>();
    int[] positions;
    int[] chromosomes;
    String[] hgnc;

    public GeneAnnotation(String file) throws IOException {
        load(file, -1, null);
    }

    public GeneAnnotation(String file, int chr) throws IOException {
        load(file, chr, null);
    }

    public GeneAnnotation(String file, int chr, Set<String> geneLimitSet) throws IOException {
        load(file, chr, geneLimitSet);
    }

    private void load(String file, int chrLimit, Set<String> geneLimitSet) throws IOException {
        System.out.println("Loading gene annotation from: " + file);
        System.out.println("Limiting to genes on chr: " + chrLimit);
        System.out.println("Max number of genes: " + geneLimitSet.size());
        TextFile tf = new TextFile(file, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        ArrayList<String> geneids = new ArrayList<>();
        ArrayList<Integer> chrs = new ArrayList<>();
        ArrayList<Integer> poss = new ArrayList<>();

        int ctr = 0;
        while (elems != null) {
            String ensg = Strings.cache(elems[1]);
            int chr = Chromosome.parseChr(elems[3]).getNumber();
            if ((chrLimit == -1 || chrLimit == chr) && (geneLimitSet == null || geneLimitSet.contains(ensg))) {
                Integer pos = Integer.parseInt(elems[4]);
                String geneid = Strings.cache(elems[2]);
                geneids.add(geneid);
                chrs.add(chr);
                poss.add(pos);
                geneToId.put(ensg, ctr);
                ctr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        positions = Primitives.toPrimitiveArr(poss);
        chromosomes = Primitives.toPrimitiveArr(chrs);
        hgnc = geneids.toArray(new String[0]);
        System.out.println("Annotations loaded for " + geneToId.size() + " genes.");
    }

    public void reorder(ArrayList<String> ids) {
        int[] newpos = new int[ids.size()];
        int[] newchr = new int[ids.size()];
        String[] newhgnc = new String[ids.size()];
        HashMap<String, Integer> newmap = new HashMap<>();
        for (int i = 0; i < ids.size(); i++) {
            String ensg = ids.get(i);
            Integer id = geneToId.get(ensg);
            if (id == null) {
                newpos[i] = -1;
                newchr[i] = -1;
                newhgnc[i] = Strings.cache("-");
            } else {
                newpos[i] = positions[id];
                newchr[i] = chromosomes[id];
                newhgnc[i] = hgnc[id];
            }
            newmap.put(ensg, i);
        }
        positions = newpos;
        chromosomes = newchr;
        hgnc = newhgnc;
        geneToId = newmap;
    }

    public Integer getGeneId(String gene) {
        return geneToId.get(gene);
    }

    public int getChr(int id) {
        return chromosomes[id];
    }

    public int getPos(int id) {
        return positions[id];
    }

    public Set<String> getAllGenes() {
        return geneToId.keySet();
    }
}
