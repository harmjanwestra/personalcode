package nl.harmjanwestra.playground.cis;


import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.Gene;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.AssociationPanel;
import umcg.genetica.graphics.panels.GenePanel;
import umcg.genetica.graphics.panels.LDPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;

public class EQTLRegionPlot {


    public void run(String eqtlfile, String ldmatrix, String gene, String geneAnnotation) throws Exception {


        ArrayList<EQTL> eqtls = new ArrayList<EQTL>();
        EQTL topeffect = null;
        Feature region = null;
        {
            QTLTextFile tf = new QTLTextFile(eqtlfile, TextFile.R);
            Iterator<EQTL> eqtlit = tf.getEQtlIterator();


            int minpos = Integer.MAX_VALUE;
            int maxPos = 0;
            while (eqtlit.hasNext()) {
                EQTL e = eqtlit.next();
                if (e.getProbe().equals(gene)) {
                    if (topeffect == null || Math.abs(e.getZscore()) > Math.abs(topeffect.getZscore())) {
                        topeffect = e;
                    }
                    if (e.getRsChrPos() > maxPos) {
                        maxPos = e.getRsChrPos();
                    }
                    if (e.getRsChrPos() < minpos) {
                        minpos = e.getRsChrPos();
                    }
                    eqtls.add(e);
                }
            }
            tf.close();
            Chromosome chr = Chromosome.parseChr("" + topeffect.getRsChr());
            region = new Feature(chr, minpos, maxPos);
        }

        DoubleMatrixDataset<String, String> ds = DoubleMatrixDataset.loadDoubleData(ldmatrix);
        ArrayList<Pair<Integer, Integer>> ldpos = new ArrayList<>();
        ArrayList<Double> ldvals = new ArrayList<>();

        for (int i = 0; i < eqtls.size(); i++) {
            Integer id1 = ds.getHashRows().get(eqtls.get(i).getRsName());
            if (id1 != null) {
                for (int j = i + 1; j < eqtls.size(); j++) {
                    Integer id2 = ds.getHashRows().get(eqtls.get(j).getRsName());
                    if (id2 != null) {
                        ldpos.add(new Pair<Integer, Integer>(eqtls.get(j).getRsChrPos(), eqtls.get(j).getRsChrPos()));
                        ldvals.add(ds.getMatrix().getQuick(id1, id2));
                    }
                }
            }
        }




        Grid grid = new Grid(800, 300, 3, 1, 100, 100);
        GenePanel gp = new GenePanel(1, 1);
        {
            GTFAnnotation gtf = new GTFAnnotation(geneAnnotation);

            Collection<Gene> allgenes = gtf.getGenes();
            ArrayList<Gene> genesubset = new ArrayList<>();
            for (Gene g : allgenes) {
                if (g.overlaps(region)) {
                    genesubset.add(g);
                }
            }
            gp.setData(region, genesubset);
            grid.addPanel(gp);
        }


        Integer topFxId = ds.getHashRows().get(topeffect.getRsName());

        double[] lddata = new double[eqtls.size()];
        ArrayList<Pair<Integer, Double>> pvals = new ArrayList<Pair<Integer, Double>>();
        int ctr = 0;
        for (EQTL e : eqtls) {
            pvals.add(new Pair<Integer, Double>(e.getRsChrPos(), -Math.log10(e.getPvalue())));

            if (topFxId != null) {
                Integer id = ds.getHashRows().get(e.getRsName());
                if (id != null) {
                    lddata[ctr] = ds.getMatrix().getQuick(topFxId, id);
                }
                ctr++;
            }
        }


        AssociationPanel assocP = new AssociationPanel(1, 1);
        assocP.setDataSingleDs(region, null, pvals, "eQTLGen");
        assocP.setLDData(lddata);
        grid.addPanel(assocP);

        LDPanel ldp = new LDPanel(1, 1);
        ldp.setData(region, ldpos, ldvals);
        grid.addPanel(ldp);


    }

}
