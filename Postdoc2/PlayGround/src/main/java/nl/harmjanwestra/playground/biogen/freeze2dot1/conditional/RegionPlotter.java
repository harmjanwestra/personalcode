package nl.harmjanwestra.playground.biogen.freeze2dot1.conditional;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.playground.legacy.GTFAnnotation;
import umcg.genetica.containers.Pair;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.Gene;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.AssociationPanel;
import umcg.genetica.graphics.panels.GenePanel;
import umcg.genetica.io.text.TextFile;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class RegionPlotter {

    public static void main(String[] args) {

        if (args.length < 6) {
            System.out.println("Usage: eqtlfiletemplate gtf.gz querygene startiter stopiter output.pdf/output.png [limitRegionTo100kbAroundTopSNP:false] [maxp:300] [makeYAxisEqual:false]\nFor example: IterationITER Gencodev32.gtf.gz ENSG0000001 1 1 output.png true");
        } else {
            RegionPlotter p = new RegionPlotter();
            String efilename = args[0];
            String gtf = args[1];
            String query = args[2];
            Integer start = Integer.parseInt(args[3]);
            Integer stop = Integer.parseInt(args[4]);
            String output = args[5];
            boolean limitRegionTo100kbAroundTopSNP = false;
            if (args.length >= 7) {
                limitRegionTo100kbAroundTopSNP = Boolean.parseBoolean(args[6]);
            }

            Double maxP = null;
            if (args.length >= 8) {
                maxP = Double.parseDouble(args[7]);
            }

            boolean makeYAxisEqual = false;
            if (args.length >= 9) {
                makeYAxisEqual = Boolean.parseBoolean(args[8]);
            }

            try {
                p.plot(efilename, gtf, query, start, stop, output, limitRegionTo100kbAroundTopSNP, maxP, makeYAxisEqual);
            } catch (IOException | DocumentException e) {
                e.printStackTrace();
            }
        }

    }


    public void plot(String eqtlfiletemplate, String gtf, String querygene, int startiter, int stopiter, String output, boolean limitRegionTo100kbAroundTopSNP, Double presetMaxP, boolean makeYAxisEqual) throws IOException, DocumentException {


        System.out.println(eqtlfiletemplate);
        System.out.println(gtf);
        System.out.println(querygene);
        System.out.println(startiter);
        System.out.println(stopiter);
        System.out.println(output);
        System.out.println(limitRegionTo100kbAroundTopSNP);
        System.out.println(presetMaxP);
        System.out.println(makeYAxisEqual);

        GTFAnnotation f = new GTFAnnotation(gtf);
        Grid grid = new Grid(300, 150, 10, 1, 100, 100);
        GenePanel genePanel = new GenePanel(1, 1);
        grid.addPanel(genePanel);

        AssociationPanel assocPanel = new AssociationPanel(1, 1);
        grid.addPanel(assocPanel);


        Chromosome chr = null;
        int regionstart = Integer.MAX_VALUE;
        int regionstop = 0;
        double maxp = 0;
        int maxPPos = 0;

        ArrayList<ArrayList<Pair<Integer, Double>>> allPValues = new ArrayList<>();
        ArrayList<String> datasetNames = new ArrayList<>();
        for (int i = startiter; i < stopiter + 1; i++) {
            String file = eqtlfiletemplate.replaceAll("ITER", "" + i);
            datasetNames.add("Iteration" + i);
            TextFile tf = new TextFile(file, TextFile.R);
            System.out.println("Reading file " + file + " looking for " + querygene);
            tf.readLine();
            int ctr = 0;
            String[] elems = tf.readLineElems(TextFile.tab);
            ArrayList<Pair<Integer, Double>> pvals = new ArrayList<Pair<Integer, Double>>();


            while (elems != null) {

                String gene = elems[4];
                if (gene.equals(querygene)) {
                    Double p = -Math.log10(Double.parseDouble(elems[0]));
                    Integer pos = Integer.parseInt(elems[3]);
                    if (p > maxp) {
                        maxp = p;
                        maxPPos = pos;
                        if (limitRegionTo100kbAroundTopSNP) {

                            regionstop = maxPPos + 100000;
                            regionstart = maxPPos - 100000;
                        }
                    }
                    if (!limitRegionTo100kbAroundTopSNP) {
                        if (pos > regionstop) {
                            regionstop = pos + 100;
                        }
                        if (pos < regionstart) {
                            regionstart = pos - 100;
                        }
                    }
                    if (chr == null) {
                        chr = Chromosome.parseChr(elems[2]);
                    }
                    pvals.add(new Pair<>(pos, p));
                }
                ctr++;
                if (ctr % 1000000 == 0) {
                    System.out.print(ctr + " lines parsed. " + pvals.size() + " pvals sofar. MAX: " + maxPPos + ", " + maxp + "\tRegion: " + chr + " " + regionstart + " - " + regionstop + "\r");
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();

            System.out.println();
            allPValues.add(pvals);
            System.out.println("Iteration " + i + " has " + pvals.size() + " pvalues ");
        }

        if (limitRegionTo100kbAroundTopSNP) {
            System.out.println("Pruning results.");

            ArrayList<ArrayList<Pair<Integer, Double>>> tmp = new ArrayList<>();
            for (int i = 0; i < allPValues.size(); i++) {
                ArrayList<Pair<Integer, Double>> curP = allPValues.get(i);
                ArrayList<Pair<Integer, Double>> nextP = new ArrayList<>();
                for (Pair<Integer, Double> pval : curP) {
                    if (pval.getLeft() > regionstart && pval.getLeft() < regionstop) {
                        nextP.add(pval);
                    }
                }
                tmp.add(nextP);
            }
            allPValues = tmp;
        }


        Feature region = new Feature(chr, regionstart, regionstop);


        ArrayList<Gene> genes = new ArrayList<>();
        for (Gene gobj : f.getGenes()) {
            if (gobj.overlaps(region)) {
                if (gobj.getName().equals(querygene)) {
                    genes.add(gobj);
                }
            }
        }
        genePanel.setData(region, genes);


        int idx = 0;

        Color[] colors = new Color[5];
        for (int i = 0; i < colors.length; i++) {
            switch (i) {
                case 0:
                    colors[i] = new Color(70, 67, 58, 150);
                    break;
                case 1:
                    colors[i] = new Color(174, 164, 140, 150);
                    break;
                case 2:
                    colors[i] = new Color(231, 79, 19, 150);
                    break;
                case 3:
                    colors[i] = new Color(98, 182, 177, 150);
                    break;
                case 4:
                    colors[i] = new Color(116, 156, 80, 150);
                    break;

            }
        }
        assocPanel.setData(region, null, allPValues, datasetNames.toArray(new String[0]));
        if (makeYAxisEqual) {
            if (presetMaxP != null) {
                assocPanel.setMaxPVal(presetMaxP);
            } else {
                assocPanel.setMaxPVal(maxp);
            }
        }
        assocPanel.setColors(colors);

        for (int i = startiter; i < stopiter + 1; i++) {
            AssociationPanel iterAssocPanel = new AssociationPanel(1, 1);
            iterAssocPanel.setDataSingleDs(region, null, allPValues.get(idx), "Iteration " + i);
            iterAssocPanel.setColors(new Color[]{colors[idx]});
            grid.addPanel(iterAssocPanel);
            if (makeYAxisEqual) {
                if (presetMaxP != null) {
                    assocPanel.setMaxPVal(presetMaxP);
                } else {
                    iterAssocPanel.setMaxPVal(maxp);
                }
            }


            idx++;
        }

        grid.draw(output);
        System.out.println("Output is here: " + output);
    }

}
