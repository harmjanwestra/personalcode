/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import java.util.*;
import java.util.concurrent.Callable;
import umcg.genetica.io.ExpressionDataset;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.WilcoxonMannWhitney;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ZScoreTableSNPPathwayTTestTask implements Callable<String> {

    private final int snp;
    private final ArrayList<String> pathways;
    private final HashMap<String, HashSet<String>> genesInPathways;
//    private final DoubleMatrixDataset<String, String> ds;
    private final DoubleMatrixDataset<String,String> ds;
    private final GWASCatalog c;

//    public ZScoreTableSNPPathwayTTestTask(int row, DoubleMatrixDataset<String, String> ds, ArrayList<String> pathways, HashMap<String, HashSet<String>> genesInPathways, GWASCatalog c) {
    public ZScoreTableSNPPathwayTTestTask(int row, DoubleMatrixDataset<String,String> ds, ArrayList<String> pathways, HashMap<String, HashSet<String>> genesInPathways, GWASCatalog c) {
        this.c = c;
        this.snp = row;
        this.ds = ds;
        this.pathways = pathways;
        this.genesInPathways = genesInPathways;
    }

    @Override
    public String call() throws Exception {

//	System.out.println("Running task :D");
        List<String> header = ds.colObjects;
        double[] zscores = ds.rawData[snp];// ds.rawData[snp];
//	String snpname = ds.rowObjects.get(snp);

        String snpname = ds.rowObjects.get(snp);

//	System.out.println(snpname);
        GWASSNP snpObj = c.getSnpToObj().get(snpname);
        GWASTrait[] traits = snpObj.getAssociatedTraitsArray();
        String[] traitnames = new String[traits.length];
        int ctr = 0;
        for (GWASTrait t : traits) {
            traitnames[ctr] = t.getName();
            ctr++;
        }

        String traitsassoc = Strings.concat(traitnames, Strings.comma);
        String output = null;
        StringBuilder sb = new StringBuilder();
        int pathwaysSignificant = 0;
        for (String pathway : pathways) {

            HashSet<String> genesInThisPathway = genesInPathways.get(pathway);
            HashSet<Integer> columnsInPathway = new HashSet<Integer>();
            for (int i = 0; i < header.size(); i++) {
                if (genesInThisPathway.contains(header.get(i))) {
                    columnsInPathway.add(i);
                }
            }

//	    System.out.println("Testing pathway: " + pathway + "\t" + columnsInPathway.size());

            if (columnsInPathway.size() >= 25) {

                double[] selectedGenes = new double[columnsInPathway.size()];
//		double[] otherGenes = new double[ds.colObjects.size() - columnsInPathway.size()];
                double[] otherGenes = new double[ds.colObjects.size() - columnsInPathway.size()];

                double[] selectedGenesAbs = new double[columnsInPathway.size()];
//		double[] otherGenesAbs = new double[ds.colObjects.size() - columnsInPathway.size()];
                double[] otherGenesAbs = new double[ds.colObjects.size() - columnsInPathway.size()];
                int genesIn = 0;
                int genesOut = 0;
                //System.out.println("ZScore For RNF: "+elems[colForRNF]);
                for (int i = 0; i < ds.colObjects.size(); i++) {

                    if (columnsInPathway.contains(i)) {
                        selectedGenes[genesIn] = zscores[i];
                        selectedGenesAbs[genesIn] = Math.abs(zscores[i]);
//                        System.out.println(header.get(i) + "\t" + selectedGenes[genesIn] + "\t" + selectedGenesAbs[genesIn] + "\t1");
                        genesIn++;
                    } else {
                        otherGenes[genesOut] = zscores[i];
                        otherGenesAbs[genesOut] = Math.abs(zscores[i]);
//                        System.out.println(header.get(i) + "\t" + otherGenes[genesOut] + "\t" + otherGenesAbs[genesOut] + "\t0");
                        genesOut++;
                    }
                }




                WilcoxonMannWhitney w = new WilcoxonMannWhitney();
                double wpval = w.returnWilcoxonMannWhitneyPValue(selectedGenes, otherGenes);
                double auc = w.getAUC();

                w = new WilcoxonMannWhitney();
                double wpvalAbs = w.returnWilcoxonMannWhitneyPValue(selectedGenesAbs, otherGenesAbs);
                double aucAbs = w.getAUC();
//		if (wpval < 1e-5 || wpvalAbs < 1e-5) {


                sb.append(snpname).append("\t");
                sb.append(traitsassoc).append("\t");
                sb.append(pathway).append("\t");
                sb.append(genesInThisPathway.size()).append("\t");
                sb.append(columnsInPathway.size()).append("\t");
                sb.append(otherGenes.length).append("\t");
                sb.append(wpval).append("\t");
                sb.append(auc).append("\t");
                sb.append(wpvalAbs).append("\t");
                sb.append(aucAbs);
                sb.append("\n");
//		    pathwaysSignificant++;
//		System.out.println(sb.toString());
//		}
            }
        }
//	if (pathwaysSignificant > 0) {
        return sb.toString();
//	}
//	return output;
    }
}
