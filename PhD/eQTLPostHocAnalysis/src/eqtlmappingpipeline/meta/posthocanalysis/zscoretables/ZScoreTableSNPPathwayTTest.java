/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import umcg.genetica.io.ExpressionDataset;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ZScoreTableSNPPathwayTTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {

            //TextFile tfpathway = new TextFile("/Volumes/iSnackHD/JuhaPathways/go_bp/GO_BP+.gmt", TextFile.R);
            TextFile tfpathway = new TextFile("/Volumes/iSnackHD/PathwayData/2012-06-01-PathwayGMTFiles/reactome/Reactome+.gmt", TextFile.R);
            //TextFile tfpathway = new TextFile("/Volumes/iSnackHD/JuhaPathways/kegg/KEGG+.gmt", TextFile.R);

//	    TextFile tfpathway = new TextFile("/Volumes/iSnackHD/PathwayData/2012-05-08-PathwayGMTFiles/OriginalAssignment-GMTs/Reactome/ReactomePathwaysENSG.gmt", TextFile.R);
//	    HashSet<String> genesInPathway = new HashSet<String>();
            String[] pwelems = tfpathway.readLineElems(TextFile.tab);

            ArrayList<String> pathways = new ArrayList<String>();
            HashMap<String, HashSet<String>> genesInPathways = new HashMap<String, HashSet<String>>();

            GWASCatalog gwasCatalog = new GWASCatalog();
            gwasCatalog.read("/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-24-gwascatalog.txt");

            

            System.out.println(pathways.size() + " pathways");
            System.out.println(genesInPathways.size() + " genes detected in pathway");

            String indir = "/Volumes/iSnackHD/MetaAnalysisFinal/cistrans/2012-05-30-Cistrans/";
            String outdir = indir + "PathwayAnnotation/AllPathwaysVsPCSNP/";
            Gpio.createDir(outdir);
            for (int i = 2; i < 10; i++) {
                DoubleMatrixDataset<String,String> data = null;
                TextFile output = null;
                if (i == 0) {
                    if (!Gpio.exists(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.binary.dat")) {
                        data = new DoubleMatrixDataset<String,String>(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.txt");

                        // data = new DoubleMatrixDataset<String, String>(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.txt");
                        data.save(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.binary");
                    } else {
                        data = new DoubleMatrixDataset<String,String>(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.binary");
                        //data = new DoubleMatrixDataset<String, String>(indir + "metazscoretable.txt.gz-filtered.txt-ens.txt-collapsed.binary");
                    }

                    output = new TextFile(outdir + "Reactome-Plus-TTest.txt", TextFile.W);

                } else {
                    //

                    // /Volumes/iSnackHD/MetaAnalysisFinal/cistrans/2012-05-30-Cistrans/PathwayAnnotation/ProstateCancerPW/
//		    data = new ExpressionDataset(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.txt");
//		    data = new DoubleMatrixDataset<String, String>(indir+"metazscoretable-Permutation"+i+".txt.gz-filtered.txt-ens.txt-collapsed.txt");
                    if (!Gpio.exists(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.binary.dat")) {
//			data = new DoubleMatrixDataset<String, String>(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.txt");
                        data = new DoubleMatrixDataset<String,String>(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.txt");
                        data.save(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.binary");
                    } else {
                        data = new DoubleMatrixDataset<String,String>(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.binary");
//                        data = new DoubleMatrixDataset<String, String>(indir + "metazscoretable-Permutation" + i + ".txt.gz-filtered.txt-ens.txt-collapsed.binary");
                    }
                    output = new TextFile(outdir + "Reactome-Plus-TTest-Permutation" + i + ".txt", TextFile.W);
                }

                System.out.println("Writing output to: " + output.getFileName());

                int nrprocs = Runtime.getRuntime().availableProcessors();
                ExecutorService threadPool = Executors.newFixedThreadPool(nrprocs);
                CompletionService<String> pool = new ExecutorCompletionService<String>(threadPool);

                // rows contain the SNPs
                int tasksSubmitted = 0;
                for (int row = 0; row < data.rowObjects.size(); row++) {
                    String snpname = data.rowObjects.get(row);
//		    System.out.println(snpname);
//                    if (snpname.equals("rs10187424")) {
                        ZScoreTableSNPPathwayTTestTask task = new ZScoreTableSNPPathwayTTestTask(row, data, pathways, genesInPathways, gwasCatalog);
                        pool.submit(task);
                        tasksSubmitted++;
//                    }
                }

                System.out.println(tasksSubmitted + " tasks submitted..");
                output.writeln("snp\tassocTraits\tpathway\tgenesInPathway\tpathwayGenesTested\totherGenesTested\twilcoxonP\twilcoxonAUC\twilcoxonPAbs\twilcoxonAUCAbs");

                int returnedResults = 0;

                while (returnedResults < tasksSubmitted) {
                    try {
                        String result = pool.take().get();
                        if (result != null) {
//                            System.out.print(result);
                            output.write(result);
                            returnedResults++;
                        }


                        double perc = Math.ceil((double) returnedResults / data.rowObjects.size());
                        if (perc % 10 == 0) {
                            System.out.println(returnedResults + "/" + data.rowObjects.size() + "\t" + perc + "%");

                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }

                threadPool.shutdown();

                output.close();
            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
