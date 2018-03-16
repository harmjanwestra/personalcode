/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.celltypespecificmetaanalysis.forpaper;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.probeannotation.ProbeTranslation;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harm-jan
 */
public class DetermineNumberOfSignificantEffectsPerCovariate {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String fdrfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/MetaAnalysisZScoreMatrix.binary-FDR.binary";
        String metamatrix = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/MetaAnalysis/MetaAnalysisZScoreMatrix.binary";
        String probetranslationfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-18-ProbeAnnotationFile.txt";
        double fdrthreshold = 0.05;
        String output = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/MetaAnalysisResults/2013-08-20-MetaAnalysis7DatasetsEffectsNotFlipped/2013-09-27-FDREstimates/FDRSignificantEffectsPerCovariate-FDR" + fdrthreshold + ".txt";
        String gwascatalog = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-09-23-gwascatalog.txt";
        double threshold = 5E-8;
        String efile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2013-07-11-HT12v3SNPProbeCombos.txt";
        String eqtlfile = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/DataFiles/2012-10-23-eQTLsFDR-AllProbeLevelQTLs.txt.gz";
        try {
            DetermineNumberOfSignificantEffectsPerCovariate cov = new DetermineNumberOfSignificantEffectsPerCovariate();
            cov.run(fdrfile, metamatrix, fdrthreshold, probetranslationfile, gwascatalog, threshold, output, efile, eqtlfile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String fdrfile, String metamtrix, double fdrthreshold, String probetranslationfile, String gwascatalog, double gwasthreshold, String output, String efile, String eqtlfile) throws IOException {
        DoubleMatrixDataset<String, String> meta = new DoubleMatrixDataset<String, String>(metamtrix);
        DoubleMatrixDataset<String, String> fdrmatrix = new DoubleMatrixDataset<String, String>(fdrfile);

        GWASCatalog catalog = new GWASCatalog(gwascatalog);
        HashSet<String> significantSNPs = new HashSet<String>();
        HashSet<String> allGWASSNPs = new HashSet<String>();
        for (GWASTrait t : catalog.getTraits()) {
            GWASSNP[] snps = t.getSNPs(gwasthreshold);
            GWASSNP[] snps2 = t.getSNPs();
            
            for (GWASSNP snp : snps) {
                significantSNPs.add(snp.getName());
            }
            for (GWASSNP snp : snps2) {
                allGWASSNPs.add(snp.getName());
            }
        }
        HashSet<String> uniqueGWASSNPs = new HashSet<String>();
        ProbeTranslation pb = new ProbeTranslation();
        HashMap<String, String> ht12toGene = pb.getProbeTranslation(probetranslationfile, "HT12v3.txt", "Gene");
        TextFile outfile = new TextFile(output, TextFile.W);
        outfile.writeln("probe\tgene\tzneg\tzother\tzpos\tzneg gwas\tzother gwas\tzpos gwas");
        for (int covariate = 0; covariate < meta.nrRows; covariate++) {

            String probe = meta.rowObjects.get(covariate);
            String gene = ht12toGene.get(probe);
            int zpos = 0;
            int zneg = 0;
            int zoth = 0;

            int zposg = 0;
            int znegg = 0;
            int zothg = 0;
            for (int eqtl = 0; eqtl < meta.nrCols; eqtl++) {
                String snp = meta.colObjects.get(eqtl).split("-")[0];
                if (significantSNPs.contains(snp)) {
                    uniqueGWASSNPs.add(snp);
                }
                double fdr = fdrmatrix.rawData[covariate][eqtl];
                double z = meta.rawData[covariate][eqtl];
                if (fdr < fdrthreshold) {
                    if (z > 0) {
                        zpos++;
                        if (significantSNPs.contains(snp)) {
                            zposg++;
                        }
                    } else {
                        zneg++;
                        if (significantSNPs.contains(snp)) {
                            znegg++;
                        }
                    }
                } else {
                    zoth++;
                    if (significantSNPs.contains(snp)) {
                        zothg++;
                    }
                }
            }
            // System.out.println(uniqueGWASSNPs.size() + " GWAS SNPs!");
            outfile.writeln(probe + "\t" + gene + "\t" + zneg + "\t" + zoth + "\t" + zpos
                    + "\t" + znegg + "\t" + zothg + "\t" + zposg);

        }
        outfile.close();

        HashSet<String> visitedProbes = new HashSet<String>();
        HashSet<String> topSNPs = new HashSet<String>();
        TextFile tf2 = new TextFile(eqtlfile, TextFile.R);
        String[] elems2 = tf2.readLineElems(TextFile.tab);
        while(elems2!=null){
            String snp = elems2[1];
            String probe = elems2[4];
            if(!visitedProbes.contains(probe)){
                visitedProbes.add(probe);
                topSNPs.add(snp);
            }
            elems2 = tf2.readLineElems(TextFile.tab);
        }
        tf2.close();
        
        HashSet<String> otherSNPs = new HashSet<String>();
        TextFile tf1 = new TextFile(efile, TextFile.R);
        String[] elems = tf1.readLineElems(TextFile.tab);
        HashSet<String> uniqueGWASSNPs2 = new HashSet<String>();
        HashSet<String> uniqueSNPs = new HashSet<String>();
        HashSet<String> uniqueLines = new HashSet<String>();
        HashSet<String> uniqueTopSNPs= new HashSet<String>();
        
        while (elems != null) {

            uniqueLines.add(elems[0]+"\t"+elems[1]);
            String snp = elems[0];
            uniqueSNPs.add(snp);
            if(allGWASSNPs.contains(snp)){
                uniqueGWASSNPs2.add(snp);
            } else if(topSNPs.contains(snp)){
                uniqueTopSNPs.add(snp);
            } else {
                otherSNPs.add(snp);
            }
            
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();
        
        System.out.println("Unique: "+uniqueSNPs.size());
        System.out.println("GWAS: "+uniqueGWASSNPs2.size());
        System.out.println("Top SNPs: "+uniqueTopSNPs.size());
        System.out.println("Other: "+otherSNPs.size());
        System.out.println("uniqueLines: "+uniqueLines.size());

    }
}
