/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;
import umcg.genetica.ncbi.dbsnp.SNPAnnotation;

/**
 *
 * @author harmjan
 */
public class GetFunctionalAnnotationForSNPs {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {


            String queryfilename = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-10-06-TestedGWASSNPs.txt";
            String contiglocidfile = "/Volumes/iSnackHD/Data/SNPReferenceData/dbSNP/b130/b130_SNPContigLocusId_36_3.bcp";
            String snpfuncclassfile = "/Volumes/iSnackHD/Data/SNPReferenceData/dbSNP/b130/SnpFunctionCode.bcp";
            String reference = "/Volumes/iSnackHD/Data/GWAS/HapMapCEU/TriTyper/";
            // String reference = "/Volumes/iSnackHD/Data/SNPReferenceData/1000G-20110521-TriTyper/Merged/"
            String outdir = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/GWASSNPDBSNPAnnotations/";

            
            
            HashSet<String> querySNPs = new HashSet<String>();
            TextFile queryfile = new TextFile(queryfilename, TextFile.R);
            querySNPs.addAll(queryfile.readAsArrayList());
            queryfile.close();


            TriTyperGenotypeData ds = new TriTyperGenotypeData(reference);
            SNPLoader loader = ds.createSNPLoader();
            String[] snps = ds.getSNPs();
            DetermineLD ldcalc = new DetermineLD();
            HashSet<String> finalSNPList = new HashSet<String>();
            HashMap<String, HashSet<String>> proxiesForSNPs = new HashMap<String, HashSet<String>>();
            
            for (String query : querySNPs) {
                Integer snpId = ds.getSnpToSNPId().get(query);

                if (snpId == null) {
                    System.err.println("ERROR: snp " + query + " not in reference");
                } else {
                    byte snpchr = ds.getChr(snpId);
                    int pos = ds.getChrPos(snpId);
                    SNP snp1 = ds.getSNPObject(snpId);
                    HashSet<String> proxies = new HashSet<String>();
                    if (snpchr != -1) {
                        loader.loadGenotypes(snp1);
                        for (int i = 0; i < snps.length; i++) {
                            if (i != snpId && ds.getChr(i) == snpchr) {

                                int pos2 = ds.getChrPos(i);
                                if (Math.abs(pos - pos2) < 5000000) {
                                    SNP snp2 = ds.getSNPObject(i);
                                    loader.loadGenotypes(snp2);

                                    double r2 = ldcalc.getRSquared(snp1, snp2, ds, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                                    if (r2 >= 0.8) {
                                        proxies.add(snps[i]);
                                    }
                                    snp2.clearGenotypes();
                                }
                            }
                        }

                        proxiesForSNPs.put(query, proxies);


                        snp1.clearGenotypes();
                    }
                    System.out.println(query + "\t" + snpchr + "\t" + proxies.size());
                    finalSNPList.add(query);
                    finalSNPList.addAll(proxies);
                }
            }
            loader.close();

            HashSet<String> originalQuerySNPs = querySNPs;
            querySNPs = finalSNPList;
            System.out.println("Final SNP list: " + querySNPs.size());

            SNPAnnotation s = new SNPAnnotation("/Volumes/iSnackHD/Data/SNPReferenceData/dbSNP/b130/b130_ContigInfo_36_3.bcp", "/Volumes/iSnackHD/Data/SNPReferenceData/dbSNP/b130/b130_SNPContigLoc_36_3.bcp", "reference");
            HashMap<Integer, Integer> rsToChrPos = s.rsToChrPos;
            HashSet<String> contigs = s.getContigs();

            HashMap<String, String> snpfuncclassses = new HashMap<String, String>();
            TextFile funcclassfile = new TextFile(snpfuncclassfile, TextFile.R);

            String[] elems = funcclassfile.readLineElems(TextFile.tab);
            while (elems != null) {
                String classid = elems[0];
                String funcdesc = elems[1];
                snpfuncclassses.put(classid, funcdesc);
                elems = funcclassfile.readLineElems(TextFile.tab);
            }
            funcclassfile.close();


            TextFile tf = new TextFile(contiglocidfile, TextFile.R);
            elems = tf.readLineElems(TextFile.tab);
            HashSet<String> snpsWAnnot = new HashSet<String>();
            TextFile outfile = new TextFile(outdir+"Annotation.txt", TextFile.W);
            while (elems != null) {
                String snp = "rs" + elems[0];
                String contig = elems[10];

                if (querySNPs.contains(snp)) {
                    String funcid = elems[11];
                    String func = snpfuncclassses.get(funcid);
                    String geneId = elems[5];
                    String gene = elems[6];
                    if (!func.equals("cds-reference")) {
                        System.out.println(snp + "\t" + contig + "\t" + funcid + "\t" + func + "\t" + geneId + "\t" + gene);
                        outfile.writeln(snp + "\t" + contig + "\t" + funcid + "\t" + func + "\t" + geneId + "\t" + gene);
                    }
                    snpsWAnnot.add(snp);
                }
                elems = tf.readLineElems(TextFile.tab);
            }
            outfile.close();
            tf.close();

            System.out.println(snpsWAnnot.size());
            System.out.println("");
            System.out.println("");

            TextFile proxyout = new TextFile(outdir+"proxies.txt", TextFile.W);
            for (String snp : originalQuerySNPs) {
                HashSet<String> proxies = proxiesForSNPs.get(snp);
                if (proxies != null) {
                    for (String snp2 : proxies) {
                        System.out.println(snp + "\t" + snp2);
                        proxyout.writeln(snp+"\t"+snp2);
                    }
                }
            }
            proxyout.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }
}
