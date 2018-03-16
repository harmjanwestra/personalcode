/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class eQTLSNPFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            String qfile = "/Volumes/iSnackHD/Data/Projects/MarkDaly/2013-01-16-IBD/IBC_loci.txt";
            String efile = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/cis/2012-07-18/Sorted/2012-10-23-PostQC/eQTLs.txt-FilteredForProbeLevelFDR.txt.gz";
            String ofile = "/Volumes/iSnackHD/Data/Projects/MarkDaly/2013-01-16-IBD/IBD_loci-EQTLsFromLargeMetaAnalysis-FDR0.05.txt";
            eQTLSNPFilter s = new eQTLSNPFilter();
            s.run(qfile, efile, ofile);
        } catch (IOException e){
            e.printStackTrace();
        }
    }
    
    public void run(String query, String efile, String ofile) throws IOException {
        TextFile tf = new TextFile(query, TextFile.R);
        HashSet<String> querySNPs = new HashSet<String>();
        querySNPs.addAll(tf.readAsArrayList());
        tf.close();
        
        System.out.println(querySNPs.size()+" SNPs loaded");
        
        TextFile eqtls = new TextFile(efile, TextFile.R);
        TextFile out = new TextFile(ofile, TextFile.W);
        out.writeln(eqtls.readLine());
        String[] elems = eqtls.readLineElems(TextFile.tab);
        HashSet<String> detectedSNPs = new HashSet<String>();
        while(elems!=null){
            String snp = elems[1];
            if(querySNPs.contains(snp)){
                out.writeln(Strings.concat(elems, Strings.tab));
                detectedSNPs.add(snp);
            }
            elems = eqtls.readLineElems(TextFile.tab);
        }
        eqtls.close();
        out.close();
        
        int n=0;
        System.out.println("SNPs without eQTL");
        for(String s: querySNPs){
            if(!detectedSNPs.contains(s)){
                System.out.println(s);
                n++;
            }
        }
        System.out.println(n);
        
    }
}
