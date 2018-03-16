/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class DetermineTwoSNPFileOverlap {

    public static void main(String[] args) {
        try {
            DetermineTwoSNPFileOverlap.run(
                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscle/LiverCyto/Individuals.txt",
                    "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscle/LiverOmni/Individuals.txt");
//            DetermineTwoSNPFileOverlap.snpOverlapinEQTLFile(
//                    "/Volumes/iSnackHD/Patrick/GonlImputed-eQTLsFDR0.05.txt",
//                    "/Volumes/iSnackHD/Patrick/NonImputed-eQTLsFDR0.05.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void run(String file1, String file2) throws IOException {
        HashSet<String> elemsInFile1 = new HashSet<String>();
        TextFile tf1 = new TextFile(file1, TextFile.R);
        elemsInFile1.addAll(tf1.readAsArrayList());
        tf1.close();

        TextFile tf2 = new TextFile(file2, TextFile.R);
        String ln = tf2.readLine();
        HashSet<String> elemsInFile2 = new HashSet<String>();
        int shared = 0;
        while (ln != null) {
            if (elemsInFile1.contains(ln.trim())) {
                shared++;
                System.out.println(ln);
            }
            elemsInFile2.add(ln);
            ln = tf2.readLine();
        }
        tf2.close();

        System.out.println("Elems in file 1: " + elemsInFile1.size());
        System.out.println("Elems in file 2: " + elemsInFile2.size());
        System.out.println("Shared: " + shared);
    }

    public static void snpOverlapinEQTLFile(String f1, String f2) throws IOException {
        TextFile tf1 = new TextFile(f1, TextFile.R);

        String[] elems = tf1.readLineElems(TextFile.tab);
        
        HashSet<String> snps = new HashSet<String>();
        elems = tf1.readLineElems(TextFile.tab);
        while(elems!=null){
            snps.add(elems[1]);
            elems = tf1.readLineElems(TextFile.tab);
        }
        tf1.close();

        TextFile in = new TextFile(f2, TextFile.R);
        elems = in.readLineElems(TextFile.tab);
        HashSet<String> snpsin2 = new HashSet<String>();
        while (elems != null) {
            snpsin2.add(elems[1]);
            elems = in.readLineElems(TextFile.tab);
        }
        
        in.close();
        
        int ctr = 0;
        for(String snp: snpsin2){
            if(snps.contains(snp)){
                ctr++;
                
            }
        }
        System.out.println(ctr+" shared");
    }
}
