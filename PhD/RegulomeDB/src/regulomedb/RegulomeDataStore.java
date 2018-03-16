/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package regulomedb;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class RegulomeDataStore {

    private static final String MAGIC = "\n";
    HashMap<String, Integer> snpIndex = null;
    HashMap<String, Integer> classIndex = null;

    public RegulomeDataStore(String db) throws IOException {
        snpIndex = new HashMap<String, Integer>();
        classIndex = new HashMap<String, Integer>();
        init(db);
    }

    private void init(String dir) throws IOException {
    }

    public void convertToMatrix(String regdbdir) throws IOException {
        // determine unique classes

        HashSet<String> uniqueClasses = new HashSet<String>();
        String[] regdbfiles = Gpio.getListOfFiles(regdbdir, "gz");
        HashSet<Pair<String, String>> eqtls = new HashSet<Pair<String, String>>();

        HashSet<String> uniqueSNPs = new HashSet<String>();
        for (String f : regdbfiles) {
            String filename = f;
            System.out.println("Parsing file: " + filename);
            TextFile tf = new TextFile(filename, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {

                String chr = elems[0];
                String snp = elems[2];
                String pos = elems[1];
                String annot = elems[3];

                uniqueSNPs.add(snp);

//                System.out.println(Strings.concat(elems, Strings.tab));
                String[] annotelems = annot.split(", ");
                for (String s : annotelems) {
//                    System.out.println(s);
                    if (s.contains("eQTL")) {
                        uniqueClasses.add("eQTL");
                        // add eQTL
//                        System.out.println(s);
                        String gene = s.split("|")[2];
                        Pair<String, String> eqtl = new Pair<String, String>(snp, gene);
                        eqtls.add(eqtl);

                    } else {
                        uniqueClasses.add(s);
                    }

                }

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }

        System.out.println("");
        System.out.println("Total unique classes: " + uniqueClasses.size());
        System.out.println("Total unique snps: " + uniqueSNPs.size());

        // magic character is null      
        BinaryFile bf = new BinaryFile(regdbdir + "RegulomeDB.dat", BinaryFile.W);

        // first write all eQTLs
        for (Pair<String, String> eqtl : eqtls) {
            bf.writeString(eqtl.getLeft());
            bf.writeString(eqtl.getRight());
        }


        bf.writeString(MAGIC);
        // now write the boolean matrix.. first write the categories
        HashMap<String, Integer> classIndex = new HashMap<String, Integer>();
        HashMap<String, Integer> snpIndex = new HashMap<String, Integer>();
        int ctr = 0;
        for (String s : uniqueClasses) {
            bf.writeString(s);
            classIndex.put(s, ctr);
            ctr++;
        }
        bf.writeString(MAGIC);
        // write the SNP names
        ctr = 0;
        for (String s : uniqueSNPs) {
            bf.writeString(s);
//            snpIndex.put(s, ctr);
            ctr++;
        }
        bf.writeString(MAGIC);




        // load the matrix into memory
        // now write the matrix.

        HashMap<String, boolean[]> snpToAnnotation = new HashMap<String, boolean[]>();
        for (String f : regdbfiles) {
            String filename = f;
            System.out.println("Parsing file: " + filename);
            TextFile tf = new TextFile(filename, TextFile.R);
            String[] elems = tf.readLineElems(TextFile.tab);

            while (elems != null) {

                String chr = elems[0];
                String snp = elems[2];
                String pos = elems[1];
                String annot = elems[3];

                boolean[] byteAnnotation = snpToAnnotation.get(snp);
                if (byteAnnotation == null) {
                    byteAnnotation = new boolean[uniqueClasses.size()];
                }


//                System.out.println(Strings.concat(elems, Strings.tab));
                String[] annotelems = annot.split(", ");
                for (String s : annotelems) {
//                    System.out.println(s);
                    if (s.contains("eQTL")) {
                        uniqueClasses.add("eQTL");
                        // add eQTL
//                        System.out.println(s);
                        String gene = s.split("|")[2];
                        Pair<String, String> eqtl = new Pair<String, String>(snp, gene);
                        eqtls.add(eqtl);

                    } else {
                        uniqueClasses.add(s);
                        int index = classIndex.get(s);
                        byteAnnotation[index] = true;
                    }

                }
                snpToAnnotation.put(snp, byteAnnotation);

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
        }

        for (String s : uniqueSNPs) {
            boolean[] arr = snpToAnnotation.get(s);
            for (boolean b : arr) {
                bf.writeBool(b);
            }


        }

        bf.close();
    }
}
