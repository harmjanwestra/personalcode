/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ListCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String f1 = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNPs-Rerun/eQTLs.txt";
            String f2 = "/Volumes/iSnackHD/Data/Projects/Isis/eQTLResults/2012-10-30-eQTLs-Cis/SNPs-Rerun/SNPQCLog.txt.gz";
            int col1 = 1;
            int col2 = 0;
            boolean hasheader1 = true;
            boolean hasheader2 = true;
            compare(f1, col1, hasheader1, f2, col2, hasheader2);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void compare(String f1, int col1, boolean file1HasHeader, String f2, int col2, boolean file2HasHeader) throws IOException {
        HashSet<String> elems1 = new HashSet<String>();
        HashSet<String> elems2 = new HashSet<String>();
        TextFile in = new TextFile(f1, TextFile.R);
        String[] elems = in.readLineElems(Strings.whitespace);
        if (file1HasHeader) {
            elems = in.readLineElems(Strings.whitespace);
        }
        while (elems != null) {
            elems1.add(elems[col1]);
            elems = in.readLineElems(Strings.whitespace);
        }
        in.close();


        TextFile in2 = new TextFile(f2, TextFile.R);
        elems = in2.readLineElems(Strings.whitespace);
        if (file1HasHeader) {
            elems = in2.readLineElems(Strings.whitespace);
        }
        while (elems != null) {
            elems2.add(elems[col2]);
            elems = in2.readLineElems(Strings.whitespace);
        }
        in2.close();

        System.out.println("list1: " + elems1.size());
        System.out.println("list2: " + elems2.size());

        int ctr = 0;
        for (String s : elems1) {
            if (elems2.contains(s)) {
                ctr++;
            }
        }
        System.out.println("Overlap: " + ctr);
    }
}
