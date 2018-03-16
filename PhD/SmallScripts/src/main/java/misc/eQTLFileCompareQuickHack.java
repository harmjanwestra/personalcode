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
public class eQTLFileCompareQuickHack {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        // TODO code application logic here
        String file1 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD4-eQTLs.txt.gz";
        String file2 = "/Volumes/iSnackHD/AeroFS/cellTypeeQTL/Replication/2013-09-04-CD4AndCD8Cells/CD8-eQTLs.txt.gz";
        
        try {
            TextFile tf = new TextFile(file1, TextFile.R);
            tf.readLine();
            String[] elems = tf.readLineElems(TextFile.tab);
            HashSet<String> probesInFile1 = new HashSet<String>();
            HashSet<String> snpsInFile1 = new HashSet<String>();
            
            
            HashSet<String> probesInFile2 = new HashSet<String>();
            HashSet<String> snpsInFile2 = new HashSet<String>();
            
            while (elems != null) {

                String probe = elems[4];
                String snp = elems[1];
                
                snpsInFile1.add(snp);
                probesInFile1.add(probe);
                

                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();
            
            
            tf = new TextFile(file2, TextFile.R);
            tf.readLine();
            elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {

                
                String probe = elems[4];
                String snp = elems[1];
                
                snpsInFile2.add(snp);
                probesInFile2.add(probe);
                

                elems = tf.readLineElems(TextFile.tab);
            }
            tf.close();
            
            for(String snp: snpsInFile1){
                if(!snpsInFile2.contains(snp)){
                    System.out.println("SNP missing in CD8 file: "+snp);
                }
            }
            
            for(String snp: snpsInFile2){
                if(!snpsInFile1.contains(snp)){
                    System.out.println("SNP missing in CD4 file: "+snp);
                }
            }
            
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
