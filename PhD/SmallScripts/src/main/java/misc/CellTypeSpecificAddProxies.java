/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package misc;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class CellTypeSpecificAddProxies {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try{
            String file = "";
            String[] proxyfiles = new String[]{"/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/100GenomesCEU-Proxies-0.99.txt","/Volumes/iSnackHD/AeroFS/2013-12-09-ProxiesForASE/HapMap2CEU-Proxies-0.99.txt"};
            TextFile tf = new TextFile(file, TextFile.R);
            Set<String> query = tf.readAsSet(0, TextFile.tab);
            HashMap<String, String> proxiesPerSNP = new HashMap<String, String>();
            for(String proxy: proxyfiles){
                TextFile tf2 = new TextFile(proxy, TextFile.R);
                String[] elems =tf2.readLineElems(TextFile.tab);
                
                
                tf2.close();
            }
            tf.close();
        } catch (IOException e){
            e.printStackTrace();
        }
    }
    
}
