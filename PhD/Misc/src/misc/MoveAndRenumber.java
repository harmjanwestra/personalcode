/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.File;
import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class MoveAndRenumber {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            String dir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/RandomlySNPSelectionsForFuncAnalysisSet2/";
            for (int perm = 10; perm < 20; perm++) {
                String newDir = "/Volumes/iSnackHD/Data/Projects/2012-eQTLMetaAnalysis/MetaAnalysisResults/traitassociatedsnps/RandomlySNPSelectionsForFuncAnalysisSet2/Set" + perm + "/";
                Gpio.createDir(newDir);

//Set-0.txt-WithProxies-List.txt-SNPNexus.txt
                Gpio.copyFile(new File(dir + "Set-" + (perm - 10) + ".txt"), new File(newDir + "Set-" + (perm) + ".txt"));
                Gpio.copyFile(new File(dir + "Set-" + (perm - 10) + "-Dist.txt"), new File(newDir + "Set-" + (perm) + "-Dist.txt"));
                Gpio.copyFile(new File(dir + "Set-" + (perm - 10) + ".txt-WithProxies-List.txt"), new File(newDir + "Set-" + (perm) + ".txt-WithProxies-List.txt"));
                Gpio.copyFile(new File(dir + "Set-" + (perm - 10) + ".txt-WithProxies-List.txt-SNPNexus.txt"), new File(newDir + "Set-" + (perm) + ".txt-WithProxies-List.txt-SNPNexus.txt"));
                Gpio.copyFile(new File(dir + "Set-" + (perm - 10) + ".txt-WithProxies.txt"), new File(newDir + "Set-" + (perm) + ".txt-WithProxies.txt"));

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
