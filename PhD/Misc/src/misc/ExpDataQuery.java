/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import eqtlmappingpipeline.util.ExpressionDataQuery;
import java.io.IOException;

/**
 *
 * @author harmjan
 */
public class ExpDataQuery {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {

            String probeAnnotationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-HT12v3.txt";
            String snpList = "/Volumes/iSnackHD/Data/Projects/TonuEsko/2013-01-16-HeightSNPs.txt";
            String snpAnnotation = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2011-09-27-SNPMappings-dbSNP130.txt.gz";
            int eqtlwindow = 1000000;
            String outfile = "/Volumes/iSnackHD/Data/Projects/TonuEsko/2013-01-16-HeightSNPsProbesWithin1MB.txt";
            ExpressionDataQuery.getProbesAroundSNPs(probeAnnotationFile, snpList, snpAnnotation, eqtlwindow, outfile);
            
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
