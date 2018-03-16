/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.eqtlfile;

import eqtlmappingpipeline.util.eQTLFileCompare;
import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harm-jan
 */
public class eQTLFileBatchCompare {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        eQTLFileCompare eqfc = new eQTLFileCompare();
                
        String transFile = "d:\\SkyDrive\\latesteQTLs\\transFDR0.05.txt.gz";
        try{
            int[] arrr = new int[]{25,30,35,40};
            for(int q: arrr){
                String fairfaxfile = "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\bcellmetaFiltered\\"+q+"PC\\eQTLs.txt.gz";
                String outdir = "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\bcellmetaCompare\\"+q+"PC\\";
                Gpio.createDir(outdir);
                eqfc.compareOverlapAndZScoreDirectionTwoEQTLFiles(transFile, fairfaxfile, outdir+"CompNoFDR", false);
            }
            
            arrr = new int[]{20,30,35,40,50};
            for(int q: arrr){
                String fairfaxfile  ="d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometaFiltered\\"+q+"PC\\eQTLs.txt.gz";    
                String outdir = "d:\\SkyDrive\\latesteqtls\\Replication\\FaifFaxProbeSubsetNorm\\monometaCompare\\"+q+"\\";
                Gpio.createDir(outdir);
                eqfc.compareOverlapAndZScoreDirectionTwoEQTLFiles(transFile, fairfaxfile, outdir+"CompNoFDR", false);
            }
        } catch (IOException e){
            e.printStackTrace();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}
