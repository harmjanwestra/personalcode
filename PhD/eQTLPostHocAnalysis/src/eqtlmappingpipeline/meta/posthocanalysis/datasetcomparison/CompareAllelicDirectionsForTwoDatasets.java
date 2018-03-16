/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.datasetcomparison;

import java.io.IOException;
import java.util.zip.DataFormatException;
import umcg.genetica.io.trityper.bin.BinaryResultDataset;
import umcg.genetica.io.trityper.bin.BinaryResultProbe;
import umcg.genetica.io.trityper.bin.BinaryResultSNP;
import umcg.genetica.io.trityper.util.BaseAnnot;

/**
 *
 * @author harmjan
 */
public class CompareAllelicDirectionsForTwoDatasets {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {
	    String[] datasets = new String[2];
	    String[] platform = new String[2];
	    datasets[0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-24-EGCUT-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemove/";
	    platform[0] = "HT12v3";

	    //datasets[2][0] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-14-Groningen-BloodHT12-TRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    datasets[1] = "/Volumes/ADATA NH03/Marjolein/Results/2011-11-14-Groningen-BloodHT12-CIS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/";
	    platform[1] = "HT12v3";

	    BinaryResultDataset ds1 = new BinaryResultDataset(datasets[0], "Dataset", 0);
	    BinaryResultDataset ds2 = new BinaryResultDataset(datasets[1], "Dataset", 0);


	    BinaryResultSNP snp1 = ds1.getStringToSNP().get("rs903603");
	    if(snp1 == null){
		System.out.println("ERROR: snp not tested in Dataset1");
		System.exit(0);
	    }
	    Float[] z1 = ds1.readSNPZScores(snp1);
	    BinaryResultProbe probe1 = ds1.getStringToProbe().get("7570022");
	    if(probe1 == null){
		System.out.println("Probe not found for platform!");
		System.exit(0);
	    }

	    BinaryResultSNP snp2 = ds2.getStringToSNP().get("rs903603");
	    if(snp2 == null){
		System.out.println("ERROR: snp not tested in Dataset2");
		System.exit(0);
	    }
	    Float[] z2 = ds2.readSNPZScores(snp2);
	    BinaryResultProbe probe2 = ds2.getStringToProbe().get("7570022");

	    if(probe2 == null){
		System.out.println("Probe not found for platform!");
		System.exit(0);
	    }

	    System.out.println(
		    z1[probe1.getId()]+"\t"
		    +BaseAnnot.toString(snp1.getAssessedAllele())+"\t"
		    +BaseAnnot.toString(snp1.getAlleles()[0])+"/"+BaseAnnot.toString(snp1.getAlleles()[1])+"\t"
		    +BaseAnnot.toString(snp1.getMinorAllele())+"\t"
		    +snp1.getMaf());

	    System.out.println(
		    z2[probe2.getId()]+"\t"
		    +BaseAnnot.toString(snp2.getAssessedAllele())+"\t"
		    +BaseAnnot.toString(snp2.getAlleles()[0])+"/"+BaseAnnot.toString(snp2.getAlleles()[1])+"\t"
		    +BaseAnnot.toString(snp2.getMinorAllele())+"\t"
		    +snp2.getMaf());



	} catch (IOException e) {
	    e.printStackTrace();
	} catch (DataFormatException e) {
	    e.printStackTrace();
	}
    }
}
