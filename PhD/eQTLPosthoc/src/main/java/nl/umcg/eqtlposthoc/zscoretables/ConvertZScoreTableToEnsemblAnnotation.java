/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.zscoretables;

import java.io.IOException;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class ConvertZScoreTableToEnsemblAnnotation {

    static String probeListFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-05-02-ProbesThatHaveMappingInBothAnnotationFilesAndHT12v3.txt";
    static String ensemblAnnotationFile = "/Volumes/Data2/MarjoleinHomeAccount/marjolein/AnnotationFiles/2012-04-23-IlluminaAll96PercentIdentity-ProbeAnnotation-ProbesWithWrongMappingLengthFilteredOut-OnlyEnsemblAnnotation.txt";

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	for (int i = 8; i < 11; i++) {
	    String in = null;
	    if (i == 0) {
//		in = "/Volumes/BackupDisk/MetaAnalysisFinal/cistrans/2012-05-30-Cistrans/metazscoretable";
	    } else {
		in = "/Volumes/iSnackHD/MetaAnalysisFinal/cistrans/2012-05-30-Cistrans/metazscoretable-Permutation"+i;
//		in = "/Data/tmp2/metazscoretable-Permutation" + i;
	    }
	    try {
		String infile = in + ".txt.gz";

		ZScoreTableManipulation manipulator = new ZScoreTableManipulation();
		String filename = infile;

		String finalFileName = filename + "-filtered.txt" + "-ens.txt" + "-collapsed.txt";
//		if (!Gpio.exists(finalFileName)) {
//		    if (!Gpio.exists(filename + "-filtered.txt")) {
			manipulator.filterZScoreTableForSetOfProbes(probeListFile, infile);
//		    }

		    filename += "-filtered.txt";
//		    if (!Gpio.exists(filename + "-ens.txt")) {
			manipulator.convertProbeNumbersToEnsemblId(ensemblAnnotationFile, filename);
//		    }
		    filename += "-ens.txt";

//		    if (!Gpio.exists(filename + "-collapsed.txt")) {
			manipulator.collapseDuplicateProbes(filename);
//		    }
		    filename += "-collapsed.txt";
//		}
	    } catch (IOException e) {

		e.printStackTrace();
	    }
	}

    }
}
