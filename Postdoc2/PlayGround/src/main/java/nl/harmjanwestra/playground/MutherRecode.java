package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class MutherRecode {
	public static void main(String[] args) {
		
		MutherRecode r = new MutherRecode();
		if (args[0].equals("eprs")) {
			r.muthereprs(args[1], args[2]);
		} else {
			r.muthereqtl(args[1], args[2]);
		}
		
	}
	
	public void muthereprs(String inf, String outf) {
		try {
			TextFile tf = new TextFile(inf, TextFile.R);
			TextFile out = new TextFile(outf, TextFile.W);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			/*
			PValue
			SNPName
			SNPChr
			SNPChrPos
			ProbeName
			ProbeChr
			ProbeCenterChrPos
			CisTrans
			SNPType
			AlleleAssessed
			OverallZScore
			DatasetNames
			DatasetZ
			DatasetN
			IncludedDatasetsCorrelationCoefficient
			Meta-Beta (SE)
			Beta (SE)
			FoldChange
			FDR
			 */
			
			String header = "PValue" +
					"\tSNPName" +
					"\tSNPChr" +
					"\tSNPChrPos" +
					"\tProbeName" +
					"\tProbeChr" +
					"\tProbeCenterChrPos" +
					"\tCisTrans" +
					"\tSNPType" +
					"\tAlleleAssessed" +
					"\tOverallZScore" +
					"\tDatasetNames" +
					"\tDatasetZ" +
					"\tDatasetN" +
					"\tIncludedDatasetsCorrelationCoefficient" +
					"\tMeta-Beta (SE)" +
					"\tBeta (SE)" +
					"\tFoldChange" +
					"\tFDR";
			out.writeln(header);
			while (elems != null) {
				
				String[] elemsout = new String[19];
				for (int q = 0; q < elemsout.length; q++) {
					elemsout[q] = "-";
				}
				
				// GWAS    Gene    Pval    Beta    SE(Beta)
				
				elemsout[0] = elems[2];
				elemsout[1] = elems[0];
				elemsout[4] = elems[1];
				elemsout[QTLTextFile.DATASETNAMES] = "Muther";
				elemsout[QTLTextFile.DATASETSIZE] = "" + 535;
				
				double beta = Double.parseDouble(elems[3]);
				Double p = Double.parseDouble(elems[2]);
				double z = ZScores.zToP(p);
				if (beta < 0) {
					z *= -1;
				}
				
				elemsout[QTLTextFile.METAZ] = "" + z;
				elemsout[QTLTextFile.DATASETZSCORE] = "" + z;
				elemsout[QTLTextFile.ASESSEDALLELE - 1] = "T/C"; // assessed
				elemsout[QTLTextFile.ASESSEDALLELE] = "C";

//				if (!elems[13].equals("NA")) {
				out.writeln(Strings.concat(elemsout, Strings.tab));
//				}
//
				elems = tf.readLineElems(TextFile.tab);
				
			}
			
			out.close();
			tf.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void muthereqtl(String inf, String outf) {
		try {
			TextFile tf = new TextFile(inf, TextFile.R);
			TextFile out = new TextFile(outf, TextFile.W);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			/*
			PValue
			SNPName
			SNPChr
			SNPChrPos
			ProbeName
			ProbeChr
			ProbeCenterChrPos
			CisTrans
			SNPType
			AlleleAssessed
			OverallZScore
			DatasetNames
			DatasetZ
			DatasetN
			IncludedDatasetsCorrelationCoefficient
			Meta-Beta (SE)
			Beta (SE)
			FoldChange
			FDR
			 */
			
			String header = "PValue" +
					"\tSNPName" +
					"\tSNPChr" +
					"\tSNPChrPos" +
					"\tProbeName" +
					"\tProbeChr" +
					"\tProbeCenterChrPos" +
					"\tCisTrans" +
					"\tSNPType" +
					"\tAlleleAssessed" +
					"\tOverallZScore" +
					"\tDatasetNames" +
					"\tDatasetZ" +
					"\tDatasetN" +
					"\tIncludedDatasetsCorrelationCoefficient" +
					"\tMeta-Beta (SE)" +
					"\tBeta (SE)" +
					"\tFoldChange" +
					"\tFDR";
			out.writeln(header);
			while (elems != null) {
				
				String[] elemsout = new String[19];
				for (int q = 0; q < elemsout.length; q++) {
					elemsout[q] = "-";
				}
				elemsout[0] = elems[12];
				elemsout[1] = elems[1];
				elemsout[4] = elems[4];
				elemsout[QTLTextFile.DATASETNAMES] = "Muther";
				elemsout[QTLTextFile.DATASETSIZE] = "" + 535;
				
				elemsout[QTLTextFile.METAZ] = elems[13];
				elemsout[QTLTextFile.DATASETZSCORE] = elems[13];
				elemsout[QTLTextFile.ASESSEDALLELE - 1] = elems[8]; // assessed
				elemsout[QTLTextFile.ASESSEDALLELE] = elems[9];
				
				if (!elems[13].equals("NA")) {
					out.writeln(Strings.concat(elemsout, Strings.tab));
				}
				
				elems = tf.readLineElems(TextFile.tab);
				
			}
			
			out.close();
			tf.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
