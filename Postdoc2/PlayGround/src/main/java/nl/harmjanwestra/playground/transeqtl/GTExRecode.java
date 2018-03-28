package nl.harmjanwestra.playground.transeqtl;

import eqtlmappingpipeline.util.QTLFileSorter;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class GTExRecode {
	
	public static void main(String[] args) {
		GTExRecode gtr = new GTExRecode();
		try {
			gtr.run(args[0], args[1], args[2]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run(String in, String datasetname, String out) throws IOException {
		
		System.out.println("in:\t" + in);
		System.out.println("out:\t" + out);
		TextFile outf = new TextFile(out, TextFile.W);
		outf.writeln(QTLTextFile.header);
		TextFile tf = new TextFile(in, TextFile.R);
		
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			EQTL e = new EQTL();
			if (elems.length > 8) {
				// snp	gene	alleles	effect_allele	maf	beta	beta_se	pvalue	n_samples
				e.setRsName(elems[0]);
				e.setProbe(elems[1]);
				e.setAlleles(elems[2]);
				e.setAlleleAssessed(elems[3]);
				
				try {
					double beta = Double.parseDouble(elems[5]);
					double se = Double.parseDouble(elems[6]);
					e.setDatasets(new String[]{datasetname});
					e.setDatasetZScores(new Double[]{beta / se});
					e.setBeta(elems[5] + " (" + elems[6] + ")");
					e.setZscore(beta / se);
					e.setPvalue(Double.parseDouble(elems[7]));
					e.setDatasetsSamples(new Integer[]{Integer.parseInt(elems[8])});
					outf.writeln(e.toString());
				} catch (NumberFormatException ex) {
					
					System.out.println("Number format exception at line: " + Strings.concat(elems, Strings.tab));
				}
				
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		outf.close();
		
		// sort eQTL file
		QTLFileSorter sorter = new QTLFileSorter();
		String tmpout = out + "tmp.txt.gz";
		sorter.run(out, tmpout);
		Gpio.moveFile(tmpout, out);
		
		
	}
	
}
