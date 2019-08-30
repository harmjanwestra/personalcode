package nl.harmjanwestra.playground.cis;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.Iterator;

public class SimplifyQTLFile {

	public static void main(String[] args) {


		try {
			QTLTextFile f = new QTLTextFile(args[0], false);
			Iterator<EQTL> it = f.getEQtlIterator();

			TextFile out = new TextFile(args[1], TextFile.W);
			out.writeln("Chr\tPos\tRsId\tAlleles\tAssessed\tGene\tBeta\tSe\tN\tP\tFDR");
			EQTL e = it.next();
			int ctr = 0;
			while (it.hasNext()) {
				String betastr = e.getMetaBeta();
				String[] betaElems = betastr.split(" ");
				String beta = betaElems[0];
				String se = betaElems[1].replaceAll("\\(", "").replaceAll("\\)", "");

				Integer[] n = e.getDatasetsSamples();
				int sumn = 0;
				for (Integer i : n) {
					if (i != null) {
						sumn += i;
					}
				}
				out.writeln(e.getRsChr() + "\t" + e.getRsChrPos() + "\t" + e.getRsName() + "\t" + e.getAlleles() + "\t" + e.getAlleleAssessed()
						+ "\t" + e.getProbe() + "\t" + beta + "\t" + se + "\t" + sumn + "\t" + e.getPvalue() + "\t" + e.getFDR());
				ctr++;
				if (ctr % 1000000 == 0) {
					System.out.println(ctr + " lines processed.");
				}
				e = it.next();
			}
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
