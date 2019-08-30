package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

public class eQTLFileZscorePlot {


	public static void main(String[] args) {

//		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Trans\\eQTLsFDR0.05.txt.gz";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Trans\\eQTLsFDR0.05-ZScoreComp.png";


//		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Genes\\Primary\\eQTLProbesFDR0.05-ProbeLevel-ZScoreComp.png";

//		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\PRS\\eQTLsFDR0.05.txt.gz";
//		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\PRS\\eQTLsFDR0.05-ZScoreComp.png";

		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Transcripts\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-10-Results\\Cis-Transcripts\\eQTLProbesFDR0.05-ProbeLevel-ZScoreComp.png";
		eQTLFileZscorePlot z = new eQTLFileZscorePlot();
		try {
			z.run(in, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

	public void run(String in, String out) throws IOException, DocumentException {
		TextFile tf = new TextFile(in, TextFile.R);

		// 10 --> meta
		// 11 --> dsnames
		// 12 --> zscores
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		String[] ds = elems[11].split(";");
		int[] sizes = new int[ds.length];


		while (elems != null) {
			String[] dselem = elems[11].split(";");
			String[] size = elems[13].split(";");
			for (int i = 0; i < dselem.length; i++) {
				if (!dselem[i].equals("-")) {
					ds[i] = dselem[i];
					sizes[i] = Integer.parseInt(size[i]);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		int nrCols = ds.length / 4;
		Grid grid = new Grid(200, 200, 4, nrCols, 50, 50);
		ArrayList<ScatterplotPanel> panels = new ArrayList<>();
		Range range = new Range(-20, -20, 20, 20);
		for (int i = 0; i < ds.length; i++) {
			panels.add(new ScatterplotPanel(1, 1));
//			panels.get(i).setLabels("MetaAnalsys", ds[i]);
			panels.get(i).setTitle(ds[i] + " (n=" + sizes[i] + ")");
			panels.get(i).setDataRange(range);
			panels.get(i).setAlpha(0.2f);
		}


		tf.open();
		tf.readLine();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			Double z = Double.parseDouble(elems[10]);

			String[] dselem = elems[12].split(";");
			for (int i = 0; i < dselem.length; i++) {
				if (!dselem[i].equals("-")) {
					double z2 = Double.parseDouble(dselem[i]);
					panels.get(i).addData(z, z2);
				}
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		for (ScatterplotPanel p : panels) {
			grid.addPanel(p);
		}

		grid.draw(out);

	}
}
