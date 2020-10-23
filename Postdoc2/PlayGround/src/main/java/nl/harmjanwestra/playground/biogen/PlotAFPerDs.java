package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.Set;

public class PlotAFPerDs {

	public static void main(String[] args) {
		String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\SNPQCLog-test.txt.gz-AFPerDS.txt.gz";
		String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\SNPQCLog-test.txt.gz-AFPerDS.png";
		PlotAFPerDs d = new PlotAFPerDs();
		String snpset = "U:\\IEUGWAS\\2020-06-25-2020-05-03-topHitsUniqueRSids-wgwascatalog-wALS-wAlzheimer-wMetaBrain-MetaBrain2dot1IDs.txt.gz";
		snpset = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-04-28-ALSGWASAndMSGWAS\\IEUGWAS\\2020-06-25-2020-05-03-topHitsUniqueRSids-wgwascatalog-wALS-wAlzheimer-wMetaBrain-MetaBrain2dot1IDs.txt.gz";
		try {
//            d.plot(input, output);
			output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\SNPQCLog-test.txt.gz-AFPerDS-allVsAll.png";
			d.plotAllvsAll(input, snpset, output);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}

	public void plot(String input, String output) throws IOException, DocumentException {
		TextFile tf = new TextFile(input, TextFile.R);
		Grid g = new Grid(300, 300, 10, 10, 100, 100);
		String[] header = tf.readLineElems(TextFile.tab);
		ScatterplotPanel[] p = new ScatterplotPanel[header.length - 3];
		for (int i = 0; i < p.length; i++) {
			p[i] = new ScatterplotPanel(1, 1);
			g.addPanel(p[i]);
			p[i].setTitle(header[i + 3]);
			p[i].setLabels("Ref", header[i + 3]);
			p[i].setDataRange(new Range(0, 0, 1, 1));
		}

		String[] elems = tf.readLineElems(TextFile.tab);
		int r = 0;
		while (elems != null) {
			// SNP     Alleles AF-A    EUR-AMPAD-MAYO-V2
			double ref = Double.parseDouble(elems[2]);
			if (!Double.isNaN(ref)) {
				for (int i = 3; i < elems.length; i++) {
					double v = Double.parseDouble(elems[i]);
					if (!Double.isNaN(v)) {
						p[i - 3].addData(ref, v);
//                        if (ref - v > 0.5 && i > 4) {
//                            System.out.println("got one");
//                        }
					}
				}
			}
			r++;
			if (r % 1000000 == 0) {
				System.out.println(r);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		g.draw(output);

	}

	public void plotAllvsAll(String input, String snpfilter, String output) throws IOException, DocumentException {

		TextFile tf2 = new TextFile(snpfilter, TextFile.R);
		Set<String> snpset = tf2.readAsSet(0, TextFile.tab);
		tf2.close();

		TextFile tf = new TextFile(input, TextFile.R);
		String[] header = tf.readLineElems(TextFile.tab);
		Grid g = new Grid(300, 300, header.length - 3, header.length - 3, 100, 100);
		g = new Grid(150, 150, 4, 4, 100, 100);

		ScatterplotPanel[][] p = new ScatterplotPanel[header.length - 3][header.length - 3];
//		for (int i = 0; i < p.length; i++) {
		int i = 0;
		for (int j = i + 1; j < p.length; j++) {
			ScatterplotPanel panel = new ScatterplotPanel(1, 1);
			p[i][j] = panel;
//			g.addPanel(panel, i, j);
			g.addPanel(panel);
			//panel.setTitle(header[i + 3]);
			panel.setLabels(header[i + 3], header[j + 3]);
			panel.setDataRange(new Range(0, 0, 1, 1));
		}


//		}

		String[] elems = tf.readLineElems(TextFile.tab);
		int r = 0;
		while (elems != null) {
			// SNP     Alleles AF-A    EUR-AMPAD-MAYO-V2
			if (snpset.contains(elems[0])) {
//				for (int i = 3; i < elems.length; i++) {
				i=3;
					double v = Double.parseDouble(elems[i]);

					if (!Double.isNaN(v)) {
						for (int j = i + 1; j < elems.length; j++) {
							double v2 = Double.parseDouble(elems[j]);
							if (!Double.isNaN(v2)) {
								p[i - 3][j - 3].addData(v, v2);
							}
						}

//                        if (ref - v > 0.5 && i > 4) {
//                            System.out.println("got one");
//                        }
					}
//				}


			}
			r++;
			if (r % 1000000 == 0) {
				System.out.println(r);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		g.draw(output);

	}
}
