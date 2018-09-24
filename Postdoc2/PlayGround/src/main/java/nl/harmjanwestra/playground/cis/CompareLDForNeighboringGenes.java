package nl.harmjanwestra.playground.cis;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class CompareLDForNeighboringGenes {
	
	
	public static void main(String[] args) {

//		if (args.length < 3) {
//			System.out.println("Usage: ref file out");
//		} else {
		
		
		String[] ludefile = new String[]{
				"D:\\ld\\EstimatedLDBetweenPairsOfSNPs.txt",
				"D:\\ld\\EstimatedLDBetweenPairsOfSNPs-alldatasets.txt",
				"D:\\ld\\EstimatedLDBetweenPairsOfSNPs-WOCOEXPANDNONEUR.txt",
			
		};
//		String hjfile = "d:\\ld\\ENSG00000002933-2018-07-11-QueryGenes-InnerProduct.-eur-txt.gz";
		String hjfile = "d:\\ld\\ENSG00000002933-2018-07-11-QueryGenes-Pearson.txt.gz";
		String out = "D:\\ld\\ldcomp3.txt";
		
		
		CompareLDForNeighboringGenes g = new CompareLDForNeighboringGenes();
		try {
//			g.gtfToProbeAnnotationFile(args[0], args[1], args[2]);
			g.run(ludefile, hjfile, out);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}

//		}
		
	}
	
	public void run(String[] reference, String file, String outf) throws IOException, DocumentException {
		
		Grid g = new Grid(500, 500, 1, reference.length, 100, 100);
		
		for (String rds : reference) {
			HashMap<String, Double> comboToId = new HashMap<String, Double>();
			
			TextFile tf = new TextFile(rds, TextFile.R);
			tf.readLine();
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				String snp1 = elems[0];
				String snp2 = elems[1];
				
				if (!snp1.equals(snp2)) {
					comboToId.put(snp1 + "-" + snp2, Double.parseDouble(elems[4]));
				}
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			
			ArrayList<Double> x = new ArrayList<Double>();
			ArrayList<Double> y = new ArrayList<Double>();
			TextFile tf2 = new TextFile(file, TextFile.R);
			TextFile out = new TextFile(outf, TextFile.W);
			tf2.readLine();
			elems = tf2.readLineElems(TextFile.tab);
			while (elems != null) {
				String snp1 = elems[0];
				String snp2 = elems[1];
				
				if (!snp1.equals(snp2)) {
					String combo = snp1 + "-" + snp2;
					if (comboToId.containsKey(combo)) {
						double r = Double.parseDouble(elems[6]);
						double ref = comboToId.get(combo);
						x.add(ref);
						y.add(r);
						out.writeln(combo + "\t" + ref + "\t" + r);
					}
				}
				elems = tf2.readLineElems(TextFile.tab);
			}
			tf2.close();
			out.close();
			
			
			PearsonsCorrelation c = new PearsonsCorrelation();
			double correl = c.correlation(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
			
			ScatterplotPanel p = new ScatterplotPanel(1, 1);
			p.setTitle("Correlation of R-squareds: " + correl);
			p.setDataRange(new Range(0, 0, 1, 1));
			p.setLabels(rds, "1KG");
			p.setData(Primitives.toPrimitiveArr(x), Primitives.toPrimitiveArr(y));
			p.setPlotElems(true, false);
			g.addPanel(p);
		}
		
		g.draw(outf + ".png");
		
		
	}
	
}
