package nl.harmjanwestra.playground.transeqtl;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class MixAndMatch {
	
	public static void main(String[] args) {
//		String file1 = "D:\\Work\\TMP\\eqtlinvalidate\\top200fxgtex.txt";
//		String file1 = "D:\\Work\\TMP\\eqtlinvalidate\\2018-03-27-GTExReplication.txt";
//		String file2 = "D:\\Work\\TMP\\eqtlinvalidate\\eQTLsFDR-Significant-0.05-1mbcismapping.txt";
//		String outf = "D:\\Work\\TMP\\eqtlinvalidate\\comp2.txt";
//
		MixAndMatch m = new MixAndMatch();
//		try {
//			m.combineEPRSWithCisAndTransFX(file1, file2, outf);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
		String gtexreplfile = "D:\\Work\\TMP\\eqtlinvalidate\\2018-03-27-GTExReplication.txt";
		String mappercfile = "D:\\Work\\TMP\\eqtlinvalidate\\2018-03-27-eQTLCrossMap-5mb-35bp-t2.txt";
		String crossmapfile = "D:\\Work\\TMP\\eqtlinvalidate\\scripts-wtih-data\\replication_Whole_Blood_0_formatted_crossmap.txt";
		String out = "D:\\Work\\TMP\\gtexrepl\\2018-03-27-GTExReplication-crossmap-5mb.txt";
		
		try {
			m.run2(gtexreplfile, crossmapfile, mappercfile, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	public void run(String efile, String efile2, String outf) throws IOException {
		TextFile tf = new TextFile(efile, TextFile.R);
		
		HashMap<String, String> strGene = new HashMap<String, String>();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			strGene.put(elems[0] + "_" + elems[1], Strings.concat(elems, Strings.tab));
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile out = new TextFile(outf, TextFile.W);
		TextFile tf2 = new TextFile(efile2, TextFile.R);
		
		out.writeln(tf2.readLine() + "\t" + strGene.get("snp_gene"));
		elems = tf2.readLineElems(TextFile.tab);
		while (elems != null) {
			String query = elems[0] + "_" + elems[3];
			
			if (strGene.containsKey(query)) {
				String output = Strings.concat(elems, Strings.tab) + "\t" + strGene.get(query);
				out.writeln(output);
			}
			
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		out.close();
		
		
	}
	
	public void run2(String gtexrepfile, String crossmapfile, String mappercentagefile, String outf) throws IOException {
		
		HashMap<String, String> mapperc = new HashMap<String, String>();
		HashMap<String, String> mapclosest = new HashMap<String, String>();
		TextFile tf3 = new TextFile(mappercentagefile, TextFile.R);
		String[] elems = tf3.readLineElems(TextFile.tab);
		while (elems != null) {
			if (elems.length > 5) {
				String snp = elems[0];
				String gene = elems[1];
				String map = elems[4];
				String closest = elems[5];
				String eqtl = snp + "_" + gene;
				
				mapclosest.put(eqtl, closest);
				mapperc.put(eqtl, map);
			}
			elems = tf3.readLineElems(TextFile.tab);
		}
		tf3.close();
		
		
		TextFile tf = new TextFile(crossmapfile, TextFile.R);
		HashMap<String, String> crossmap = new HashMap<String, String>();
		elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			
			if (elems.length > 22) {
				String snp = elems[0];
				String gene = elems[13];
				String map = elems[22];
				String eqtl = snp + "_" + gene;
				crossmap.put(eqtl, map);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile out = new TextFile(outf, TextFile.W);
		TextFile tf2 = new TextFile(gtexrepfile, TextFile.R);
		
		out.writeln(tf2.readLine() + "\tcrossmap\tmapperc");
		elems = tf2.readLineElems(TextFile.tab);
		
		ArrayList<Double> x1 = new ArrayList<Double>();
		ArrayList<Double> x2 = new ArrayList<Double>();
		ArrayList<Double> y1 = new ArrayList<Double>();
		ArrayList<Double> y2 = new ArrayList<Double>();
		
		while (elems != null) {
			String query = elems[3] + "_" + elems[0];
			
			String output = Strings.concat(elems, Strings.tab) + "\t" + crossmap.get(query) + "\t" + mapperc.get(query) + "\t" + mapclosest.get(query);
			out.writeln(output);
			
			String dsteststr = elems[elems.length - 3];
			String dssigstr = elems[elems.length - 2];
			
			String perccrossmap = mapperc.get(query);
			if (perccrossmap != null) {
				
				double nrds = Double.parseDouble(dsteststr);
				if (nrds >= 10) {
					double x = Double.parseDouble(dssigstr) / nrds;
					
					if(x>1){
						System.out.println("?");
					}
					
					if (crossmap.get(query) != null && crossmap.get(query).startsWith("ENSG")) {
						// match by ashis
						x1.add(x);
						y1.add(Double.parseDouble(perccrossmap));
					} else {
						x2.add(x);
						y2.add(Double.parseDouble(perccrossmap));
					}
				}
			}
			
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		out.close();
		
		Grid g = new Grid(500, 500, 1, 1, 50, 50);
		
		
		double[][] x = new double[2][];
		double[][] y = new double[2][];
		
		x[0] = Primitives.toPrimitiveArr(x2);
		y[0] = Primitives.toPrimitiveArr(y2);
		x[1] = Primitives.toPrimitiveArr(x1);
		y[1] = Primitives.toPrimitiveArr(y1);
//
		ScatterplotPanel p = new ScatterplotPanel(1, 1);
		p.setData(x, y);
		p.setDataRange(new Range(0, 0, 1, 1));
		g.addPanel(p);
		try {
			g.draw(outf + "-plot.pdf");
			System.out.println(outf + "-plot.pdf");
		} catch (DocumentException e) {
			e.printStackTrace();
		}
		
		
	}
}
