package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import nl.harmjanwestra.utilities.graphics.Grid;
import nl.harmjanwestra.utilities.graphics.Range;
import nl.harmjanwestra.utilities.graphics.panels.ScatterplotPanel;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.BaseAnnot;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Replicate {
	
	
	public static void main(String[] args) {
		
		String significantCisEqtlgen = "";
		String allTransMetabrain = "";
		String allCisMetabrain = "";
		
		Replicate r = new Replicate();
		
		// metabrain to eqtlgen
		try {
			String significantTransMetabrain = "D:\\biogen\\metabrainfx\\trans\\trans-eQTLsFDR0.05.txt.gz";
			String allTransEqtlgen = "D:\\biogen\\eqtlgenfx\\trans\\2018-09-04-trans-eQTLsFDR-CohortInfoRemoved.txt.gz";
			String out = "D:\\biogen\\comparisons\\metabrainToEqtlgen\\transComp";
			r.metabrainToEqtlGen(significantTransMetabrain, allTransEqtlgen, out);
			
			String significantCisMetabrain = "D:\\biogen\\metabrainfx\\cis\\eQTLProbesFDR0.05-ProbeLevel.txt.gz";
			String allCisEqtlgen = "D:\\Sync\\SyncThing\\Postdoc2\\2018-05-eQTLMeta\\data\\2018-05-22-assoc\\cis\\2018-01-31-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved.txt.gz";
			out = "D:\\biogen\\comparisons\\metabrainToEqtlgen\\cisComp";
//			r.metabrainToEqtlGen(significantCisMetabrain, allCisEqtlgen, out);
			
			// eqtlgen to metabrain
			String significantTransEqtlgen = "D:\\biogen\\eqtlgenfx\\trans\\2018-04-03-trans-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CohortInfoRemoved-b38-snpupdate-geneupdate.txt.gz";
			String replicationTransMetabrain = "D:\\biogen\\eqtlgenfx-metabrainrepl\\10pcs\\eQTLsFDR.txt.gz";
			out = "D:\\biogen\\comparisons\\metabrainToEqtlgen\\transComp-eqtlgenrepl";
			r.metabrainToEqtlGen(replicationTransMetabrain, significantTransEqtlgen, out);
			
			String significantTransHighConfidenceEqtlgen = "D:\\biogen\\eqtlgenfx\\trans\\2018-04-03-trans-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05-CohortInfoRemoved-b38-snpupdate-geneupdate.txt.gz";
			String replicationTransHighConfidenceMetabrain = "D:\\biogen\\eqtlgenfx-metabrainrepl\\10pcs-highconfidence\\eQTLsFDR.txt.gz";
			out = "D:\\biogen\\comparisons\\metabrainToEqtlgen\\transComp-eqtlgenrepl-highconfidence";
			r.metabrainToEqtlGen(replicationTransHighConfidenceMetabrain, significantTransHighConfidenceEqtlgen, out);
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (DocumentException e) {
			e.printStackTrace();
		}
	}
	
	class EQTL {
		String alleles;
		String assessed;
		Double z;
		String snp;
		String gene;
		Double fdr;
		Double p;
		
		@Override
		public String toString() {
			return p + "\t" + snp + "\t" + gene + "\t" + alleles + "\t" + assessed + "\t" + z + "\t" + fdr;
		}
	}
	
	// significant metabrain into eqtlgen full
	public void metabrainToEqtlGen(String metabrain, String eqtlgen, String out) throws IOException, DocumentException {
		
		HashMap<String, EQTL> stringToEQTL = new HashMap<String, EQTL>();
		TextFile tf = new TextFile(metabrain, TextFile.R);
		tf.readLine();
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			EQTL e = lineAsEQTL(elems);
			stringToEQTL.put(e.snp + "/" + e.gene, e);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		System.out.println(stringToEQTL.size() + " eqtls from: " + metabrain);
		
		ArrayList<Pair<EQTL, EQTL>> pairs = new ArrayList<Pair<EQTL, EQTL>>();
		TextFile tf2 = new TextFile(eqtlgen, TextFile.R);
		tf2.readLine();
		elems = tf2.readLineElems(TextFile.tab);
		int ctr = 0;
		
		while (elems != null) {
			String gene = elems[4];
			String[] geneelems = gene.split("\\.");
			gene = geneelems[0];
			String combo = elems[1] + "/" + gene;
			EQTL other = stringToEQTL.get(combo);
			if (other != null) {
				
				EQTL e = lineAsEQTL(elems);
				pairs.add(new Pair<>(other, e));
				
				Boolean flip = BaseAnnot.flipalleles(other.alleles, other.assessed, e.alleles, e.assessed);
				if (flip != null && flip) {
					e.z *= -1;
				}
				
			}
			ctr++;
			if (ctr % 1000000 == 0) {
				System.out.println(pairs.size() + " matched, " + ctr + " lines read.");
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		
		// make plots
		Grid g = new Grid(500, 500, 1, 2, 100, 100);
		ArrayList<Double> x1 = new ArrayList<>();
		ArrayList<Double> x2 = new ArrayList<>();
		ArrayList<Double> y1 = new ArrayList<>();
		ArrayList<Double> y2 = new ArrayList<>();
		int x1samedir = 0;
		int x2samedir = 0;
		for (Pair<EQTL, EQTL> p : pairs) {
			EQTL e1 = p.getLeft();
			EQTL e2 = p.getRight();
			
			boolean samedir = true;
			if (e1.z > 0 && e2.z < 0 || e1.z < 0 && e2.z > 0) {
				samedir = false;
			}
			x1.add(e1.z);
			y1.add(e2.z);
			if (samedir) {
				x1samedir++;
			}
			
			if (e1.fdr < 0.05 && e2.fdr < 0.05) {
				x2.add(e1.z);
				y2.add(e2.z);
				if (samedir) {
					x2samedir++;
				}
			}
		}
		
		
		ScatterplotPanel p1 = new ScatterplotPanel(1, 1);
		p1.setData(Primitives.toPrimitiveArr(x1), Primitives.toPrimitiveArr(y1));
		
		p1.setLabels("Metabrain (Z-score)", "eQTLGen (Z-score)");
		p1.setTitle("Shared: " + x1.size() + ", same direction: " + x1samedir);
		p1.setPlotElems(true, false);
		
		Range r = new Range(Primitives.toPrimitiveArr(x1), Primitives.toPrimitiveArr(y1));
		r.round();
		
		ScatterplotPanel p2 = new ScatterplotPanel(1, 1);
		p2.setData(Primitives.toPrimitiveArr(x2), Primitives.toPrimitiveArr(y2));
		p2.setTitle("FDR < 0.05 in both - Shared: " + x2.size() + ", same direction: " + x2samedir);
		p2.setLabels("Metabrain (Z-score)", "eQTLGen (Z-score)");
		p2.setPlotElems(true, false);
		p2.setDataRange(r);
		
		g.addPanel(p1);
		g.addPanel(p2);
		
		g.draw(out + "-Plot.pdf");
		
		// write uniquely significant hits from metabrain to disk
		TextFile outf = new TextFile(out + "-UniqueInMetabrain.txt", TextFile.W);
		for (Pair<EQTL, EQTL> p : pairs) {
			if (p.getLeft().fdr < 0.05 && p.getRight().fdr > 0.05) {
				outf.writeln(p.getLeft().toString());
			}
		}
		outf.close();
		
	}
	
	private EQTL lineAsEQTL(String[] elems) {
		Double p = Double.parseDouble(elems[0]);
		String snp = elems[1];
		String gene = elems[4];
		String[] geneelems = gene.split("\\.");
		gene = geneelems[0];
		String alleles = elems[8];
		String assessed = elems[9];
		Double z = Double.parseDouble(elems[10]);
		Double fdr = Double.parseDouble(elems[elems.length - 1]);
		
		EQTL e = new EQTL();
		e.p = p;
		e.alleles = alleles;
		e.assessed = assessed;
		e.z = z;
		e.snp = snp;
		e.gene = gene;
		e.fdr = fdr;
		return e;
	}
	
	
}
