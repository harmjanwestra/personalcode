package nl.harmjanwestra.playground.meqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class EQTMFlipper {
	
	
	public static void main(String[] args) {
		
		
		double fdrthreshold = 0.0;
//		String infile = "D:\\Work\\Projects\\2017-11-eQTLMeta\\eqtm\\eQTLsFDR0.05-SNPLevel.txt.gz";
//		String infile = "D:\\Work\\Projects\\2018-eQTMPredict\\bonder\\2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt";
//		String outfile = "D:\\Work\\Projects\\2017-11-eQTLMeta\\eqtm\\eQTLsFDR" + fdrthreshold + "-SNPLevel-flipped.txt.gz";
//		String outfile = "D:\\Work\\Projects\\2018-eQTMPredict\\bonder\\2015_09_02_cis_eQTMsFDR0.0-CpGLevel-filter.txt";
		
		String infile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\eQTLsFDR0.05.txt";
		
		
		EQTMFlipper f = new EQTMFlipper();
		try {
			String outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-gt0.0-flipped.txt";
//			f.combineEPRSWithCisAndTransFX(infile, outfile, fdrthreshold, MODE.ABOVE);
			
			outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-et0.0-flipped.txt";
//			f.combineEPRSWithCisAndTransFX(infile, outfile, fdrthreshold, MODE.BELOW);
			
			infile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\eQTLsFDR.txt.gz";
			outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-gt0.05-flipped.txt.gz";
			fdrthreshold = 0.05;
			System.out.println(fdrthreshold);
//			f.flipnfilter(infile, outfile, fdrthreshold, MODE.ABOVE);
			
			outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-gt0.5-flipped.txt.gz";
			fdrthreshold = 0.5;
			System.out.println(fdrthreshold);
			
//			f.flipnfilter(infile, outfile, fdrthreshold, MODE.ABOVE);
			
			infile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\eQTLsFDR.txt.gz";
			outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-lt0.5-flipped.txt.gz";
			fdrthreshold = 0.5;
			f.flipnfilter(infile, outfile, fdrthreshold, MODE.BELOW);

			outfile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-02-eQTMPredict\\2017-12-09-eQTMs-250k\\2017-12-09-eQTLsFDR-lt0.05-flipped.txt.gz";
			fdrthreshold = 0.05;
			f.flipnfilter(infile, outfile, fdrthreshold, MODE.BELOW);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	enum MODE {
		ABOVE,
		BELOW
	}
	
	public void flipnfilter(String infile, String outfile, double fdrthreshold, MODE m) throws IOException {
		QTLTextFile f = new QTLTextFile(infile, QTLTextFile.R);
		Iterator<EQTL> it = f.getEQtlIterator();
		
		
		
		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln(QTLTextFile.header);
		
		int nrpos = 0;
		int nrneg = 0;
		
		
		HashSet<String> visitedCGs = new HashSet<String>();
		while(it.hasNext()){
			EQTL me = it.next();
			boolean flip = false;
			double fdr = me.getFDR();
			boolean include = false;
			if (m.equals(MODE.ABOVE)) {
				if (fdr > fdrthreshold) {
					include = true;
				}
			} else {
				if (fdr <= fdrthreshold) {
					include = true;
				}
			}
			if (include) {

				if (!me.getAlleleAssessed().equals("C")) {
					flip = true;
				}
				
				if (flip) {
					double z = me.getZscore();
					z *= -1;
					me.setZscore(z);
				}
				
				me.setAlleleAssessed("C");
				me.setAlleles("C/T");
				
				out.writeln(me.toString());
			}
		}
		out.close();
		
		
	}
	
	public void run(String infile, String outfile, double fdrthreshold, MODE m) throws IOException {
		QTLTextFile f = new QTLTextFile(infile, QTLTextFile.R);
		ArrayList<EQTL> meqtl = f.readList();
		f.close();
		
		TextFile out = new TextFile(outfile, TextFile.W);
		out.writeln(QTLTextFile.header);
		
		int nrpos = 0;
		int nrneg = 0;
		
		ArrayList<EQTL> split = new ArrayList<>();
		HashSet<String> visitedCGs = new HashSet<String>();
		for (EQTL me : meqtl) {
			boolean flip = false;
			double fdr = me.getFDR();
			boolean include = false;
			if (m.equals(MODE.ABOVE)) {
				if (fdr > fdrthreshold) {
					include = true;
				}
			} else {
				if (fdr <= fdrthreshold) {
					include = true;
				}
			}
			if (include) {
				if (!me.getAlleleAssessed().equals("C")) {
					flip = true;
				}
				
				if (flip) {
					double z = me.getZscore();
					z *= -1;
					me.setZscore(z);
				}
				
				me.setAlleleAssessed("C");
				me.setAlleles("C/T");
				
				double z = me.getZscore();
				
				
				if (!visitedCGs.contains(me.getRsName())) {
					split.add(me);
					visitedCGs.add(me.getRsName());
					if (z < 0) {
						nrneg++;
					} else {
						nrpos++;
					}
				}
				
				out.writeln(me.toString());
			}
		}
		out.close();
		
		
		System.out.println(nrneg + " negative");
		System.out.println(nrpos + " positive");
		
		// get top negative effects
		int nroutputneg = 0;
		int nroutputpos = 0;
		
		int nrTooutput = Math.min(nrneg, nrpos);
		TextFile out2 = new TextFile(outfile + "split.txt", TextFile.W);
		out2.writeln(QTLTextFile.header);
		HashMap<String, ArrayList<EQTL>> cgToQTL = new HashMap<>();
		
		for (EQTL me : split) {
			
			double z = me.getZscore();
			if (z < 0) {
				if (nroutputneg < nrTooutput) {
					out2.writeln(me.toString());
					ArrayList<EQTL> mes = cgToQTL.get(me.getRsName());
					if (mes == null) {
						mes = new ArrayList<>();
					}
					mes.add(me);
					cgToQTL.put(me.getRsName(), mes);
					nroutputneg++;
				}
			} else {
				if (nroutputpos < nrTooutput) {
					ArrayList<EQTL> mes = cgToQTL.get(me.getRsName());
					if (mes == null) {
						mes = new ArrayList<>();
					}
					mes.add(me);
					cgToQTL.put(me.getRsName(), mes);
					out2.writeln(me.toString());
					nroutputpos++;
				}
			}
		}
		out2.close();
		
		int nrcgsoppositedirs = 0;
		for (String cg : cgToQTL.keySet()) {
			int prevdir = 0;
			ArrayList<EQTL> mes = cgToQTL.get(cg);
			
			if (mes.size() > 1) {
				boolean positive = mes.get(0).getZscore() > 0;
				boolean samedir = true;
				for (int i = 1; i < mes.size(); i++) {
					boolean alsopositive = mes.get(i).getZscore() > 0;
					if (positive && !alsopositive) {
						samedir = false;
					} else if (!positive && alsopositive) {
						samedir = false;
					}
				}
				if (!samedir) {
					nrcgsoppositedirs++;
				}
			}
			
			
		}
		System.out.println(nrcgsoppositedirs + " out of " + cgToQTL.size() + " have QTMs with opposite dirs");
		
	}
	
}
