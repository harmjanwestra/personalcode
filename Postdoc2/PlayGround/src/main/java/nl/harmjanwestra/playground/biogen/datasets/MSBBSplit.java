package nl.harmjanwestra.playground.biogen.datasets;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class MSBBSplit {
	
	public static void main(String[] args) {
		
		String wtob = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\MSBB-BrainBankIDToWGS.txt";
		String btor = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\MSBB-RNASeqToBrainBankID.txt";
		
		
		try {
			
			
			String infile = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\MSBB_Mapping-hjw.txt";
			String allRNAseq = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-10-09-AMPAD\\mssbrnaseqsampleIds.txt";
			
			{
				TextFile tf = new TextFile(allRNAseq, TextFile.R);
				String ln = tf.readLine();
				HashSet<String> samples = new HashSet<String>();
				while (ln != null) {
					samples.add(ln);
					ln = tf.readLine();
				}
				tf.close();
				
				// link to dna
				HashMap<String, String> rnaseqToBrainBank = new HashMap<String, String>();
				TextFile tf1 = new TextFile(btor, TextFile.R);
				String[] elems1 = tf1.readLineElems(TextFile.tab);
				while (elems1 != null) {
					if (samples.contains(elems1[0])) {
						rnaseqToBrainBank.put(elems1[0], elems1[1]);
					}
					elems1 = tf1.readLineElems(TextFile.tab);
				}
				tf1.close();
				
				HashMap<String, String> brainbanktodna = new HashMap<String, String>();
				TextFile tf2 = new TextFile(wtob, TextFile.R);
				String[] elems2 = tf2.readLineElems(TextFile.tab);
				while (elems2 != null) {
					brainbanktodna.put(elems2[0], elems2[1]);
					elems2 = tf2.readLineElems(TextFile.tab);
				}
				tf2.close();
				
				// link
				TextFile linkout = new TextFile(infile, TextFile.W);
				for (String sample : samples) {
					String bbid = rnaseqToBrainBank.get(sample);
					if (bbid != null) {
						String dna = brainbanktodna.get(bbid);
						if (dna != null) {
							linkout.writeln(dna + "\t" + sample);
						}
						
					}
				}
				linkout.close();
				
				// FILTER EUROPEANS
				String eur = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-Europeans.txt";
				TextFile tf3 = new TextFile(eur, TextFile.R);
				HashSet<String> europeans = new HashSet<String>();
				String ln3 = tf3.readLine();
				while (ln3 != null) {
					europeans.add(ln3);
					ln3 = tf3.readLine();
				}
				tf3.close();
				
				linkout = new TextFile(infile + "-EUR.txt", TextFile.W);
				for (String sample : samples) {
					String bbid = rnaseqToBrainBank.get(sample);
					if (bbid != null) {
						String dna = brainbanktodna.get(bbid);
						if (dna != null) {
							if (europeans.contains(dna)) {
								linkout.writeln(dna + "\t" + sample);
							}
						}
						
					}
				}
				linkout.close();
				
				// split per tissue
				String sampleToBm = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\MSSBSampleToBM.txt";
				TextFile tf4 = new TextFile(sampleToBm, TextFile.R);
				HashMap<String, String> sampleToBM = new HashMap<String, String>();
				String[] elems4 = tf4.readLineElems(TextFile.tab);
				while (elems4 != null) {
					
					sampleToBM.put(elems4[0], elems4[1]);
					elems4 = tf4.readLineElems(TextFile.tab);
				}
				tf4.close();
				
				String[] regions = new String[]{
						"BM10", "BM22", "BM36", "BM44"
				};
				
				for (String region : regions) {
					
					linkout = new TextFile(infile + "-EUR-" + region + ".txt", TextFile.W);
					for (String sample : samples) {
						String sampleregion = sampleToBM.get(sample);
						if (sampleregion.equals(region)) {
							String bbid = rnaseqToBrainBank.get(sample);
							if (bbid != null) {
								String dna = brainbanktodna.get(bbid);
								if (dna != null) {
									if (europeans.contains(dna)) {
										linkout.writeln(dna + "\t" + sample);
									}
								}
								
							}
						}
					}
					linkout.close();
				}
				
				
			}
			System.exit(-1);


//            {
//                HashMap<String, String> wtobh = new HashMap<>();
//                TextFile tf = new TextFile(wtob, TextFile.R);
//                String[] elems = tf.readLineElems(TextFile.tab);
//                while (elems != null) {
//                    if (elems.length > 1) {
//                        wtobh.put(elems[0], elems[1]);
//                        System.out.println(elems[1] + "\t" + elems[0]);
//                    }
//                    elems = tf.readLineElems(TextFile.tab);
//                }
//                tf.close();
//                System.out.println(wtobh.size() + " read");
//                TextFile out = new TextFile(infile, TextFile.W);
//                TextFile tf2 = new TextFile(btor, TextFile.R);
//                elems = tf2.readLineElems(TextFile.tab);
//                while (elems != null) {
//                    if (elems.length > 1) {
//                        String rna = elems[0];
//                        String bbid = elems[1];
//
//                        String dna = wtobh.get(bbid);
//                        System.out.println(rna + "\t" + bbid + "\t" + dna);
//                        if (dna != null) {
//                            out.writeln(dna + "\t" + rna);
//                        }
//                    }
//                    elems = tf2.readLineElems(TextFile.tab);
//                }
//                tf2.close();
//                out.close();
////                System.exit(0);
//            }
			
			String sampleToBm = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\MSSBSampleToBM.txt";
			String eur = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\genotypepca\\AMPAD-Europeans.txt";
			
			String out = "D:\\Sync\\SyncThing\\Postdoc2\\2018-04-BioGen\\data\\2018-08-02-AMPAd\\linksToRNA\\MSBB_final.txt";
			
			TextFile tf3 = new TextFile(eur, TextFile.R);
			HashSet<String> eursamples = new HashSet<>();
			String[] e3 = tf3.readLineElems(TextFile.tab);
			while (e3 != null) {
				System.out.println(e3[0]);
				eursamples.add(e3[0]);
				e3 = tf3.readLineElems(TextFile.tab);
			}
			tf3.close();
			
			TextFile tf = new TextFile(sampleToBm, TextFile.R);
			HashMap<String, String> sampleToBmH = new HashMap<String, String>();
			HashSet<String> bms = new HashSet<String>();
			
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {
				sampleToBmH.put(elems[0], elems[1]);
				bms.add(elems[1]);
				elems = tf.readLineElems(TextFile.tab);
			}
			tf.close();
			
			for (String bm : bms) {
				System.out.println("loading BM " + bm);
				TextFile tf2 = new TextFile(infile, TextFile.R);
				String[] elems2 = tf2.readLineElems(TextFile.tab);
				TextFile tfo = new TextFile(out + "-" + bm + ".txt", TextFile.W);
				HashSet<String> visited = new HashSet<>();
				while (elems2 != null) {
//                    System.out.println(elems2[0]);
					if (elems2.length >= 2) {
						String rna = elems2[1];
						String dna = elems2[0];
						
						if (!visited.contains(dna) && eursamples.contains(dna)) {
//							System.out.println(dna + "\t" + eursamples.contains(dna));
							String rnabm = sampleToBmH.get(rna);
							if (rnabm == null) {
								System.out.println(rna + "\tno bm");
							} else if (rnabm.equals(bm)) {
								System.out.println(dna + "\t" + eursamples.contains(dna) + "\t" + rna);
								tfo.writeln(dna + "\t" + rna);
							}
							visited.add(dna);
						}
					}
					elems2 = tf2.readLineElems(TextFile.tab);
				}
				tf.close();
				tfo.close();
			}
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
}
