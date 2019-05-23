package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class TriTyperMerge {


	public class SNP implements Comparable<SNP> {
		int chr;
		int pos;
		int nr;
		String name;

		public SNP(int i, int i1, String s) {
			this.chr = i;
			this.pos = i1;
			this.name = s;
		}


		@Override
		public int compareTo(SNP o) {
			if (o.equals(this)) {
				return 0;
			}
			if (chr == o.chr) {
				if (this.chr > o.chr) {
					return -1;
				} else if (this.chr < o.chr) {
					return 1;
				} else {
					return 0;
				}
			} else if (chr > o.chr) {
				return -1;
			} else {
				return 1;
			}
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			SNP snp = (SNP) o;
			return chr == snp.chr &&
					pos == snp.pos &&
					name.equals(snp.name);
		}

		@Override
		public int hashCode() {
			return Objects.hash(chr, pos, name);
		}
	}

	public void run(String[] folders, String[] names) throws IOException {
		// read SNPs
		SNP[] snps = readAndFilterSNPs(folders);

		// read individuals
		ArrayList<String> individuals = new ArrayList<>();
		for (int f = 0; f < folders.length; f++) {
			ArrayList<String> inds = readAsList(folders[f] + "Individuals.txt");
			for (String s : inds) {
				individuals.add(names[f] + "-" + s);
			}
		}


	}

	private SNP[] readAndFilterSNPs(String[] folders) {
		SNP[] snps = null;
		{

			ArrayList<SNP> allsnpTMP = new ArrayList<>();

			// list snps
			HashSet<String> allsnpstr = new HashSet<>();
			IntStream.range(0, folders.length).parallel().forEach(i -> {

				try {
					ArrayList<String> list = readAsList(folders[i] + "SNPs.txt.gz");
					synchronized (allsnpTMP) {
						allsnpstr.addAll(list);
					}
				} catch (IOException e) {
					e.printStackTrace();
				}

			});

			// index snps
			HashMap<String, SNP> snpmap = new HashMap<String, SNP>();
			for (String s : allsnpstr) {
				SNP obj = new SNP(0, 0, s);
				allsnpTMP.add(obj);
				snpmap.put(s, obj);
			}

			// read annotation
			IntStream.range(0, folders.length).parallel().forEach(i -> {
				try {
					TextFile tf = new TextFile(folders[i] + "SNPMappings.txt.gz", TextFile.R);
					String[] elems = tf.readLineElems(TextFile.tab);
					while (elems != null) {

						SNP snpobj = snpmap.get(elems[2]);
						if (snpobj != null) {
							Integer chr = Integer.parseInt(elems[0]);
							Integer pos = Integer.parseInt(elems[1]);
							if (chr > 0 && pos > 0) {
								synchronized (snpobj) {
									if (snpobj.chr > 0) {
										if (snpobj.chr != chr.intValue() && snpobj.pos != pos.intValue()) {
											snpobj.chr = -1;
											snpobj.pos = -1;
										} else {
											snpobj.nr++;
										}
									} else if (snpobj.chr != -1) {
										snpobj.chr = chr;
										snpobj.pos = pos;
										snpobj.nr++;
									}
								}
							}
						}
						elems = tf.readLineElems(TextFile.tab);
					}
					tf.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			});

			// strip variants with missing annotation
			ArrayList<SNP> allsnpstmp2 = new ArrayList<>();
			for (SNP obj : allsnpTMP) {
				if (obj.chr > 0) {
					allsnpstmp2.add(obj);
				}
			}
			snps = allsnpstmp2.toArray(new SNP[0]);
			Arrays.parallelSort(snps);
		}
		return snps;
	}

	private ArrayList<String> readAsList(String s) throws IOException {
		System.out.println("Reading: " + s);
		ArrayList<String> output = new ArrayList<>();
		TextFile tf = new TextFile(s, TextFile.R);
		String ln = tf.readLine();
		while (ln != null) {
			output.add(ln);
			ln = tf.readLine();
		}
		tf.close();
		System.out.println("File: " + s + " has " + output.size() + " lines.");
		return output;
	}
}
