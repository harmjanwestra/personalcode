package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.BaseAnnot;

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

	private byte convertDosageToByte(double dosageValue) {
		int dosageInt = (int) Math.round(dosageValue * 100d);
		byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
		return dosageByte;
	}

	public void run(String[] folders, String[] names, String outdir) throws IOException {
		// read SNPs
		SNP[] snps = readAndFilterSNPs(folders);

		// read individuals
		ArrayList<String> individuals = new ArrayList<>();
		TextFile indo = new TextFile(outdir + "Individuals.txt", TextFile.W);
		TextFile pheno = new TextFile(outdir + "Individuals.txt", TextFile.W);

		for (int f = 0; f < folders.length; f++) {
			ArrayList<String> inds = readAsList(folders[f] + "Individuals.txt");
			for (String s : inds) {
				individuals.add(names[f] + "-" + s);
				indo.writeln(names[f] + "-" + s);
				pheno.writeln(names[f] + "-" + s );
			}
		}
		pheno.close();
		indo.close();

		// now we have the dimensions of the new matrix
		TriTyperGenotypeData[] datasets = new TriTyperGenotypeData[folders.length];
		SNPLoader[] loaders = new SNPLoader[folders.length];
		IntStream.range(0, folders.length).parallel().forEach(i -> {
			try {
				datasets[i] = new TriTyperGenotypeData(folders[i]);
				loaders[i] = datasets[i].createSNPLoader(1);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});


		byte[] genotypesA1 = new byte[individuals.size()];
		byte[] genotypesA2 = new byte[individuals.size()];
		byte[] dosages = new byte[individuals.size()];
		int[] startint = new int[datasets.length];
		for (int i = 1; i < startint.length; i++) {
			startint[i] = datasets[i - 1].getIndividuals().length;
		}


		BinaryFile bfgt = new BinaryFile(outdir + "GenotypeMatrix.dat", BinaryFile.W);
		BinaryFile bfds = new BinaryFile(outdir + "ImputedDosageMatrix.dat", BinaryFile.W);
		TextFile snpo = new TextFile(outdir + "SNPs.txt.gz", TextFile.W);
		TextFile snpmo = new TextFile(outdir + "SNPMappings.txt.gz", TextFile.W);

		for (int i = 0; i < snps.length; i++) {

			String snpname = snps[i].name;
			for (int j = 0; j < individuals.size(); j++) {
				genotypesA1[j] = -1;
				genotypesA2[j] = -1;
				dosages[j] = -1;
			}


			umcg.genetica.io.trityper.SNP referenceSNP = null;
			String refalleles = null;
			byte[] refallelesb = null;
			String refminor = null;

			for (int j = 0; j < datasets.length; j++) {
				int id = datasets[i].getSnpToSNPId().get(snpname);
				int starti = startint[j];
				if (id > -1) {
					umcg.genetica.io.trityper.SNP snpobj = datasets[i].getSNPObject(id);
					loaders[i].loadGenotypes(snpobj);
					if (loaders[i].hasDosageInformation()) {
						loaders[i].loadDosage(snpobj);
					}

					if (referenceSNP == null) {
						// set reference snp
						referenceSNP = snpobj;

						// copy alleles and dosage (if present)
						System.arraycopy(snpobj.getAllele1(), 0, genotypesA1, starti, datasets[i].getIndividuals().length);
						System.arraycopy(snpobj.getAllele2(), 0, genotypesA2, starti, datasets[i].getIndividuals().length);
						refalleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
						refallelesb = snpobj.getAlleles();
						refminor = BaseAnnot.toString(snpobj.getMinorAllele());

						if (snpobj.hasDosageInformation()) {
							double[] dsdosge = snpobj.getDosageValues();

							for (int d = 0; d < dsdosge.length; d++) {
								dosages[starti + d] = convertDosageToByte(dsdosge[d]);
							}
						} else {
							// convert genotypes do dosage
							byte[] genotypes = snpobj.getGenotypes();
							for (int d = 0; d < genotypes.length; d++) {
								dosages[starti + d] = convertDosageToByte(genotypes[d]);
							}

						}

					} else {
						// compare alleles
						String alleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
						String minor = BaseAnnot.toString(snpobj.getMinorAllele());

						Boolean flipallele = BaseAnnot.flipalleles(refalleles, refminor, alleles, minor);
						if (flipallele != null) {
							if (flipallele) {
								// use genotypes to define alleles
								byte[] genotypes = snpobj.getGenotypes();
								double[] dsdosge = snpobj.getDosageValues();
								for (int d = 0; d < genotypes.length; d++) {
									if (genotypes[d] > -1) {
										if (genotypes[d] == 0) {
											genotypesA1[starti + d] = refallelesb[1];
											genotypesA2[starti + d] = refallelesb[1];
										} else if (genotypes[d] == 1) {
											genotypesA1[starti + d] = refallelesb[0];
											genotypesA2[starti + d] = refallelesb[1];
										} else {
											genotypesA1[starti + d] = refallelesb[0];
											genotypesA2[starti + d] = refallelesb[0];
										}
									}
									if (snpobj.hasDosageInformation()) {
										if (dsdosge[d] != -1) {
											dosages[starti + d] = convertDosageToByte(Math.abs(dsdosge[d] - 2));
										}
									} else {
										if (genotypes[d] != -1) {
											dosages[starti + d] = convertDosageToByte(Math.abs(genotypes[d] - 2));
										}
									}
								}
							} else {
								// need to take into account allele encoding...
								byte[] genotypes = snpobj.getGenotypes();
								double[] dsdosge = snpobj.getDosageValues();
								for (int d = 0; d < genotypes.length; d++) {
									if (genotypes[d] > -1) {
										if (genotypes[d] == 0) {
											genotypesA1[starti + d] = refallelesb[0];
											genotypesA2[starti + d] = refallelesb[0];
										} else if (genotypes[d] == 1) {
											genotypesA1[starti + d] = refallelesb[0];
											genotypesA2[starti + d] = refallelesb[1];
										} else {
											genotypesA1[starti + d] = refallelesb[1];
											genotypesA2[starti + d] = refallelesb[1];
										}
									}
									if (snpobj.hasDosageInformation()) {
										if (dsdosge[d] != -1) {
											dosages[starti + d] = convertDosageToByte(dsdosge[d]);
										}
									} else {
										if (genotypes[d] != -1) {
											dosages[starti + d] = convertDosageToByte(genotypes[d]);
										}
									}
								}
							}
						}
					}
				}
			}

			// write result to disk
			bfgt.write(genotypesA1);
			bfgt.write(genotypesA2);
			bfds.write(dosages);
			snpo.writeln(snpname);
			snpmo.writeln(snps[i].chr + "\t" + snps[i].pos + "\t" + snpname + "\t" + BaseAnnot.toString(refallelesb[0]) + "," + BaseAnnot.toString(refallelesb[1]));

		}

		// close writers
		bfgt.close();
		bfds.close();
		snpo.close();

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
