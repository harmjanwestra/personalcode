package nl.harmjanwestra.playground.conv;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.Gpio;
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
				if (this.pos > o.pos) {
					return 1;
				} else if (this.pos < o.pos) {
					return -1;
				} else {
					return 0;
				}
			} else if (chr > o.chr) {
				return 1;
			} else {
				return -1;
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


	public static void main(String[] args) {

		TriTyperMerge m = new TriTyperMerge();
		if (args.length < 2) {
			System.out.println("Usage: inputfile.txt outdir [fraction of datasets snps must be present in] [listofsnps]");
			System.out.println("Input file format (tab-sep):\n" +
					"Location\tDatasetName\t[SNPs.txt.gz loc]\t[SNPMappings.txt.gz loc]\t[list of samples to exclude]");
		} else {
			try {
				double fractionpresent = 1d;
				String snplist = null;
				if (args.length >= 3) {
					fractionpresent = Double.parseDouble(args[2]);
				}
				if (args.length >= 4) {
					snplist = args[3];
				}
				m.runFromList(args[0], args[1], fractionpresent, snplist);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	private void runFromList(String arg, String outdir, double fractionpresent, String snplist) throws IOException {

		System.out.println("TriTyper Merger. ");
		System.out.println("Input: " + arg);
		System.out.println("Output: " + outdir);
		System.out.println("% datasets in which SNPs must be present: " + fractionpresent);
		System.out.println("SNPlist: " + snplist);
		ArrayList<Dataset> datasets = new ArrayList<>();

		TextFile tf = new TextFile(arg, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		int ctr = 0;
		boolean cont = true;
		while (elems != null) {

			String loc = elems[0];
			if (!Gpio.exists(loc)) {
				System.out.println("Could not find: " + loc);
				cont = false;
			}
			String name = "Dataset_" + ctr;
			if (elems.length > 1) {
				name = elems[1];

			}

			String snps = null;
			String snpmap = null;
			String samplesExclude = null;

			if (elems.length > 2 && elems[2].length() > 0) {
				snps = elems[2];
				if (snps.length() >= 0 && !Gpio.exists(snps)) {
					System.out.println("Error: could not find " + snps);
					cont = false;
				}
			}
			if (elems.length > 3 && elems[3].length() > 0) {
				snpmap = elems[3];
				if (!Gpio.exists(snpmap)) {
					System.out.println("Error: could not find " + snpmap);
					cont = false;
				}
			}

			if (elems.length > 4 && elems[4].length() > 0) {
				samplesExclude = elems[4];
				if (!Gpio.exists(samplesExclude)) {
					System.out.println("Error: could not find: " + samplesExclude);
					cont = false;
				}
			}

			ctr++;
			Dataset ds = new Dataset();
			ds.loc = loc;
			ds.name = name;
			if (snps == null) {
				snps = loc + "SNPs.txt.gz";
			}
			ds.snps = snps;
			if (snpmap == null) {
				snpmap = loc + "SNPMappings.txt.gz";
			}
			ds.snpmap = snpmap;
			datasets.add(ds);
			ds.samplesExclude = samplesExclude;
			System.out.println(ds);

			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();

		if (!cont) {
			System.exit(-1);
		}

		run(datasets, outdir, fractionpresent, snplist);
	}

	public class Dataset {
		public String samplesExclude;
		public Boolean[] sampleInclude;
		public ArrayList<String> loadedIndividuals;
		String name;
		String loc;
		String snpmap;
		String snps;

		@Override
		public String toString() {
			return "Dataset{" +
					"samplesExclude='" + samplesExclude + '\'' +
					", name='" + name + '\'' +
					", loc='" + loc + '\'' +
					", snpmap='" + snpmap + '\'' +
					", snps='" + snps + '\'' +
					'}';
		}
	}

	public void run(ArrayList<Dataset> datasets, String outdir, double fractionpresent, String snplist) throws IOException {
		// read SNPs


		SNP[] snps = readAndFilterSNPs(datasets, fractionpresent, snplist);

		// read individuals
		ArrayList<String> individuals = new ArrayList<>();
		TextFile indo = new TextFile(outdir + "Individuals.txt", TextFile.W);
		TextFile pheno = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);

		boolean cont = true;
		for (int f = 0; f < datasets.size(); f++) {

			HashSet<String> excludeSamples = null;
			if (datasets.get(f).samplesExclude != null) {
				excludeSamples = readHash(datasets.get(f).samplesExclude);
			}

			ArrayList<String> inds = readAsList(datasets.get(f).loc + "Individuals.txt");

			Boolean[] sampleInclude = new Boolean[inds.size()];

			ArrayList<String> loadedInds = new ArrayList<>();
			for (int i = 0; i < inds.size(); i++) {
				String s = inds.get(i);
				if (excludeSamples == null || !excludeSamples.contains(s)) {
					sampleInclude[i] = true;
					individuals.add(datasets.get(f).name + "-" + s);
					indo.writeln(datasets.get(f).name + "-" + s);
					pheno.writeln(datasets.get(f).name + "-" + s + "\tunknown\tinclude\tunknown");
					loadedInds.add(s);
				} else {
					sampleInclude[i] = false;
				}
			}
			datasets.get(f).loadedIndividuals = loadedInds;
			if (loadedInds.isEmpty()) {
				System.out.println("Error: no samples remain in dataset " + datasets.get(f).name);
				cont = false;
			}
			datasets.get(f).sampleInclude = sampleInclude;

		}
		pheno.close();
		indo.close();

		if (!cont) {
			System.exit(-1);
		}

		// now we have the dimensions of the new matrix
		TriTyperGenotypeData[] genotypeDatasets = new TriTyperGenotypeData[datasets.size()];
		SNPLoader[] loaders = new SNPLoader[datasets.size()];
		IntStream.range(0, datasets.size()).parallel().forEach(i -> {
			try {
				genotypeDatasets[i] = new TriTyperGenotypeData();
				genotypeDatasets[i].load(datasets.get(i).loc, datasets.get(i).snpmap, datasets.get(i).snps);
				genotypeDatasets[i].setIsIncluded(datasets.get(i).sampleInclude);
				loaders[i] = genotypeDatasets[i].createSNPLoader(100);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});

		byte[] genotypesA1 = new byte[individuals.size()];
		byte[] genotypesA2 = new byte[individuals.size()];
		byte[] dosages = new byte[individuals.size()];
		int[] startint = new int[datasets.size()];

		System.out.println("Dataset 0: " + 0);
		for (int i = 1; i < startint.length; i++) {
			startint[i] = startint[i - 1] + datasets.get(i - 1).loadedIndividuals.size();
			System.out.println("Dataset " + i + ": " + startint[i]);
		}


		BinaryFile bfgt = new BinaryFile(outdir + "GenotypeMatrix.dat", BinaryFile.W);
		BinaryFile bfds = new BinaryFile(outdir + "ImputedDosageMatrix.dat", BinaryFile.W);
		TextFile snpo = new TextFile(outdir + "SNPs.txt.gz", TextFile.W);
		TextFile snpmo = new TextFile(outdir + "SNPMappings.txt.gz", TextFile.W);

		TextFile tfsnpqc = new TextFile(outdir + "SNPQC.txt.gz", TextFile.W);

		byte nullbyte = (byte) 0;
		ProgressBar pb = new ProgressBar(snps.length, "Merrrrrrging!");
		for (int snpNr = 0; snpNr < snps.length; snpNr++) {

			if (snps[snpNr].chr > 0) {
				Arrays.fill(genotypesA1, nullbyte);
				Arrays.fill(genotypesA2, nullbyte);
				Arrays.fill(dosages, nullbyte);

				String snpname = snps[snpNr].name;

				umcg.genetica.io.trityper.SNP referenceSNP = null;
				String refalleles = null;
				byte[] refallelesb = null;
				String refminor = null;

				for (int datasetId = 0; datasetId < datasets.size(); datasetId++) {
					TriTyperGenotypeData gtDs = genotypeDatasets[datasetId];
					Dataset ds = datasets.get(datasetId);
					int nrIndsInDs = ds.loadedIndividuals.size();
					int snpDatasetId = gtDs.getSnpToSNPId().get(snpname);
					int starti = startint[datasetId];
					Boolean flipAllele = false;
					String alleles = null;
					String minor = null;
					if (snpDatasetId > -1) {
						umcg.genetica.io.trityper.SNP snpobj = gtDs.getSNPObject(snpDatasetId);
						loaders[datasetId].loadGenotypes(snpobj);

						if (loaders[datasetId].hasDosageInformation()) {
							loaders[datasetId].loadDosage(snpobj);
						}

						if (referenceSNP == null) {
							// set reference snp
							referenceSNP = snpobj;

							// copy alleles and dosage (if present)
							System.arraycopy(snpobj.getAllele1(), 0, genotypesA1, starti, nrIndsInDs);
							System.arraycopy(snpobj.getAllele2(), 0, genotypesA2, starti, nrIndsInDs);
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
							// compare alleles with reference
							alleles = BaseAnnot.getAllelesDescription(snpobj.getAlleles());
							minor = BaseAnnot.toString(snpobj.getMinorAllele());

							flipAllele = BaseAnnot.flipalleles(refalleles, refminor, alleles, minor);


							if (flipAllele != null) {
								if (flipAllele) {
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
//						System.out.println(datasetId + "\t" + snpname + "\t" + snpDatasetId + "\t" + refalleles + "," + refminor + "\t" + alleles + "\t" + minor + "\t" + flipAllele);
						tfsnpqc.writeln(datasetId + "\t" + snpname + "\t" + snpDatasetId + "\t" + refalleles + "," + refminor + "\t" + alleles + "\t" + minor + "\t" + flipAllele);
						snpobj.clearGenotypes();
					}
				}
//				System.out.println();

				if (referenceSNP != null) {
					// write result to disk
					bfgt.write(genotypesA1);
					bfgt.write(genotypesA2);
					bfds.write(dosages);
					snpo.writeln(snpname);
					snpmo.writeln(snps[snpNr].chr + "\t" + snps[snpNr].pos + "\t" + snpname + "\t" + BaseAnnot.toString(refallelesb[0]) + "," + BaseAnnot.toString(refallelesb[1]));
				} else {
					System.out.println("Could not find " + snpname + " in any dataset?");
				}
			}
			if (snpNr % 5000 == 0) {
				pb.set(snpNr);
			}

		}
		pb.close();

		for (int d = 0; d < loaders.length; d++) {
			loaders[d].close();
		}

		tfsnpqc.close();
		// close writers
		bfgt.close();
		bfds.close();
		snpo.close();
		snpmo.close();

	}

	private HashSet<String> readHash(String samplesExclude) throws IOException {

		TextFile tf = new TextFile(samplesExclude, TextFile.R);
		HashSet<String> output = new HashSet<String>();
		output.addAll(tf.readAsArrayList());
		tf.close();

		return output;
	}

	private SNP[] readAndFilterSNPs(ArrayList<Dataset> folders, double fractionpresent, String snplist) throws IOException {

		HashSet<String> allowedSNPs = null;
		if (snplist != null) {
			TextFile tf = new TextFile(snplist, TextFile.R);

			String ln = tf.readLine();
			while (ln != null) {
				allowedSNPs.add(ln);
				ln = tf.readLine();
			}
			tf.close();
		}


		SNP[] snps = null;
		{
			ArrayList<SNP> allsnpTMP = new ArrayList<>();

			// list snps
			HashSet<String> allsnpstr = new HashSet<>();
			IntStream.range(0, folders.size()).parallel().forEach(i -> {

				try {
					ArrayList<String> list = readAsList(folders.get(i).snps);
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
				if (allowedSNPs == null || allowedSNPs.contains(s)) {
					SNP obj = new SNP(0, 0, s);
					allsnpTMP.add(obj);
					snpmap.put(s, obj);
				}
			}

			// read annotation
			IntStream.range(0, folders.size()).parallel().forEach(i -> {
				try {
					TextFile tf = new TextFile(folders.get(i).snpmap, TextFile.R);
					String[] elems = tf.readLineElems(TextFile.tab);
					while (elems != null) {

						SNP snpobj = snpmap.get(elems[2]);
						if (snpobj != null) {
							try {
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
							} catch (NumberFormatException e) {
								System.out.println("Could not parse " + elems[0] + " or " + elems[1]);
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
			System.out.println(allsnpTMP.size() + " SNPs initially found");
			ArrayList<SNP> allsnpstmp2 = new ArrayList<>();
			for (SNP obj : allsnpTMP) {
				if (obj.chr > 0 && ((double) obj.nr / folders.size()) >= fractionpresent) {
					allsnpstmp2.add(obj);
				}
			}
			snps = allsnpstmp2.toArray(new SNP[0]);
			Arrays.parallelSort(snps);
			System.out.println(snps.length + " SNPs found in total (located on autosome, present in >1 dataset.");
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
