package nl.harmjanwestra.playground.methylation;

import umcg.genetica.enums.Gender;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

public class GEOGatherAnnotation {

	public static void main(String[] args) {
		String in = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\GSEIds.txt";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2019-methylation\\2019-05-10-SampleAnnotation.txt";
		GEOGatherAnnotation r = new GEOGatherAnnotation();
		try {
			r.run(in, out);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void run(String in, String out) throws IOException {
		TextFile tf = new TextFile(in, TextFile.R);

		String ln = tf.readLine();
		int ctr = 0;
		TextFile tfout = new TextFile(out, TextFile.W);
		String header = "GSE\tGSM\tgender\tage\tchildage\tchildgender\tcellline\ttissue\tsource\tdiseasestate\ttaxonomy";
		tfout.writeln(header);
		while (ln != null) {

			System.out.println(ctr + "\t" + ln);
			process(ln, tfout);
			ln = tf.readLine();

			ctr++;
		}
		tfout.close();
		tf.close();
//
//		String ln = "GSE62219";
//
//		TextFile tfout = new TextFile(out + ln + ".txt", TextFile.W);
//		process(ln, tfout);
//		ln = tf.readLine();
//		tfout.close();


	}

	public static final String ANSI_BLACK = "\u001B[30m";
	public static final String ANSI_RED = "\u001B[31m";
	public static final String ANSI_GREEN = "\u001B[32m";
	public static final String ANSI_YELLOW = "\u001B[33m";
	public static final String ANSI_BLUE = "\u001B[34m";
	public static final String ANSI_PURPLE = "\u001B[35m";
	public static final String ANSI_CYAN = "\u001B[36m";
	public static final String ANSI_WHITE = "\u001B[37m";
	public static final String ANSI_RESET = "\u001B[0m";

	enum SAMPLETYPE {
		CELLINE,
		SOMATIC,
		UNKNOWN
	}

	private void process(String id, TextFile tfout) throws IOException {

		String platformStr = "GPL13534";
		String gse = id;
		String gsennn = id.substring(0, id.length() - 3) + "nnn";
		// req =

		String loc = "https://ftp.ncbi.nlm.nih.gov/geo/series/" + gsennn + "/" + gse + "/matrix/" + gse + "_series_matrix.txt.gz";

		URL url = new URL(loc);
		HttpURLConnection con = null;
		int len = -1;
		try {
			con = (HttpURLConnection) url.openConnection();
			con.setRequestProperty("Accept-Encoding", "gzip");
			len = con.getContentLength();
		} catch (FileNotFoundException e) {
			System.out.println("File not found: " + loc);
			con.disconnect();
			con = null;
		}
		if (con == null) {
			loc = "https://ftp.ncbi.nlm.nih.gov/geo/series/" + gsennn + "/" + gse + "/matrix/" + gse + "-" + platformStr + "_series_matrix.txt.gz";
			System.out.println("Trying: ");
			try {
				url = new URL(loc);
				con = (HttpURLConnection) url.openConnection();
				con.setRequestProperty("Accept-Encoding", "gzip");
				len = con.getContentLength();
			} catch (FileNotFoundException e) {
				System.out.println("File not found: " + loc);
				con.disconnect();
				con = null;
			}
		}


		if (con != null && len > 0) {
			ArrayList<String> accids = new ArrayList<>();
			ArrayList<Gender> genders = new ArrayList<>();
			ArrayList<String> ages = new ArrayList<>();
			ArrayList<String> childages = new ArrayList<>();
			ArrayList<Gender> childgenders = new ArrayList<>();
			ArrayList<SAMPLETYPE> celllines = new ArrayList<>();
			ArrayList<String> tissues = new ArrayList<>();
			ArrayList<String> sources = new ArrayList<>();
			ArrayList<String> diseasestate = new ArrayList<>();
			ArrayList<String> taxids = new ArrayList<>();


			System.out.println(gse + "\t" + gsennn + "\tLength : " + len);
			BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(con.getInputStream())));

			String ln = reader.readLine();
			while (ln != null) {
//			System.out.println(ln);
				if (ln.startsWith("!") || ln.trim().length() == 0) {
					if (ln.trim().length() > 0) {

						ln = ln.replaceAll("\"", "");


						String[] elems = ln.toLowerCase().split("\t");

						if (ln.startsWith("!Series_summary")) {
//							System.out.println(ln);
						} else if (ln.startsWith("!Sample_source_name_ch1")) {
							boolean gender = false;
							boolean tissue = false;
							boolean disease = false;
							boolean age = false;
							tissue = true;
							int nrtissuesset = 0;
							for (int i = 1; i < elems.length; i++) {
								sources.set(i - 1, elems[1]);
								nrtissuesset++;
							}
							if (tissue) {

//								System.out.println(ln);
								System.out.println(ANSI_PURPLE + Strings.concat(sources, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_PURPLE + "Found sample source for " + nrtissuesset + " samples" + ANSI_RESET);
							}

						} else if (ln.startsWith("!Sample_characteristics_ch1")) {

							boolean gender = false;
							boolean tissue = false;
							boolean disease = false;
							boolean age = false;

							int nrgendersset = 0;
							int nrtissuesset = 0;
							int nrdiseaseset = 0;
							int nrageset = 0;
							for (int i = 1; i < elems.length; i++) {
								if (elems[i].contains("cell line")) {
									celllines.set(i - 1, SAMPLETYPE.CELLINE);
								}
								if (elems[i].startsWith("gender:") || elems[i].startsWith("sex:") || elems[i].startsWith("cell sex:")) {
									gender = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										if (comps[1].startsWith("male") || comps[1].startsWith("m")) {
											genders.set(i - 1, Gender.MALE);
											nrgendersset++;
										} else if (comps[1].startsWith("female") || comps[1].startsWith("f")) {
											genders.add(i - 1, Gender.FEMALE);
											nrgendersset++;
										}
									}


								} else if (elems[i].startsWith("child gender:")
										|| elems[i].startsWith("child sex:")) {
									gender = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										if (comps[1].startsWith("male") || comps[1].startsWith("m")) {
											childgenders.set(i - 1, Gender.MALE);
											nrgendersset++;
										} else if (comps[1].startsWith("female") || comps[1].startsWith("f")) {
											childgenders.add(i - 1, Gender.FEMALE);
											nrgendersset++;
										}
									}
								} else if (elems[i].startsWith("tissue:")
										|| elems[i].startsWith("cell type:")
										|| elems[i].startsWith("tissue_type:")
										|| elems[i].startsWith("histology:")
										|| elems[i].startsWith("organ:")
										|| elems[i].startsWith("cell:")
										|| elems[i].startsWith("site:")
										|| elems[i].startsWith("cell line:")
										|| elems[i].startsWith("cell lineage:")) {
									tissue = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										tissues.set(i - 1, comps[1]);
										nrtissuesset++;
									}

								} else if (elems[i].startsWith("disease state:")
										|| elems[i].startsWith("disease status:")) {
									disease = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										diseasestate.set(i - 1, comps[1]);
										nrdiseaseset++;
									}

								} else if (elems[i].startsWith("adult age:")
										|| elems[i].startsWith("age:")
										|| elems[i].startsWith("subject age:")
										|| elems[i].startsWith("age(years):")
										|| elems[i].startsWith("age(yrs):")
										|| elems[i].startsWith("gestational age:")) {
									age = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										ages.set(i - 1, comps[1]);
										nrageset++;
									}

								} else if (elems[i].startsWith("child age:")) {
									age = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										childages.set(i - 1, comps[1]);
										nrageset++;
									}

								} else if (elems[i].startsWith("age (months):")) {
									age = true;
									String[] comps = elems[i].split(": ");
									if (comps.length > 0) {
										double aged = -1;
										try {
											aged = (Double.parseDouble(comps[1]) / 12);
										} catch (NumberFormatException e) {

										}
										childages.set(i - 1, "" + aged);
										nrageset++;
									}

								}
							}
							if (gender) {
//								System.out.println(ln);
								System.out.println(ANSI_BLUE + ln + ANSI_RESET);
								ArrayList<String> genderstr = new ArrayList<>();
								for (int i = 0; i < genders.size(); i++) {
									genderstr.add(genders.get(i).toString());
								}
								System.out.println(ANSI_BLUE + Strings.concat(genderstr, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_BLUE + "Found genders for " + nrgendersset + " samples" + ANSI_RESET);
							} else if (tissue) {

//								System.out.println(ln);
								System.out.println(ANSI_GREEN + ln + ANSI_RESET);
								System.out.println(ANSI_GREEN + Strings.concat(tissues, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_GREEN + "Found tissue for " + nrtissuesset + " samples" + ANSI_RESET);
							} else if (disease) {

								System.out.println(ANSI_CYAN + ln + ANSI_RESET);
								System.out.println(ANSI_CYAN + Strings.concat(diseasestate, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_CYAN + "Found disease for " + nrdiseaseset + " samples" + ANSI_RESET);
							} else if (age) {

//								System.out.println(ln);
								System.out.println(ANSI_PURPLE + ln + ANSI_RESET);
								System.out.println(ANSI_PURPLE + Strings.concat(ages, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_PURPLE + Strings.concat(childages, Strings.tab) + ANSI_RESET);
								System.out.println(ANSI_PURPLE + "Found age for " + nrageset + " samples" + ANSI_RESET);
							} else {
								System.out.println(ln);
							}

						} else if (ln.startsWith("!Sample_description")) {
//							System.out.println(ln);
						} else if (ln.startsWith("!Sample_geo_accession")) {
							for (int i = 1; i < elems.length; i++) {
								accids.add(elems[i]);
								genders.add(Gender.UNKNOWN);
								childgenders.add(Gender.UNKNOWN);
								tissues.add("Unknown");
								diseasestate.add("Unknown");
								sources.add("");
								ages.add("-1");
								celllines.add(SAMPLETYPE.UNKNOWN);
								taxids.add("Unknown");
								childages.add("-1");
							}
							System.out.println(ANSI_CYAN + accids.size() + " samples found.");
						} else if (ln.startsWith("!Sample_taxid_ch1")) {
							for (int i = 1; i < elems.length; i++) {
								taxids.set(i - 1, elems[i]);
							}
						} else if (ln.startsWith("!Sample_title")) {
//							System.out.println(ln);
						}
						//
					}
				} else {
					break;
				}
				ln = reader.readLine();

			}
			System.out.println();
			reader.close();
			con.disconnect();

			for (int i = 0; i < accids.size(); i++) {
				String lnout = gse
						+ "\t" + accids.get(i).toUpperCase()
						+ "\t" + genders.get(i)
						+ "\t" + ages.get(i)
						+ "\t" + childages.get(i)
						+ "\t" + childgenders.get(i)
						+ "\t" + celllines.get(i)
						+ "\t" + tissues.get(i)
						+ "\t" + sources.get(i)
						+ "\t" + taxids.get(i);
				tfout.writeln(lnout);
			}

		}

//		System.exit(-1);
	}


}
