package nl.harmjanwestra.playground.biogen.sqtl;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;

public class Recluster {
	private ArrayList<Feature> allFeatures;
	private HashMap<String, ArrayList<Feature>> featuresPerCluster;

	//	private static int[] minDarr = new int[]{5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
	private static int[] minDarr = new int[]{5};

	public static void main(String[] args) {

		String s = "d:\\cluster.txt";
		Recluster r = new Recluster();
		try {
			r.readClusters(s);
//			r.measureDistance();
			for (int d : minDarr) {
				r.recluster(d);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}


	public void readClusters(String clusterTextFile) throws IOException {
		TextFile tf = null;

		tf = new TextFile(clusterTextFile, TextFile.R);
		String ln = tf.readLine();

		allFeatures = new ArrayList<Feature>();
		featuresPerCluster = new HashMap<>();

		while (ln != null) {
			String[] elems = ln.split(":");
			Chromosome chr = Chromosome.parseChr(elems[0]);
			Integer sta = Integer.parseInt(elems[1]);
			Integer sto = Integer.parseInt(elems[2]);
			String cluster = elems[3];
			Feature f = new Feature(chr, sta, sto);
			f.setName(ln);


			allFeatures.add(f);

			ArrayList<Feature> fc = featuresPerCluster.get(cluster);
			if (fc == null) {
				fc = new ArrayList<>();
			}
			fc.add(f);
			featuresPerCluster.put(cluster, fc);

			ln = tf.readLine();
		}
		tf.close();
	}


	public class Cluster {
		Feature mergedFeature = null;
		ArrayList<Feature> includedFeatures = new ArrayList<>();

		public int getStart() {
			return mergedFeature.getStart();
		}

		public int getStop() {
			return mergedFeature.getStop();
		}

		public void merge(Cluster otherCluster, int startDistance, int stopDistance, int maxDistance) {

			Feature newMergedFeature = new Feature(mergedFeature);
			if (startDistance > 0 && startDistance < maxDistance) {
				if (this.getStart() > otherCluster.getStart()) {
					newMergedFeature.setStart(otherCluster.getStart());
				}
			}
			if (stopDistance > 0 && stopDistance < maxDistance) {
				if (this.getStop() < otherCluster.getStop()) {
					newMergedFeature.setStop(otherCluster.getStop());
				}
			}
			includedFeatures.addAll(otherCluster.includedFeatures);
			// System.out.println("Merging: " + mergedFeature.getName() + " with " + otherCluster.mergedFeature.getName() + " new position: " + newMergedFeature + " - startdistance: " + startDistance + " - stopdistance: " + stopDistance + " | junctions: " + includedFeatures.size());
			mergedFeature = newMergedFeature;
		}

	}

	public void recluster(int maxDistance) {


		for (int chr = 1; chr < 23; chr++) {
			ArrayList<Feature> features = new ArrayList<>();
			for (Feature f : allFeatures) {
				if (f.getChromosome().getNumber() == chr) {
					features.add(f);
				}
			}
			features.sort(new FeatureComparator());

			ArrayList<Cluster> clusters = new ArrayList<>();
			for (Feature f : features) {
				Cluster c = new Cluster();
				c.includedFeatures.add(f);
				c.mergedFeature = new Feature(f);
				clusters.add(c);
			}


			int originalNumberOfClusters = clusters.size();
			boolean thingsWereMerged = true;
			int iteration = 0;
			while (thingsWereMerged) {

				System.out.println("--- Iteration " + iteration + " - chr" + chr + " - current number of clusters: " + clusters.size());

				thingsWereMerged = false;
				ArrayList<Cluster> newClusters = new ArrayList<>();

				boolean[] alreadyMerged = new boolean[clusters.size()];
				for (int c1 = 0; c1 < clusters.size(); c1++) {
					if (!alreadyMerged[c1]) {
						Cluster currentCluster = clusters.get(c1);
						newClusters.add(currentCluster);

						for (int c2 = c1 + 1; c2 < clusters.size(); c2++) {
							if (!alreadyMerged[c2]) {
								Cluster otherCluster = clusters.get(c2);

								int startDistance = Math.abs(currentCluster.getStart() - otherCluster.getStart());
								int stopDistance = Math.abs(currentCluster.getStop() - otherCluster.getStop());

								// don't merge things with exactly the same start ór end position (they might be different actual splice events)
								if (startDistance < maxDistance && stopDistance < maxDistance) {
									// choose what to merge
									currentCluster.merge(otherCluster, startDistance, stopDistance, maxDistance);
									thingsWereMerged = true;
									alreadyMerged[c2] = true;
								}
							}
						}
					}
				}
				clusters = newClusters;
				iteration++;
				if (!thingsWereMerged) {
					double perc = ((double) clusters.size() / originalNumberOfClusters) * 100;
					System.out.println("Chr" + chr + " done - " + clusters.size() + " remain out of " + originalNumberOfClusters + " we started with (" + perc + ")");
					// System.exit(0);
				}
			}

		}
	}


	public void measureDistance() {


		// check whether there is overlap within clusters
		Set<String> clusterset = featuresPerCluster.keySet();
		ArrayList<String> clusters = new ArrayList<>();
		clusters.addAll(clusterset);
		Collections.sort(clusters);

		System.out.println("Within cluster");

		int[] minDarrCtrWithinClusters = new int[minDarr.length];
		int[] minDarrCtrBetweenClusters = new int[minDarr.length];
		for (int c = 0; c < clusters.size(); c++) {
			String cluster = clusters.get(c);
			ArrayList<Feature> features = featuresPerCluster.get(cluster);
			Collections.sort(features, new FeatureComparator());
			for (int f = 0; f < features.size(); f++) {
				Feature fObj1 = features.get(f);
				for (int f2 = f + 1; f2 < features.size(); f2++) {
					Feature fObj2 = features.get(f2);
					boolean overlap = fObj2.overlaps(fObj1);

					int d1 = Math.abs(fObj1.getStart() - fObj2.getStart());
					int d2 = Math.abs(fObj1.getStart() - fObj2.getStop());
					d2 = -1;
					int d3 = Math.abs(fObj1.getStop() - fObj2.getStop());
					int d4 = Math.abs(fObj1.getStop() - fObj2.getStart());
					d4 = -1;

//						System.out.println();

					for (int q = 0; q < minDarr.length; q++) {
						int minD = minDarr[q];
						if (d1 > 0 && d1 < minD
								|| d2 > 0 && d2 < minD
								|| d3 > 0 && d3 < minD
								|| d4 > 0 && d4 < minD
						) {
							minDarrCtrWithinClusters[q]++;
//								System.out.println(cluster + "\t" + f + "\t" + fObj1.getName() + "\t" + f2 + "\t" + fObj2.getName() + "\t" + overlap);
//								System.out.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4);
							// System.exit(0);
						}
					}


				}

			}
		}


		System.out.println("Between clusters");

		for (int c = 1; c < 23; c++) {
			ArrayList<Feature> chrFeatures = new ArrayList<>();
			for (Feature f : allFeatures) {
				if (f.getChromosome().getNumber() == c) {
					chrFeatures.add(f);
				}
			}
			Collections.sort(chrFeatures, new FeatureComparator());
			for (int f = 0; f < chrFeatures.size(); f++) {
				Feature fObj1 = chrFeatures.get(f);
				for (int f2 = f + 1; f2 < chrFeatures.size(); f2++) {
					Feature fObj2 = chrFeatures.get(f2);
					boolean overlap = fObj2.overlaps(fObj1);

					int d1 = Math.abs(fObj1.getStart() - fObj2.getStart());
					int d2 = Math.abs(fObj1.getStart() - fObj2.getStop());
					d2 = -1;
					int d3 = Math.abs(fObj1.getStop() - fObj2.getStop());
					int d4 = Math.abs(fObj1.getStop() - fObj2.getStart());
					d4 = -1;


//						System.out.println();

					for (int q = 0; q < minDarr.length; q++) {
						int minD = minDarr[q];
						if (d1 > 0 && d1 < minD
								|| d2 > 0 && d2 < minD
								|| d3 > 0 && d3 < minD
								|| d4 > 0 && d4 < minD
						) {
							minDarrCtrBetweenClusters[q]++;
//								System.out.println("chr" + c + "\t" + f + "\t" + fObj1.getName() + "\t" + f2 + "\t" + fObj2.getName() + "\t" + overlap);
//								System.out.println(d1 + "\t" + d2 + "\t" + d3 + "\t" + d4);
							// System.exit(0);
						}
					}


				}

			}
		}

		System.out.println();
		System.out.println();
		System.out.println("------");

		for (int q = 0; q < minDarr.length; q++) {
			int minD = minDarr[q];
			System.out.println(minD + "\t" + minDarrCtrWithinClusters[q] + "\t" + minDarrCtrBetweenClusters[q]);
		}


	}
}
