/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.zscoretables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author harmjan
 */
public class ZScoreTableManipulation {

    public String filterZScoreTableForSetOfProbes(String probenrFile, String matrixIn) throws IOException {
	System.out.println("Filtering probes for " + probenrFile);
	// load probe translation
	TextFile tf = new TextFile(probenrFile, TextFile.R);

	String line = tf.readLine();
	HashSet<String> includeProbes = new HashSet<String>();
	while (line != null) {

	    includeProbes.add(line.trim());
	    line = tf.readLine();
	}

	tf.close();

	TextFile tf2 = new TextFile(matrixIn, TextFile.R);
	String[] headerelems = tf2.readLineElems(TextFile.tab);
	boolean[] includecols = new boolean[headerelems.length];
	includecols[0] = false;
	int ctr = 0;

	for (int i = 0; i < headerelems.length; i++) {
	    if (includeProbes.contains(headerelems[i])) {
		includecols[i] = true;
		ctr++;
	    }
	}

	System.out.println(ctr + " cols selected");
//	TextFile out = new TextFile("/Volumes/BackupDisk/Meta0/trans/2012-01-18-Groningen-CISTRANS-40PCs-4GWASPCs-GeneticVectorsNotRemoved/metazscoretable_HT12v3.txt", TextFile.W);
//	out.writeln(header);
	int lnctr = 0;
	String[] elems = tf2.readLineElems(TextFile.tab);
	int[] colNullCtr = new int[headerelems.length];
	while (elems != null) {
	    for (int i = 1; i < elems.length; i++) {
		if (includecols[i]) {
		    if (elems[i].equals("null")) {
			includecols[i] = false;
			colNullCtr[i]++;
		    }

		}
	    }

	    elems = tf2.readLineElems(TextFile.tab);
//	    System.out.println(lnctr);
	    lnctr++;
	}

//	for(int i=0; i<colNullCtr.length; i++){
//	    if(colNullCtr[i]>0){
//		System.out.println(headerelems[i]+"\t"+colNullCtr[i]);
//	    }
//	}

	String header = "snp";
	for (int i = 0; i < headerelems.length; i++) {
	    if (includecols[i]) {
		header += "\t" + headerelems[i];
	    }
	}

	tf2.open();
	elems = tf2.readLineElems(TextFile.tab); // header

	elems = tf2.readLineElems(TextFile.tab);

	System.out.println("Writing to file");
	String filenameout = matrixIn + "-filtered.txt";
	TextFile out = new TextFile(filenameout, TextFile.W);

	out.writeln(header);
	while (elems != null) {

	    StringBuilder lineout = new StringBuilder();
	    lineout.append(elems[0]);
	    for (int i = 1; i < elems.length; i++) {
		if (includecols[i]) {
		    lineout.append("\t").append(elems[i]);
		}
	    }

	    out.writeln(lineout.toString());

	    elems = tf2.readLineElems(TextFile.tab);
//	    System.out.println(lnctr);
	    lnctr++;
	}

	out.close();
	tf2.close();


	return filenameout;

    }

    public String convertProbeNumbersToEnsemblId(String ensemblannot, String file) throws IOException {
	System.out.println("Converting probe names to ensembl ids");
	TextFile etf = new TextFile(ensemblannot, TextFile.R);

	String[] elems = etf.readLineElems(TextFile.tab);
	HashMap<String, String> probeToEns = new HashMap<String, String>();
	while (elems != null) {

	    if (elems.length >= 5) {
		String probe = elems[0].trim();
		String ens = elems[4].trim();

		probeToEns.put(probe, ens);
	    }
	    elems = etf.readLineElems(TextFile.tab);
	}
	etf.close();

	TextFile tf2 = new TextFile(file, TextFile.R);


	String filenameout = file + "-ens.txt";
	TextFile tf2out = new TextFile(filenameout, TextFile.W);

	elems = tf2.readLineElems(TextFile.tab);

	String header = "snp";
	boolean[] includecol = new boolean[elems.length];
	for (int i = 1; i < elems.length; i++) {
	    String ens = probeToEns.get(elems[i]);
	    if (ens == null) {
		includecol[i] = false;
	    } else {
		header += "\t" + ens;
		includecol[i] = true;
	    }

	}

	tf2out.writeln(header);
	elems = tf2.readLineElems(TextFile.tab);
	int lnctr = 0;
	while (elems != null) {

	    StringBuilder b = new StringBuilder();
	    b.append(elems[0]);

	    for (int i = 1; i < elems.length; i++) {
		if (includecol[i]) {
		    b.append("\t").append(elems[i]);
		}
	    }

	    tf2out.writeln(b.toString());
	    elems = tf2.readLineElems(TextFile.tab);
//	    System.out.println(lnctr);
	    lnctr++;
	}

	tf2out.close();
	tf2.close();
	return filenameout;
    }

    public String collapseDuplicateProbes(String file) throws IOException {
	System.out.println("Collapsing duplicate genes");
	DoubleMatrixDataset<String, String> ds = new DoubleMatrixDataset<String, String>(file);

	System.out.println("Loaded matrix is: " + ds.nrRows + " x " + ds.nrCols);
	List<String> cols = ds.colObjects;
	HashMap<String, Integer> uniquecols = new HashMap<String, Integer>();
	HashMap<Integer, String> uniquecolsrev = new HashMap<Integer, String>();


	int itr = 0;
	System.out.println("Iterating through " + cols.size() + " columns");
	for (int i = 0; i < cols.size(); i++) {
	    if (!uniquecols.containsKey(cols.get(i))) {
		uniquecols.put(cols.get(i), itr);
		uniquecolsrev.put(itr, cols.get(i));

		itr++;
	    }
	}

	System.out.println(uniquecols.size() + " unique genes");
	double[][] sums = new double[ds.nrRows][uniquecols.size()];
	int[][] occurrences = new int[ds.nrRows][uniquecols.size()];

	System.out.println("Creating matrix: " + sums.length + " x " + sums[0].length);

	for (int r = 0; r < ds.nrRows; r++) {

	    for (int c = 0; c < ds.nrCols; c++) {

		Integer sumCol = uniquecols.get(cols.get(c));
		if (sumCol == null) {
		    System.out.println(cols.get(c));
		} else {
		    sums[r][sumCol] += ds.rawData[r][c];
		    occurrences[r][sumCol]++;
		}

	    }
	}

	for (int i = 0; i < sums.length; i++) {
	    for (int j = 0; j < sums[i].length; j++) {
		sums[i][j] /= occurrences[i][j];
	    }
	}


	ArrayList<String> colsnames = new ArrayList<String>();
	for (int i = 0; i < uniquecols.size(); i++) {
	    colsnames.add(uniquecolsrev.get(i));
	}




	DoubleMatrixDataset<String, String> ds2 = new DoubleMatrixDataset<String, String>(sums);
	ds2.rowObjects = ds.rowObjects;
	ds2.colObjects = colsnames;
	String filenameout = file + "-collapsed.txt";
	ds2.save(filenameout);


	return filenameout;

    }

    private void mergeTwoZScoreTables(String string, String string0, String outfilename) throws IOException {
	DoubleMatrixDataset<String, String> d1 = new DoubleMatrixDataset<String, String>(string);
	DoubleMatrixDataset<String, String> d2 = new DoubleMatrixDataset<String, String>(string0);


	List<String> cols1 = d1.colObjects;
	List<String> cols2 = d2.colObjects;

	HashMap<String, Integer> toNewId = new HashMap<String, Integer>();
	int itr = 0;
	ArrayList<String> cols = new ArrayList<String>();
	for (int i = 0; i < cols1.size(); i++) {

	    if (d2.hashCols.containsKey(cols1.get(i))) {
		cols.add(cols1.get(i));
		toNewId.put(cols1.get(i), itr);
		itr++;
	    }
	}

	System.out.println("Overlap: " + toNewId.size() + " elements.");


	double[][] newRawData = new double[d1.nrRows + d2.nrRows][toNewId.size()];
	ArrayList<String> rows = new ArrayList<String>();
	for (int r = 0; r < d1.nrRows; r++) {
	    rows.add(d1.rowObjects.get(r));
	    for (int c = 0; c < d1.nrCols; c++) {
		Integer colNr = toNewId.get(cols1.get(c));
		if (colNr != null) {
		    newRawData[r][colNr] = d1.rawData[r][c];
		}
	    }
	}

	for (int r = 0; r < d2.nrRows; r++) {
	    rows.add(d2.rowObjects.get(r));
	    for (int c = 0; c < d2.nrCols; c++) {
		Integer colNr = toNewId.get(cols2.get(c));
		if (colNr != null) {
		    newRawData[r + d1.nrRows][colNr] = d2.rawData[r][c];
		}
	    }
	}



	TextFile outfile = new TextFile(outfilename, TextFile.W);
	String header = "snp";
	for (int i = 0; i < cols.size(); i++) {
	    header += "," + cols.get(i);
	}
	outfile.writeln(header);
	for (int r = 0; r < newRawData.length; r++) {
	    StringBuilder out = new StringBuilder();

	    out.append(rows.get(r));
	    for (int c = 0; c < newRawData[r].length; c++) {
		out.append(",").append(newRawData[r][c]);
	    }

	    outfile.writeln(out.toString());
	}

	outfile.close();






    }

    public void convertToBinary(String string) throws IOException {
	DoubleMatrixDataset<String, String> ds1 = new DoubleMatrixDataset<String, String>(string);
	ds1.save(string + ".binary");
    }
}
