/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.eqtlfile;

import eqtlmappingpipeline.meta.posthocanalysis.network.eQTLEdge;
import eqtlmappingpipeline.meta.posthocanalysis.network.eQTLNetwork;
import eqtlmappingpipeline.meta.posthocanalysis.network.eQTLVertex;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.gwascatalog.GWASCatalog;
import umcg.genetica.io.gwascatalog.GWASSNP;
import umcg.genetica.io.gwascatalog.GWASTrait;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class eQTLSummarize {

    GWASCatalog catalog;

    public void eQTLSummarize(String gwasCatalog, String transeQTLFile, String ciseQTLFile, String geneCorrelationMatrix, String out, String query, String ldfile, Double ldthreshold) throws IOException {

	if (!out.endsWith("/")) {
	    out += "/";
	}
	if (!Gpio.exists(out)) {
	    Gpio.createDir(out);
	}

	TextFile ld = new TextFile(ldfile, TextFile.R);

	ld.readLine();
	String[] ldelems = ld.readLineElems(TextFile.tab);

	if (ldthreshold == null) {
	    ldthreshold = 0d;
	}
	System.out.println("using ld threshold: " + ldthreshold);
	HashMap<String, String> snpmergetable = new HashMap<String, String>();
	HashMap<String, ArrayList<String>> snpmergetableRev = new HashMap<String, ArrayList<String>>();
	while (ldelems != null) {

	    String source = ldelems[0];
	    String target = ldelems[1];
	    Double ldval = Double.parseDouble(ldelems[2]);




	    if (ldval >= ldthreshold) {

		if (!snpmergetable.containsKey(source)) {
		    snpmergetable.put(target, source);
		}

	    }


	    ldelems = ld.readLineElems(TextFile.tab);
	}


//	TextFile outld = new TextFile(out+"ld-"+ldthreshold, TextFile.W);
//
//	Set<String> keys = snpmergetable.keySet();
//	for(String s: keys){
//	    outld.writeln(s+"\t"+snpmergetable.get(s));
//	}
//
//	outld.close();

	catalog = new GWASCatalog();

	catalog.read(gwasCatalog);

	HashSet<GWASSNP> gsnps = catalog.getSnps();
	HashSet<String> allowedSNPs = new HashSet<String>();

	GWASSNP[] allowedSNParr = new GWASSNP[allowedSNPs.size()];
	allowedSNParr = gsnps.toArray(allowedSNParr);
	for (GWASSNP s : allowedSNParr) {
	    allowedSNPs.add(s.getName());


	}

	eQTLNetwork net = new eQTLNetwork();
	if (transeQTLFile != null) {
	    System.out.println("Adding trans eqtls from " + transeQTLFile);
	    net.addeQTLsFromFile(transeQTLFile, false);
	}
	if (ciseQTLFile != null) {
	    System.out.println("Adding cis eqtls from " + ciseQTLFile);
	    net.addeQTLsFromFile(ciseQTLFile, true);
	}

	System.out.println("Done reading net");
	GWASTrait[] traits = catalog.getTraits();
	System.out.println("num traits: " + traits.length);
	ArrayList<GWASTrait> traitsWithoutEQTLs = new ArrayList<GWASTrait>();




//	findConvergence(catalog, net, snpmergetable, out);
//	outputnovel(catalog, net, out, query, snpmergetable);

//	findConvergenceForDiseases(catalog, net, query);
	findConvergenceForDiseases(catalog, net, query, out, snpmergetable);







    }

    public void findConvergenceForDiseases(GWASCatalog catalog, eQTLNetwork net, String query, String out, HashMap<String, String> mergetable) throws IOException {


	GWASTrait[] t = catalog.getTraits();
	ArrayList<GWASTrait> traits = new ArrayList<GWASTrait>();

	for (GWASTrait t2 : t) {
	    if (t2.getName().contains(query)) {
		traits.add(t2);
	    }
	}

	TextFile outfile = new TextFile(out+"ConvergencePerTrait.txt", TextFile.W);
	TextFile outfiledisease = new TextFile(out+"ConvergenceFor-"+query+".txt", TextFile.W);

	outfile.writeln("trait\tsnp1\tsnpchr1\tsnpchrpos1\tsnp1EType\tsnp2\tsnpchr2\tsnpchrpos2\tsnp2EType\tSNP2IsAssocWithDisease\tSNP1And2AreInLD\tSNP2IsInLDWithDiseaseAssocSNP\tGene\tGeneChr\tGeneChrPos");
	//Platelet aggregation	rs2893923	10	64931190	rs10927101	1	242240495	false	false	false	AQP10
	for (GWASTrait t2 : t) {



	    GWASSNP[] snps = t2.getSNPs();

	    HashSet<String> snpsForTrait = new HashSet<String>();
	    for (GWASSNP s : snps) {
		snpsForTrait.add(s.getName());
	    }

	    for (GWASSNP s : snps) {
		eQTLVertex snp1v = net.getVertex(s.getName());
		if (snp1v != null) {
		    eQTLEdge[] edges = snp1v.getEdges();

		    for (eQTLEdge snp1E : edges) {
			eQTLVertex probev = snp1E.getA();
			if (probev == snp1v) {
			    probev = snp1E.getB();
			}

			eQTLEdge[] edges2 = probev.getEdges();
			for (eQTLEdge probeE : edges2) {
//			    if (probeE.getType().equals(eQTLEdge.TYPE.TRANS)) {
			    eQTLVertex snp2v = probeE.getA();
			    if (snp2v == probev) {
				snp2v = probeE.getB();
			    }

			    boolean snp2IsDiseaseRelated = false;
			    if (snpsForTrait.contains(snp2v.getName())) {
				snp2IsDiseaseRelated = true;
			    }

			    String snpalter1 = mergetable.get(snp1v.getName());
			    String snpalter2 = mergetable.get(snp2v.getName());

			    if (snpalter1 == null) {
				snpalter1 = snp1v.getName();
			    }

			    if (snpalter2 == null) {
				snpalter2 = snp2v.getName();
			    }


			    boolean snp2IsInLDWithDiseaseRelatedSNP = false;
			    if (snpsForTrait.contains(snpalter2)) {
				snp2IsInLDWithDiseaseRelatedSNP = true;

			    }

			    boolean snpsAreInLD = false;
			    if (snpalter1.equals(snpalter2)) {
				snpsAreInLD = true;
			    }

			    boolean output = false;
			    if (snp1v.getChr() != snp2v.getChr()) {
				// convergence

				output = true;
			    } else {
				if (Math.abs(snp1v.getChrPos() - snp2v.getChrPos()) > 5000000) {
				    output = true;
				}
			    }

			    String e1type = "trans";
			    if (snp1E.getType() == eQTLEdge.TYPE.CIS) {
				e1type = "cis";
			    }

			    String e2type = "trans";
			    if (probeE.getType() == eQTLEdge.TYPE.CIS) {
				e2type = "cis";
			    }

			    if (output) {
				String ln = t2.getName() + "\t" + snp1v.getName() + "\t" + snp1v.getChr() + "\t" + snp1v.getChrPos() + "\t" + e1type + "\t"
					+ snp2v.getName() + "\t" + snp2v.getChr() + "\t" + snp2v.getChrPos() +"\t" + e2type
					+ "\t" + snp2IsDiseaseRelated + "\t" + snpsAreInLD + "\t" + snp2IsInLDWithDiseaseRelatedSNP + "\t"
					+ probev.getName() + "\t" + probev.getChr() + "\t" + probev.getChrPos();
				System.out.println(ln);

				outfile.writeln(ln);

				if(t2.getName().toLowerCase().contains(query)){
				    outfiledisease.writeln(snp1v.getName()+"\t"+e1type+"\t"+ probev.getName()+"\t"+snp1E.getWeight());
				    outfiledisease.writeln(snp2v.getName()+"\t"+e2type+"\t"+ probev.getName()+"\t"+probeE.getWeight());
				}
			    }

//			    }

			}

		    }
		}

	    }

	}
	outfiledisease.close();
	outfile.close();


    }

    private void findConvergence(GWASCatalog catalog, eQTLNetwork net, HashMap<String, String> snpMergeTable, String out) throws IOException {

	eQTLVertex[] nodes = net.getVertexes();


	TextFile outfile = new TextFile(out, TextFile.W);
	HashSet<String> allConvergentSNPs = new HashSet<String>();

	for (eQTLVertex v : nodes) {
	    if (v.getType().equals(eQTLVertex.TYPE.PROBE)) {
		eQTLEdge[] edges = v.getEdges();

		HashSet<String> uniqueSNPs = new HashSet<String>();
		for (eQTLEdge e : edges) {

		    if (e.getType().equals(eQTLEdge.TYPE.TRANS)) {
			eQTLVertex v2 = e.getA();
			if (!v2.getType().equals(eQTLVertex.TYPE.SNP)) {
			    v2 = e.getB();
			}

			String snpname = v2.getName();
			String altersnp = snpMergeTable.get(snpname);
			if (altersnp != null) {
			    snpname = altersnp;
			}

			uniqueSNPs.add(snpname);
		    }

		}

		if (uniqueSNPs.size() > 1) {

//		    // convergence found!
//		    String[] snpsarr = new String[uniqueSNPs.size()];
//		    snpsarr = uniqueSNPs.toArray(snpsarr);
//		    String output = v.getName() + "\t" + uniqueSNPs.size();
//		    for (String s : snpsarr) {
//			output += "\t" + s;
//
//		    }
//		    System.out.println(output);


		    for (eQTLEdge e : edges) {

//			if (e.getType().equals(eQTLEdge.TYPE.TRANS)) {
			eQTLVertex v2 = e.getA();
			if (!v2.getType().equals(eQTLVertex.TYPE.SNP)) {
			    v2 = e.getB();
			}

			String snpname = v2.getName();
			String altersnp = snpMergeTable.get(snpname);
			if (altersnp != null) {
			    snpname = altersnp;
			}

//			    uniqueSNPs.add(snpname);
			String etype = null;
			if (e.getType() == eQTLEdge.TYPE.CIS) {
			    etype = "cis";
			} else {
			    etype = "trans";
			}
//	    eQTLVertex tmpv = null;
//	    if (v.getType() == eQTLVertex.TYPE.PROBE) {
//		tmpv = v;
//		v = v2;
//		v2 = tmpv;
//	    }

			String ln = snpname + "\t" + etype + "\t" + v.getName() + "\t" + e.getWeight();
			System.out.println(ln);
			outfile.writeln(ln);

			eQTLEdge[] v2edges = v2.getEdges();
			for (eQTLEdge e2 : v2edges) {

			    eQTLVertex v3 = e2.getA();
			    if (!v3.getType().equals(eQTLVertex.TYPE.PROBE)) {
				v3 = e2.getB();
			    }

			    if (e.getType() == eQTLEdge.TYPE.CIS) {
				etype = "cis";
			    } else {
				etype = "trans";
			    }
			    ln = snpname + "\t" + etype + "\t" + v3.getName() + "\t" + e.getWeight();

			    System.out.println(ln);
			    outfile.writeln(ln);

			}

//			}

		    }


		}
	    }

	}

	outfile.close();

	filteroutDuplicateEdges(out);

    }

    private void outputnovel(GWASCatalog catalog, eQTLNetwork net, String outfilename, String query, HashMap<String, String> snpmergetable) throws IOException {
	// header:
	// source target type score

	TextFile outfile = new TextFile(outfilename, TextFile.W);
	outfile.writeln("source\ttype\ttarget\tweight");
	// start with trait relationships for snps
	GWASTrait[] traits = catalog.getTraits();
	HashSet<GWASTrait> allowedTraits = new HashSet<GWASTrait>();

	int snpctr = 0;
	int snpnotpresent = 0;
	for (GWASTrait t : traits) {
	    if (t.getName().toLowerCase().contains(query)) {
		GWASSNP[] snps = t.getSNPs();
		allowedTraits.add(t);
		for (GWASSNP s : snps) {
		    if (net.getVertex(s.getName()) != null) {
			String ln = s.getName() + "\t" + t.getName() + "\tgwas\t" + 0d;
			System.out.println(ln);
//                        outfile.writeln(ln);
			snpctr++;
		    } else {
			System.out.println(s.getName() + "\t not found");
			snpnotpresent++;
		    }
		}
	    }

	}
	int total = snpctr + snpnotpresent;
	System.out.println(total + "\t" + snpctr + "\t" + snpnotpresent);
	// now the snp -> cis effects
	HashSet<GWASSNP> gsnps = catalog.getSnps();
	GWASSNP[] allowedSNParr = new GWASSNP[gsnps.size()];
	allowedSNParr = gsnps.toArray(allowedSNParr);
	HashSet<eQTLEdge> visitedEdges = new HashSet<eQTLEdge>();
	for (GWASSNP snp : allowedSNParr) {
	    GWASTrait[] assocTraits = snp.getAssociatedTraitsArray();
	    boolean includesnp = false;
	    for (GWASTrait t : assocTraits) {
		if (allowedTraits.contains(t)) {
		    includesnp = true;
		}
	    }
	    if (includesnp) {

		String alterSNP = snpmergetable.get(snp.getName());
		if (alterSNP == null) {
		    alterSNP = snp.getName();
		}

		eQTLVertex v = net.getVertex(snp.getName());
		if (v != null) {
		    iterateVertex(v, visitedEdges, outfile, alterSNP);

		}
	    }


	}
	outfile.close();

	filteroutDuplicateEdges(outfilename);

    }

    private void outputlegacy(GWASTrait[] traits, eQTLNetwork net, String outfilename) throws IOException {
	TextFile outfile = new TextFile(outfilename, TextFile.W);
	for (GWASTrait trait : traits) {

	    String traitname = trait.getName();
	    GWASSNP[] snps = trait.getSNPs();
//		System.out.println(trait.getName() + "\t" + snps.length);
	    for (GWASSNP snp : snps) {
//		    System.out.println(snp.getName());
		String snpname = snp.getName();

		eQTLVertex v = net.getVertex(snpname);
		if (v != null) {
		    eQTLEdge[] e = v.getEdges();

		    if (e.length > 0) {
			String eQTLOutput = snpname + "\t" + v.getChr() + "\t" + v.getChrPos() + "\t" + e.length;
			for (int i = 0; i < e.length; i++) {
			    String ename = e[i].getName();
			    eQTLVertex v2 = e[i].getA();

			    if (v == v2) {
				v2 = e[i].getB();
			    }

			    if (v2.getType() == eQTLVertex.TYPE.PROBE) {
				eQTLOutput += "\t" + v2.getAnnotation();
			    } else {
				System.out.println("ERROR: eQTL effect with SNP?");
			    }


			}
			System.out.println(traitname + "\t" + eQTLOutput);

		    }
		}


	    }
	}
	outfile.close();
    }

    private void iterateVertex(eQTLVertex v, HashSet<eQTLEdge> visitedEdges, TextFile outfile, String alternateSNPName) throws IOException {
	eQTLEdge[] edges = v.getEdges();
	visitedEdges = new HashSet<eQTLEdge>();

	for (eQTLEdge e : edges) {

	    visitedEdges.add(e);
	    eQTLVertex v2 = e.getA();
	    if (v == v2) {
		v2 = e.getB();
	    }
	    String etype = null;
	    if (e.getType() == eQTLEdge.TYPE.CIS) {
		etype = "cis";
	    } else {
		etype = "trans";
	    }
//	    eQTLVertex tmpv = null;
//	    if (v.getType() == eQTLVertex.TYPE.PROBE) {
//		tmpv = v;
//		v = v2;
//		v2 = tmpv;
//	    }
	    String ln = alternateSNPName + "\t" + etype + "\t" + v2.getName() + "\t" + e.getWeight();
	    System.out.println(ln);
	    outfile.writeln(ln);
//                iterateVertex(v2, visitedEdges, outfile);



	}
    }

    private void filteroutDuplicateEdges(String outfilename) throws IOException {
	// filter duplicate edges
	TextFile tf = new TextFile(outfilename, TextFile.R);
	TextFile tf2out = new TextFile(outfilename + "-filteredForDuplicateEdges.txt", TextFile.W);

	HashSet<String> visitedEdgesStr = new HashSet<String>();

	String line = tf.readLine();
	while (line != null) {

	    String[] elems = line.split("\t");

	    if (!visitedEdgesStr.contains(elems[0] + "-" + elems[2])) {

		tf2out.writeln(line);

		visitedEdgesStr.add(elems[0] + "-" + elems[2]);
	    }

	    line = tf.readLine();
	}

	tf2out.close();

	tf.close();
    }
}
