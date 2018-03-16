/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.network;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author harmjan
 */
public class eQTLNetwork {

//    ArrayList<eQTLEdge> edges = new ArrayList<eQTLEdge>();
//    ArrayList<eQTLVertex> vertices = new ArrayList<eQTLVertex>();
//    HashMap<String, eQTLVertex> vertexNameToVertexObj = new HashMap<String, eQTLVertex>();
    ArrayList<eQTLEdge> edges = new ArrayList<eQTLEdge>();
    ArrayList<eQTLVertex> vertices = new ArrayList<eQTLVertex>();
    HashMap<String, eQTLVertex> vertexNameToVertexObj = new HashMap<String, eQTLVertex>();

    public void addVertex(eQTLVertex v) {
	vertices.add(v);
	vertexNameToVertexObj.put(v.getName(), v);
    }

    public void addEdge(eQTLEdge e) {
	edges.add(e);
    }

    public eQTLVertex getVertex(String name) {
	return vertexNameToVertexObj.get(name);
    }

    public void addeQTLsFromFile(String eqtlfile, boolean cis) throws IOException {
	TextFile tf = new TextFile(eqtlfile, TextFile.R);

	int nreqtls = tf.countLines();
	System.out.println(nreqtls + " eqtls detected");
	tf.open();
	tf.readLine();
	String[] elems = tf.readLineElemsReturnReference(TextFile.tab);
	ProgressBar pb = new ProgressBar(nreqtls);
	int lnctr = 1;
	while (elems != null) {
	    String snp = elems[1];

	    byte snpchr = ChrAnnotation.parseChr(elems[2]);
	    int snpchrpos = -1;
	    try {
		snpchrpos = Integer.parseInt(elems[3]);
	    } catch (NumberFormatException e) {
	    }

	    String finalsnpname = snp;

	    eQTLVertex snpv = new eQTLVertex(finalsnpname, null, eQTLVertex.TYPE.SNP, snpchr, snpchrpos);
	    if (vertices.contains(snpv)) {
		Integer snpi = vertices.indexOf(snpv);
		snpv = vertices.get(snpi);
	    }
	    String probe = elems[4];
	    byte probechr = ChrAnnotation.parseChr(elems[5]);
	    int probechrpos = -1;
	    try {
		probechrpos = Integer.parseInt(elems[6]);
	    } catch (NumberFormatException e) {
	    }
	    String annot = elems[16];
	    probe = elems[16];
	    if (probe.equals("-")) {
		probe = elems[4];
	    }
	    eQTLVertex probev = new eQTLVertex(probe, annot, eQTLVertex.TYPE.PROBE, probechr, probechrpos);

	    if (vertices.contains(probev)) {
		Integer probei = vertices.indexOf(probev);
		probev = vertices.get(probei);
	    }

	    addVertex(snpv);
	    addVertex(probev);
	    double zscore = 0;
	    try {

		zscore = Double.parseDouble(elems[10]);

		zscore = Math.abs(zscore);
	    } catch (NumberFormatException e) {
	    }
	    eQTLEdge e = null;
	    if (probechr != -1 && snpchr != -1) {
		if (probechr == snpchr) {

		    if (Math.abs(probechrpos - snpchrpos) < 5000000) {
			e = new eQTLEdge(snpv, probev, eQTLEdge.TYPE.CIS, zscore);
		    } else {
			e = new eQTLEdge(snpv, probev, eQTLEdge.TYPE.TRANS, zscore);
		    }
		} else {
		    e = new eQTLEdge(snpv, probev, eQTLEdge.TYPE.TRANS, zscore);
		}

		if (probev.getName().equals("TAGAP")) {
		    System.out.println("TAGAP: " + snpv.getName());

		}
		snpv.addEdge(e);
		probev.addEdge(e);
		addEdge(e);
	    } else {
		// skip
//		    System.out.println("eQTL Not added: probe\t"+probe+" chr/pos: "+probechr+"-"+probechrpos+"\t| snp\t"+snp+" chr/pos "+snpchr+"-"+snpchrpos);
	    }




	    elems = tf.readLineElemsReturnReference(TextFile.tab);
	    pb.iterate();

	}

//	System.exit(0);
	pb.close();
	System.out.println(edges.size() + " edges in network");
	System.out.println(vertices.size() + " vertices in network");
	tf.close();
    }

    public eQTLVertex[] getVertexes() {
	eQTLVertex[] vertexarr = new eQTLVertex[vertices.size()];

	vertexarr = vertices.toArray(vertexarr);
	return vertexarr;
    }
}
