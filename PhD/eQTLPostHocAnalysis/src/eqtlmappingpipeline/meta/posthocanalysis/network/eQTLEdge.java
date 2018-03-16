/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.network;

/**
 *
 * @author harmjan
 */
public class eQTLEdge {

    public String getName() {
	return name;
    }

    public eQTLVertex getA(){
	return a;
    }
    
    public eQTLVertex getB(){
	return b;
    }
    
    public TYPE getType(){
	return t;
    }
    
    
    public static enum TYPE {

	CIS, TRANS
    };

    public static enum DIRECTIONALITY {

	DIRECTED, UNDIRECTED
    };
    private TYPE t;
    private double weight;
    private String[] datasets;
    private double[] datasetzscores;
    private int hashCode = 0;
    protected DIRECTIONALITY dir;
    protected eQTLVertex a;
    protected eQTLVertex b;
    String name;

    /**
     * @return the dir
     */
    public DIRECTIONALITY getDir() {
	return dir;
    }

    /**
     * @param dir the dir to set
     */
    public void setDir(DIRECTIONALITY dir) {
	this.dir = dir;
    }

    public eQTLEdge(eQTLVertex a, eQTLVertex b, TYPE t, double w) {
	this.t = t;
	this.a = a;
	this.b = b;
	this.weight = w;
	dir = DIRECTIONALITY.DIRECTED;
    }

    @Override
    public boolean equals(Object obj) {
	if (this == obj) {
	    return true;
	}
	if (!(obj instanceof eQTLEdge)) {
	    return false;
	}
	eQTLEdge otheredge = (eQTLEdge) obj;
	return this.hashCode() == otheredge.hashCode();

    }

    @Override
    public int hashCode() {
	final int multiplier = 23;
	if (hashCode == 0) {
	    int code = 133;


	    int anamehash = 0;
	    int bnamehash = 0;
	    if (a != null) {
		anamehash = a.hashCode();
	    }

	    if (b.getName() != null) {
		bnamehash = b.hashCode();
	    }

	    if (dir == DIRECTIONALITY.UNDIRECTED) {
		code = multiplier * code + (anamehash + bnamehash);

	    } else {
		code = multiplier * code + anamehash;
		code = multiplier * code + bnamehash;
	    }


	    code = multiplier * code + getDir().hashCode();
	    code = multiplier * code + t.hashCode();
	    hashCode = code;
	}
	return hashCode;
    }
    
    public double getWeight(){
        return weight;
    }
}
