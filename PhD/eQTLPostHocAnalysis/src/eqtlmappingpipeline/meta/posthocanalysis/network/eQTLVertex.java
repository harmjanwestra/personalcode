/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.network;

import java.util.HashSet;

/**
 *
 * @author harmjan
 */
public class eQTLVertex {

    public eQTLEdge[] getEdges() {
	eQTLEdge[] traitarray = new eQTLEdge[edges.size()];
	edges.toArray(traitarray);
	return traitarray;
    }

    public TYPE getType(){
	return t;
    }

    public String getAnnotation() {
	return annotation;
    }

    public String getName() {
	return name;
    }

    public byte getChr() {
	return chr;
    }

    public int getChrPos() {
	return chrpos;
    }
    

    public static enum TYPE {PROBE, SNP };
    
    private TYPE t;
    private HashSet<eQTLEdge> edges = new HashSet<eQTLEdge>();
    private String name, annotation;
    private byte chr;
    private int chrpos;
    private int hashCode = 0;
    
    public eQTLVertex(String name, String annotation, TYPE t, byte chr, int chrpos){
	this.t= t;
	this.name = name;
	this.annotation = annotation;
	this.chr = chr;
	this.chrpos = chrpos;
    }
    
    void addEdge(eQTLEdge e) {
	edges.add(e);
    }
    
    
    @Override
    public boolean equals(Object obj) {
	if (this == obj) {
	    return true;
	}
	if (!(obj instanceof eQTLVertex)) {
	    return false;
	}
	eQTLVertex othervertex = (eQTLVertex) obj;
	if (this.hashCode() != othervertex.hashCode()) {
	    return false;
	} else {
	    return true;
	}

    }

    
    @Override
    public int hashCode() {
	final int multiplier = 23;
	if (hashCode == 0) {
	    int code = 132;
	    code = multiplier * code + name.hashCode();
	    code = multiplier * code + chr;
	    code = multiplier * code + chrpos;
	    hashCode = code;
	}
	return hashCode;
    }
}
