/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis;

/**
 *
 * @author harmjan
 */
public class SortableSNP implements Comparable<SortableSNP> {

    private String name;
    private byte chr;
    private int chrpos;

    public SortableSNP(String name, byte chr, int chrpos) {
	this.name = name;
	this.chr = chr;
	this.chrpos = chrpos;
    }

    public byte getchr() {
	return chr;
    }

    public int getchrpos() {
	return chrpos;
    }

    public String getname() {
	return name;
    }

    @Override
    public int compareTo(SortableSNP t) {
	if (t.equals(this)) {
	    return 0;
	} else if (chr == t.getchr()) {
	    if (chrpos > t.getchrpos()) {
		return 1;
	    } else if(chrpos == t.getchrpos()){
		return 0;
	    } else {
		return -1;
	    }
	} else if (chr > t.getchr()) {
	    return 1;
	} else {
	    return -1;
	}
    }

    public boolean equals(SortableSNP t) {
	if (chr == t.getchr() && chr == t.getchrpos()) {
	    return true;
	} else {
	    return false;
	}
    }
}
