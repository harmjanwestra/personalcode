/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.zscoretables;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class ZScoreTableGeneAndSNPCorrelator {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	// TODO code application logic here

	try {
	    TextFile tf = new TextFile("", TextFile.R);

	    tf.close();
	} catch (IOException e) {
	    e.printStackTrace();
	}
    }
}
