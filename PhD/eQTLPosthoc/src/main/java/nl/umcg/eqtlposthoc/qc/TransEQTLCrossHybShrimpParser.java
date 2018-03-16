/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.eqtlposthoc.qc;

import java.io.IOException;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class TransEQTLCrossHybShrimpParser {

    /**
     * @param args the command line arguments
     */
    public void run(String dir) {
	// TODO code application logic here

	

	try {
	    String[] fileList = Gpio.getListOfFiles(dir, "out");

	    for (String file : fileList) {
		TextFile tf = new TextFile(file, TextFile.R);

		String line = tf.readLine();
		boolean mappingfound = false;
		while (line != null) {

		    if (line.startsWith(">")) {
			mappingfound = true;
			String[] lnElems = Strings.tab.split(line);
			// if we can parse the edit string as an integer directly, it must be a perfect match...

			String editstring = lnElems[lnElems.length - 1];
			Integer identity = null;
			try {
			    Integer.parseInt(editstring);
			} catch (NumberFormatException e) {
			}

			if (identity != null) {
			    // print the identity found...
			    System.out.println(line + "\tPerfect match: " + identity);
			} else {
			    // we have to tear the edit string apart, soul by soul.

			    boolean lastcharacterwasaninteger = false;
			    int numEdits = 0;
			    boolean[] positionsparsed = new boolean[editstring.length()];
			    boolean[] positionisinteger = new boolean[editstring.length()];

			    for (int c = 0; c < editstring.length(); c++) {
				char p = editstring.charAt(c);
				if (p == 'A' || p == 'C' || p == 'T' || p == 'G') {
				    // we found an edit
				    numEdits++;
				    positionsparsed[c] = true;
				} else {
				    // it must be a gap, or identity...

				    try {
					Integer g = Integer.parseInt("" + p);
					positionisinteger[c] = true;
				    } catch (NumberFormatException e) {
				    }

				}
			    }

			    // parse the integers in the string
			    int numIdentity = 0;
			    for (int c = 0; c < editstring.length(); c++) {
				char p = editstring.charAt(c);


				if (!positionsparsed[c]) {
				    boolean nextcharisinteger = false;
				    if (positionisinteger[c]) {
					if (c + 1 < editstring.length()) {
					    if (positionisinteger[c + 1]) {
						nextcharisinteger = true;
						// integer is two digits..
						String integerconcat = "" + editstring.charAt(c) + editstring.charAt(c + 1);
						Integer eventualInteger = Integer.parseInt(integerconcat);
						numIdentity += eventualInteger;
						positionsparsed[c] = true;
						positionsparsed[c + 1] = true;
					    }
					}
					if (!nextcharisinteger) {
					    Integer eventualInteger = Integer.parseInt("" + editstring.charAt(c));
					    numIdentity += eventualInteger;
					    positionsparsed[c] = true;
					}
				    }

				}

			    }
			    System.out.println(line + "\tId: " + numIdentity + "\tEd: " + numEdits + "\tSum: " + (numEdits + numIdentity));
			}
		    }

		    line = tf.readLine();
		}

		if (!mappingfound) {
//		    System.out.println(file + "\tSNP-Probe combination has no mapping.");
		}

		tf.close();
	    }
	} catch (IOException e) {
	    e.printStackTrace();
	}

    }

    public TransEQTLCrossHybShrimpParser() {
    }
}
