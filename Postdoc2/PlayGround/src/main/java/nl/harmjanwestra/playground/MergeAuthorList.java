package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

public class MergeAuthorList {
	
	public static void main(String[] args) {
		String aut = "D:\\Sync\\Dropbox\\FineMap\\2018-06-EditorComments\\Authors\\Authorlist.txt";
		String aff = "D:\\Sync\\Dropbox\\FineMap\\2018-06-EditorComments\\Authors\\AffiliationList.txt";
		String output = "D:\\Sync\\Dropbox\\FineMap\\2018-06-EditorComments\\Authors\\Affiliations.html";
		
		MergeAuthorList m = new MergeAuthorList();
		try {
			m.run(aut, aff, output);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public void run(String aut, String aff, String outputfile) throws IOException {
		TextFile tf2 = new TextFile(aff, TextFile.R);
		String[] elems = tf2.readLineElems(TextFile.tab);
		HashMap<String, String> map = new HashMap<String, String>();
		HashMap<String, String> reversemap = new HashMap<String, String>();
		
		while (elems != null) {
			if (elems.length >= 2) {
				map.put(elems[0], elems[1]); // id --> affiliation
				reversemap.put(elems[1], elems[0]); // affiliation --> id
			}
			elems = tf2.readLineElems(TextFile.tab);
		}
		tf2.close();
		System.out.println(map.size() + " affiliations ");
		
		
		TextFile tf = new TextFile(aut, TextFile.R);
		elems = tf.readLineElems(TextFile.tab);
		
		HashMap<String, Integer> refToId = new HashMap<>();
		ArrayList<String> refsout = new ArrayList<String>();
		int ctr = 1;
		ArrayList<String> authorOutput = new ArrayList<>();
		while (elems != null) {
			
			if (elems.length >= 2) {
				String name = elems[0];
				String refs = elems[1];
				String[] refelems = refs.split(",");
				int[] ids = new int[refelems.length];
				for (int r = 0; r < refelems.length; r++) {
					// get name
					String ref = refelems[r];
					String affiliation = map.get(ref);
					if (refToId.containsKey(affiliation)) {
						ids[r] = refToId.get(affiliation);
						
					} else {
						refToId.put(affiliation, ctr);
						refsout.add(affiliation);
						ids[r] = ctr;
						ctr++;
					}
				}
				Arrays.sort(ids);
				String out = "<div class='author'>" + name + "<sup>" + Strings.concat(ids, Strings.comma) + "</sup></div>";
				authorOutput.add(out);
			}
			elems = tf.readLineElems(TextFile.tab);
		}
		
		tf.close();
		
		TextFile output = new TextFile(outputfile, TextFile.W);
		output.writeln("<!DOCTYPE html>");
		output.writeln("<head>" +
				"<title>authorlist</title>" +
				"<style>" +
				".author { display:inline; } " +
				"" +
				"</style>" +
				"</head><body>");
		output.writeln("<div id='authors'>" + Strings.concat(authorOutput, Pattern.compile(", ")) + "</div>");
		output.writeln("<br/>");
		output.writeln("<div id='references'>");
		for (int i = 0; i < refsout.size(); i++) {
			
			String out = "<div class='reference'><sup>" + (i + 1) + "</sup> " + refsout.get(i) + "</div><br/>";
			output.writeln(out);
		}
		output.writeln("</div></body></html>");
		output.close();
	}
}
