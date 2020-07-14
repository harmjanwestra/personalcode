package nl.harmjanwestra.methylation;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.Console;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;

public class ListGSEStillToProcess {
	
	public static void main(String[] args) {
		ListGSEStillToProcess l = new ListGSEStillToProcess();
		
		if (args.length < 3) {
			System.out.println("usage: procf inputf outputf");
		} else {
			String prccfolder = args[0]; //"S:\\projects\\2018-methylation\\GPL13534_450k_OUT\\";
			String inputfolder = args[1]; //"S:\\projects\\2018-methylation\\GPL13534_450k_GEO\\";
			String outputfile = args[2]; // "";
			try {
				l.run(prccfolder, inputfolder, outputfile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void run(String processedfolder, String inputfolder, String outputfile) throws IOException {
		File[] proclist = new File(processedfolder).listFiles();
		HashSet<String> processedset = new HashSet<String>();
		for (File f : proclist) {
			if (f.isDirectory()) {
//				System.out.println(f.getName());
				processedset.add(f.getName());
			}
		}
		File[] unproclist = new File(inputfolder).listFiles();
		HashSet<String> GSEsToFix = new HashSet<>();
		TextFile out = new TextFile(outputfile, TextFile.W);
		for (File f : unproclist) {
			if (f.getName().endsWith(".tar")) {
				String name = f.getName().replaceAll("_RAW.tar", "");
				if (!processedset.contains(name)) {
					System.out.println(f.getName() + "\t" + Gpio.humanizeFileSize(f.length()));
					
					TarArchiveInputStream tarInput = new TarArchiveInputStream(new FileInputStream(f));
					TarArchiveEntry entry;
					while (null != (entry = tarInput.getNextTarEntry())) {
						System.out.println(entry.getName());
						if (entry.getName().toLowerCase().contains("idat")) {
							GSEsToFix.add(f.getName());
						}
					}
					
					
				}
			}
		}
		out.close();
	}
	
	public static void waitForEnter(String message, Object... args) {
		Console c = System.console();
		if (c != null) {
			// printf-like arguments
			if (message != null)
				c.format(message, args);
			c.format("\nPress ENTER to proceed.\n");
			c.readLine();
		}
	}
	
	
}
