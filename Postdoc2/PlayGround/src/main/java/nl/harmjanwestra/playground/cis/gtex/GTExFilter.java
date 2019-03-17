package nl.harmjanwestra.playground.cis.gtex;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.*;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

public class GTExFilter {

    public static void main(String[] args) {
        if (args.length < 3) {
            System.out.println("Usage: combofile gtextarfile outdir");
        } else {
            GTExFilter f = new GTExFilter();
            try {
                f.run(args[0], args[1], args[2]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String combofile, String gtextar, String outdir) throws IOException {
        HashSet<String> combos = loadCombos(combofile);

        System.out.println(combos.size() + " combos loaded from " + combofile);


        try {

//			{
//				InputStream fi = new FileInputStream(gtextar);
//				InputStream bi = new BufferedInputStream(fi);
//				InputStream gzi = new GzipCompressorInputStream(bi);
//				TarArchiveInputStream fin = new TarArchiveInputStream(gzi);
//				TarArchiveEntry entry;
//				while ((entry = fin.getNextTarEntry()) != null) {
//
//					System.out.println(entry.getName());
//				}
//			}

//			System.exit(-1);
            InputStream fi = new FileInputStream(gtextar);
            InputStream bi = null;
            if (gtextar.endsWith(".gz")) {
                bi = new GZIPInputStream(new BufferedInputStream(fi));
            } else {
                bi = new BufferedInputStream(fi);
            }

            InputStream gzi = new GzipCompressorInputStream(bi);
            TarArchiveInputStream fin = new TarArchiveInputStream(gzi);
            TarArchiveEntry entry;
            while ((entry = fin.getNextTarEntry()) != null) {
                if (entry.isDirectory()) {
                    continue;
                }
                System.out.println("For File = " + entry.getName());
                BufferedReader br = null;
                if (entry.getName().endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fin))); // Read directly from tarInput
                } else {
                    br = new BufferedReader(new InputStreamReader(fin)); // Read directly from tarInput
                }

                String name = entry.getName();
                String[] nameelems = name.split("/");
                if (nameelems.length > 1) {
                    name = nameelems[nameelems.length - 1];
                }
                System.out.println("Found: " + name);
                TextFile out = new TextFile(outdir + name, TextFile.W);

                String line = br.readLine();
                System.out.println(line);
                out.writeln(line);

                int ctr = 0;
                int written = 0;
                line = br.readLine();
                while (line != null) {

                    String[] elems = line.split("\t");
                    String gene = elems[0].split("\\.")[0];
                    String[] snp = elems[1].split("_");
//					System.out.println(line);
                    String pos = snp[0] + "_" + snp[1];
                    if (combos.contains(pos + "_" + gene)) {
                        out.writeln(Strings.concat(elems, Strings.tab));
                        written++;
                    }

                    ctr++;
                    if (ctr % 100000 == 0) {
                        System.out.print(ctr + " lines read from " + name + ". " + written + " written\r");
//						System.exit(-1);
                    }
                    line = br.readLine();
                }
                out.close();
                System.out.println("Done reading " + name + " " + written + " written, " + ctr + " read");

            }
            gzi.close();


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private HashSet<String> loadCombos(String combofile) throws IOException {
        HashSet<String> combos = new HashSet<>();
        TextFile tf = new TextFile(combofile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[2] + "_" + elems[3];
            String gene = elems[4];
            combos.add(snp + "_" + gene);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        return combos;
    }
}
