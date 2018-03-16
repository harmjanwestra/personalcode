/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package noortjeubertar;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.regex.Pattern;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import umcg.genetica.io.Gpio;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class NoortjeUberTar {

    final static int BUFFER = 2048;

    /**
     * Command line arguments : argv[0]-----> Source tar.gz file. argv[1]----->
     * DestarInation directory.
     *
     */
    public static void main(String[] args) {
        if (args.length < 3) {
            System.out.println("Usage: NoortjeUberTar.jar tarfile query output");

        } else {
            try {

                NoortjeUberTar.run(args[0], args[1], args[2]);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        System.exit(0);

    }

    public static void run(String tarfile, String query, String outdir) throws FileNotFoundException, IOException {
        if (!outdir.endsWith("/")) {
            outdir += "/";
        }
        Gpio.createDir(outdir);

        /**
         * create a TarArchiveInputStream object. *
         */
        FileInputStream fin = new FileInputStream(tarfile);
        BufferedInputStream in = new BufferedInputStream(fin);
        GzipCompressorInputStream gzIn = new GzipCompressorInputStream(in);
        TarArchiveInputStream tarIn = new TarArchiveInputStream(gzIn);

        TarArchiveEntry entry = null;

        /**
         * Read the tar entries using the getNextEntry method *
         */
        System.out.println("Now extracting: " + tarfile);
        while ((entry = (TarArchiveEntry) tarIn.getNextEntry()) != null) {

            if (entry.getName().contains(query)) {
                System.out.println("Extracting: " + entry.getName());

                String[] elems = entry.getName().split("/");
                if (elems.length > 1) {
                    String tardir = Strings.concat(elems, Pattern.compile("/"), 0, elems.length-2);
                    
                    Gpio.createDir(outdir + elems[0]);
                }
//                /**
//                 * If the entry is a directory, create the directory. *
//                 */
//                if (entry.isDirectory()) {
//                    Gpio.createDir(outdir+entry.getName());
//                }
            }

        }

        /**
         * Close the input stream *
         */
        tarIn.close();


        fin = new FileInputStream(tarfile);
        in = new BufferedInputStream(fin);
        gzIn = new GzipCompressorInputStream(in);
        tarIn = new TarArchiveInputStream(gzIn);

        entry = null;

        /**
         * Read the tar entries using the getNextEntry method *
         */
        System.out.println("Now extracting: " + tarfile);
        while ((entry = (TarArchiveEntry) tarIn.getNextEntry()) != null) {

            if (entry.getName().contains(query)) {
                System.out.println("Extracting: " + entry.getName());


                /**
                 * If the entry is a directory, create the directory. *
                 */
                if (!entry.isDirectory()) {

                    /**
                     * If the entry is a file,write the decompressed file to the
                     * disk and close destination stream.
                     *
                     */
                    int count;
                    byte data[] = new byte[BUFFER];

                    FileOutputStream fos = new FileOutputStream(outdir + entry.getName());
                    BufferedOutputStream dest = new BufferedOutputStream(fos, BUFFER);
                    while ((count = tarIn.read(data, 0, BUFFER)) != -1) {
                        dest.write(data, 0, count);
                    }
                    dest.close();
                }
            }

        }

        /**
         * Close the input stream *
         */
        tarIn.close();
        System.out.println("untar completed successfully!!");
    }
}
