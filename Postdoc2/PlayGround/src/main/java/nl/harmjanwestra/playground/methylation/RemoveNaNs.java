package nl.harmjanwestra.playground.methylation;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import nl.harmjanwestra.utilities.legacy.genetica.text.Strings;

import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicInteger;

public class RemoveNaNs {
    public static void main(String[] args) {
        RemoveNaNs s = new RemoveNaNs();
        try {
            s.run(args[0], args[1]);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String in, String out) throws IOException {
        TextFile tf = new TextFile(in, TextFile.R);
        TextFile tfo = new TextFile(out, TextFile.W);
        tfo.writeln(tf.readLine());
        String ln = tf.readLine();
        AtomicInteger nrw = new AtomicInteger();
        AtomicInteger nr = new AtomicInteger();

//        ExecutorService executor = Executors.newFixedThreadPool(32);

        while (ln != null) {

            String[] elems = Strings.tab.split(ln);

            int nrnan = 0;
            double sum = 0;
            int nrnotnan = 0;
            for (int i = 1; i < elems.length; i++) {
                Double d = Double.parseDouble(elems[i]);
                if (Double.isNaN(d)) {
                    nrnan++;
                } else {
                    nrnotnan++;
                    sum += d;
                }
            }
            double perc = (double) nrnan / (elems.length - 1);
            if (perc < 0.10) {
                double mean = sum / nrnotnan;
                for (int i = 1; i < elems.length; i++) {
                    Double d = Double.parseDouble(elems[i]);
                    if (Double.isNaN(d)) {
                        elems[i] = "" + mean;
                    }
                }
                try {
                    tfo.writelnsynced(Strings.concat(elems, Strings.tab));
                } catch (IOException e) {
                    e.printStackTrace();
                }
                nrw.getAndIncrement();
            }
            nr.getAndIncrement();

            if (nr.get() % 100 == 0) {
                System.out.print("\r" + nr + "\t" + nrw);
            }
            ln = tf.readLine();

        }




//        executor.shutdown();
//        while(!executor.isShutdown()){
//
//        }
        tf.close();
        tfo.close();
    }

    class Parser implements Runnable {

        private final String ln;

        private final AtomicInteger nrw;
        private final AtomicInteger nr;
        private final TextFile tfo;

        public Parser(String ln, AtomicInteger nr, AtomicInteger nrw, TextFile tfo) {
            this.ln = ln;
            this.nr = nr;
            this.nrw = nrw;
            this.tfo = tfo;
        }

        @Override
        public void run() {

        }
    }


}
