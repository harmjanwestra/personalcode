package nl.harmjanwestra.playground.conv;

import umcg.genetica.io.trityper.converters.PlinkDosageToTriTyper;

import java.io.IOException;

public class ConvertPlinkDosageToTT {

    public static void main(String[] args) {
        ConvertPlinkDosageToTT t = new ConvertPlinkDosageToTT();
        if (args.length < 3) {
            System.out.println("Usage: famfile file output");
        } else {
            try {
                t.run(args[0], args[1], args[2]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void run(String fam, String file, String output) throws IOException {
        PlinkDosageToTriTyper t = new PlinkDosageToTriTyper();
        t.convert(fam, file, output);
    }
}
