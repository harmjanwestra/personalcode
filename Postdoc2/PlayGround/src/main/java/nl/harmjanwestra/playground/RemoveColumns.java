package nl.harmjanwestra.playground;

import nl.harmjanwestra.utilities.legacy.genetica.io.text.TextFile;
import org.apache.commons.io.comparator.NameFileComparator;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

public class RemoveColumns {


    public static void main(String[] args) {
        String p = "D:\\Sync\\SyncThing\\Postdoc\\2016-03-RAT1D-Finemapping\\Data\\2017-08-16-Reimpute4Filtered\\RA-T1D-Finemap-SummaryStats\\exhaustive\\";
        String p2 = "D:\\Code\\git\\hmsprojects\\RA-T1D-Finemap-SummaryStats\\exhaustive\\";


        HashSet<String> colsToRemove = new HashSet<>();
        colsToRemove.add("#Chr1");
        colsToRemove.add("Pos1");
        colsToRemove.add("Id1");
        colsToRemove.add("Chr2");
        colsToRemove.add("Pos2");
        colsToRemove.add("Id2");
        colsToRemove.add("Alleles1");
        colsToRemove.add("Alleles2");
        colsToRemove.add("MinorAllele1");
        colsToRemove.add("MinorAllele2");
        colsToRemove.add("ImputationQualScore");
        colsToRemove.add("ImputationQualScore2");
        colsToRemove.add("MAF1");
        colsToRemove.add("AFCases1");
        colsToRemove.add("AFControls1");
        colsToRemove.add("MAF2");
        colsToRemove.add("AFCases2");
        colsToRemove.add("AFControls2");
        colsToRemove.add("HWEP");
        colsToRemove.add("HWEP2");
        colsToRemove.add("Distance");
        colsToRemove.add("LD(D')");
        colsToRemove.add("DevianceNull");
        colsToRemove.add("DevianceGeno");
        colsToRemove.add("DfNull");
        colsToRemove.add("DfAlt");
        colsToRemove.add("DiffDf");
        colsToRemove.add("Beta(Genotype)");
        colsToRemove.add("SE(Genotype)");
        colsToRemove.add("OR");
        colsToRemove.add("OR-Hi");
        colsToRemove.add("OR-Lo");
        colsToRemove.add("Pval");

        RemoveColumns c = new RemoveColumns();
        try {
            c.run(p, p2, colsToRemove);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String loc, String locout, HashSet<String> colsToRemove) throws IOException {

        File folder = new File(loc);
        System.out.println(loc + "\texists: " + folder.exists());
        if (!folder.exists()) {
            System.exit(-1);
        }
        File[] listOfFiles = folder.listFiles();

        // sort files
        Arrays.sort(listOfFiles, NameFileComparator.NAME_COMPARATOR);
        for (File file : listOfFiles) {


            TextFile tf = new TextFile(file, TextFile.R);
            TextFile out = new TextFile(locout + file.getName(), TextFile.W);
            System.out.println("Parsing: " + file);
            System.out.println("out: " + locout + file.getName());
            String[] heade = tf.readLineElems(TextFile.tab);
            boolean[] keepcol = new boolean[heade.length];
            String outheader = "";
            boolean firstcol = true;
            for (int i = 0; i < heade.length; i++) {
                if (!colsToRemove.contains(heade[i])) {
                    keepcol[i] = true;
                    if (firstcol) {
                        outheader = heade[i];
                        firstcol = false;
                    } else {
                        outheader += "\t" + heade[i];
                    }
                }
            }

            out.writeln(outheader);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                firstcol = true;
                String outln = "";
                for (int i = 0; i < elems.length; i++) {
                    if (keepcol[i]) {
                        if (firstcol) {
                            outln = elems[i];
                            firstcol = false;
                        } else {
                            outln += "\t" + elems[i];
                        }
                    }
                }
                out.writeln(outln);
                elems = tf.readLineElems(TextFile.tab);
            }

            tf.close();

        }

    }
}
