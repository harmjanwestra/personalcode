/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class SampleReplace {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        String samplecouple = "/Volumes/iSnackHD/TonuData/New_Vcodes_Old_Pcodes.txt";
        String gtefile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/GenotypeExpressionCoupling.txt";
        String expfile = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/ExpressionData/ExpressionData.txt.40PCAsOverSamplesRemoved-GeneticVectorsNotRemoved.txt";
        String oldInds = "/Volumes/iSnackHD/TonuData/TT/Individuals.txt";
        String newInds = "/Volumes/iSnackHD/TonuData/TT/IndividualsPCodes.txt";
        String oldPheno = "/Volumes/iSnackHD/TonuData/TT/PhenotypeInformation.txt";
        String origPheno = "/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/EGCUT/Hap2ImputedGenotypes/PhenotypeInformation.txt";
        String newPheno = "/Volumes/iSnackHD/TonuData/TT/PhenotypeInformationPCodes.txt";

        try {

            HashSet<String> expSamples = new HashSet<String>();
            TextFile exp = new TextFile(expfile, TextFile.R);
            String[] exelems = exp.readLineElems(TextFile.tab);
            expSamples.addAll(Arrays.asList(exelems));
            exp.close();

            System.out.println(expSamples.size() + " expression samples");

            HashMap<String, String> gte = new HashMap<String, String>();
            TextFile tfg = new TextFile(gtefile, TextFile.R);
            String[] gtel = tfg.readLineElems(TextFile.tab);
            HashSet<String> linkedSamples = new HashSet<String>();
            while (gtel != null) {
                if (expSamples.contains(gtel[1])) {
                    
                    linkedSamples.add(gtel[0]);
                }
                gtel = tfg.readLineElems(TextFile.tab);
            }
            tfg.close();

            System.out.println(linkedSamples.size() + " linked to a genotype");
            

            TextFile tf = new TextFile(samplecouple, TextFile.R);
            HashMap<String, String> dt = (HashMap<String, String>) tf.readAsHashMap(0, 1);
            tf.close();

            TextFile tf2 = new TextFile(oldInds, TextFile.R);
            TextFile tf2out = new TextFile(newInds, TextFile.W);
            String ln = tf2.readLine();
            while (ln != null) {
                String smp = dt.get(ln);
                if (smp != null) {
                    tf2out.writeln(smp);
                } else {
                    tf2out.writeln(ln);
                }
                ln = tf2.readLine();
            }
            tf2.close();
            tf2out.close();

            HashMap<String, String> phenoinfo = new HashMap<String, String>();

            TextFile tf4 = new TextFile(origPheno, TextFile.R);
            String[] telems = tf4.readLineElems(TextFile.tab);
            while (telems != null) {
                phenoinfo.put(telems[0], Strings.concat(telems, Strings.tab));
                telems = tf4.readLineElems(TextFile.tab);
            }
            tf4.close();

            TextFile tf3 = new TextFile(oldPheno, TextFile.R);
            TextFile tf3out = new TextFile(newPheno, TextFile.W);
            String[] elems = tf3.readLineElems(TextFile.tab);

            HashSet<String> availablesamples = new HashSet<String>();
            int ctr = 0;
            int ctr2 = 0;
            int ctr3 = 0;
            while (elems != null) {
                String smp = dt.get(elems[0]);
                availablesamples.add(smp);
                if (smp == null) {
                    System.out.println("Sample not in old data: " + smp);
                    ctr3++;
                    tf3out.writeln(Strings.concat(elems, Strings.tab));
                } else {
                    String pheno = phenoinfo.get(smp);
                    if (pheno == null) {
                        System.out.println("Sample: " + smp + " is in old data, but has no phenotype information: " + pheno);
                        ctr++;
                        elems = new String[]{smp, "control", "include", "unknown"};
                        tf3out.writeln(Strings.concat(elems, Strings.tab));
                    } else {
                        ctr2++;
                        tf3out.writeln(pheno);
                    }
                }
                elems = tf3.readLineElems(TextFile.tab);
            }
            System.out.println(ctr + "\t" + ctr2 + "\t" + ctr3);
            System.out.println((ctr2 + ctr));
            tf3.close();
            tf3out.close();

            int ctr4 = 0;
            for (String sample : linkedSamples) {
                if (availablesamples.contains(sample)) {
                    ctr4++;
                } else {
                    System.out.println(sample);
                }
            }
            System.out.println(ctr4 + " with expression sample");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
