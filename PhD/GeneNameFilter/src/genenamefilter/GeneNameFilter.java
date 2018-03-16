/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package genenamefilter;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author harmjan
 */
public class GeneNameFilter {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            TextFile tf = new TextFile("/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-07-25-GenesForEQTLs2.txt", TextFile.R);

            HashSet<String> selectedGenes = new HashSet<String>();
            selectedGenes.addAll(tf.readAsArrayList());
            tf.close();
            HashSet<String> detectedGenes = new HashSet<String>();

            TextFile geneExpressionFile = new TextFile("/Volumes/iSnackHD/Data/GeneticalGenomicsDatasets/SatVatLiverMuscleBloodHT12/Blood-SAT-VAT-Liver-MuscleHT12CombinedExpressionData/BloodSATVATLiverMuscleHT12ProbesCentered.txt.50PCAsOverSamplesRemoved.txt.TriTyperFormat.txt.gz", TextFile.R);
            TextFile out = new TextFile("/Volumes/iSnackHD/Data/Projects/KarinFransen/2012-07-25-ProbesForGenesForEQTLs.txt", TextFile.W);

            String[] elems = geneExpressionFile.readLineElems(TextFile.tab); // header
            elems = geneExpressionFile.readLineElems(TextFile.tab);
            while (elems != null) {
                String probe = elems[0];
                String geneName = elems[7];

                String[] splittable = geneName.split("///");
                if (splittable.length > 1) {
                    for (int i = 0; i < splittable.length; i++) {
                        geneName = splittable[i];
                        System.out.println(geneName);
                        if (selectedGenes.contains(geneName)) {
                            out.writeln(probe);
                            detectedGenes.add(geneName);
                        }
                    }
                } else {
                    System.out.println(geneName);
                    if (selectedGenes.contains(geneName)) {
                        out.writeln(probe);
                        detectedGenes.add(geneName);
                    }
                }


                elems = geneExpressionFile.readLineElems(TextFile.tab);
            }

            geneExpressionFile.close();
            out.close();
            System.out.println("----");
            String[] selectedGenesArr = selectedGenes.toArray(new String[0]);
            for(String s: selectedGenesArr){
                if(!detectedGenes.contains(s)){
                    System.out.println(s);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();

        }

    }
}
