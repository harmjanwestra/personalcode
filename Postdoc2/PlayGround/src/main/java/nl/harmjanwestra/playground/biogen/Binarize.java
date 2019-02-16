package nl.harmjanwestra.playground.biogen;

import nl.harmjanwestra.utilities.legacy.genetica.math.matrix.DoubleMatrixDataset;

import java.io.IOException;

public class Binarize {

    public static void main(String[] args){

        try {
            DoubleMatrixDataset<String, String> str = new DoubleMatrixDataset<>(args[0]);
            str.save(args[1]);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }
}
