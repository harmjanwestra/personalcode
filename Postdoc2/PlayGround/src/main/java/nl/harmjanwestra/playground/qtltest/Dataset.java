package nl.harmjanwestra.playground.qtltest;

import umcg.genetica.util.Primitives;

import java.util.ArrayList;
import java.util.Objects;

public class Dataset implements Comparable<Dataset> {
    public String name;
    int[] genotypeIds;
    int[] expressionIds;


    ArrayList<Integer> DNAIds = new ArrayList<>();
    ArrayList<Integer> RNAIds = new ArrayList<>();

    public void append(Integer dnaId, Integer rnaId) {
        DNAIds.add(dnaId);
        RNAIds.add(rnaId);
    }

    public void toArr() {
        genotypeIds = Primitives.toPrimitiveArr(DNAIds);
        expressionIds = Primitives.toPrimitiveArr(RNAIds);
        DNAIds = null;
        RNAIds = null;
    }

    public double[] select(double[] vals, int[] ids) {
        double[] data = new double[ids.length];
        for (int g = 0; g < ids.length; g++) {
            data[g] = vals[ids[g]];
        }
        return data;
    }

    public double[] select(byte[] vals, int[] ids) {
        double[] data = new double[ids.length];
        for (int g = 0; g < ids.length; g++) {
            data[g] = vals[ids[g]];
        }
        return data;
    }


    public double[] select(double[][] vals, int[] ids) {
        double[] data = new double[ids.length];
        for (int g = 0; g < ids.length; g++) {
            data[g] = vals[ids[g]][0];
        }
        return data;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Dataset dataset = (Dataset) o;
        return Objects.equals(name, dataset.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name);
    }

    @Override
    public int compareTo(Dataset dataset) {
        if (this.equals(dataset)) {
            return 0;
        }
        return this.name.compareTo(dataset.name);

    }
}
