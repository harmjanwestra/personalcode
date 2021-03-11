package nl.harmjanwestra.playground.biogen.eqtlforestplot;

import nl.harmjanwestra.playground.biogen.mrvisualisation.MREvent;

import java.util.Objects;

public class ForestPlotEQTL implements Comparable<ForestPlotEQTL> {
    public String EA;
    public String NONEA;
    public String snp;
    public boolean isMeta;
    String gene;
    public String dataset;
    public int samplesize;

    double beta;
    double se;


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ForestPlotEQTL eqtl = (ForestPlotEQTL) o;
        return Double.compare(eqtl.beta, beta) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(beta);
    }


    @Override
    public int compareTo(ForestPlotEQTL eqtl) {
        if (this.equals(eqtl)) {
            return 0;
        }

        if (this.beta > eqtl.beta) {
            return 1;
        } else {
            return -1;
        }


    }
}
