package nl.harmjanwestra.playground.biogen.mrvisualisation;

import java.util.Objects;

public class MREvent implements Comparable<MREvent> {

    public String EA;
    public String NONEA;
    public String snp;
    public String outcome;
    public double pp4bon;
    public double pp4def;
    public String pWR;
    String gene;
    double betaExposure;
    double seExposure;
    double betaOutcome;
    double seOutcome;
    double waldratio;
    double seWaldratio;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        MREvent mrEvent = (MREvent) o;
        return Double.compare(mrEvent.waldratio, waldratio) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(waldratio);
    }


    @Override
    public int compareTo(MREvent mrEvent) {
        if (this.equals(mrEvent)) {
            return 0;
        }

        if (this.waldratio > mrEvent.waldratio) {
            return 1;
        } else {
            return -1;
        }


    }
}
