package nl.harmjanwestra.parquetreader;

public class Triple<L, M, R> {
    private final L left;
    private final M middle;
    private final R right;

    public Triple(L left, M middle, R right) {
        this.left = left;
        this.middle = middle;
        this.right = right;
    }

    public L getLeft() {
        return left;
    }

    public M getMiddle() {
        return middle;
    }

    public R getRight() {
        return right;
    }

    @Override
    public int hashCode() {
        return left.hashCode() ^ middle.hashCode() ^ right.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (o == null) {
            return false;
        }
        if (!(o instanceof Triple)) {
            return false;
        }
        Triple pairo = (Triple) o;
        return this.left.equals(pairo.getLeft())
                && this.right.equals(pairo.getRight()) && this.middle.equals(pairo.getMiddle());
    }

    @Override
    public String toString() {
        return "Triple{" + "left=" + left + ", middle=" + middle + ", right=" + right + '}';
    }

}
