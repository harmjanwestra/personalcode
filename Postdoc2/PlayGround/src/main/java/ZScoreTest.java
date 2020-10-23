import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

public class ZScoreTest {

    public static void main(String[] args) {


        System.out.println(Double.MIN_NORMAL);


        System.out.println(ZScores.zToP(16.61));
        System.out.println(ZScores.zToP(8.82));
        System.out.println(ZScores.zToP(5.26));

        double r = ZScores.zToR(5.5572363,2843);
        System.out.println(r);
        r = ZScores.zToR(5.5572363,1229);
        System.out.println(r);
    }
}
