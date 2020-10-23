package nl.harmjanwestra.playground.biogen;

import com.itextpdf.text.DocumentException;
import umcg.genetica.graphics.Grid;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.ScatterplotPanel;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;

import java.io.IOException;

public class DSZscoreCompare {

    public static void main(String[] args) {
        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\eQTLDump-sortZ-FDR.txt.gz-Significant-0.05.txt.gz";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\eQTLDump-sortZ-FDR.txt.gz-Significant-0.05.png";
        input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\trans-eQTLDump.txt.gz";
        output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-08-11-MAFCompare\\trans-eQTLDump.png";
        DSZscoreCompare d = new DSZscoreCompare();
        try {
            d.plot(input, output);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (DocumentException e) {
            e.printStackTrace();
        }
    }

    public void plot(String input, String output) throws IOException, DocumentException {

        TextFile tf = new TextFile(input, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);

        ScatterplotPanel[] p = null;
        Grid g = new Grid(300, 300, 5, 5, 100, 100);
        while (elems != null) {
            String[] dsnames = elems[11].split(";");
            String[] dsZ = elems[12].split(";");
            double metaZ = Double.parseDouble(elems[10]);
            if (p == null) {
                p = new ScatterplotPanel[dsnames.length];
                for (int i = 0; i < p.length; i++) {
                    p[i] = new ScatterplotPanel(1, 1);
                    g.addPanel(p[i]);
                    p[i].setDataRange(new Range(-40, -40, 40, 40));
                }
            }
            for (int d = 0; d < dsnames.length; d++) {
                if (!dsnames[d].equals("-")) {
                    p[d].setTitle(dsnames[d]);
                    p[d].setLabels("Meta", dsnames[d]);
                    double v = Double.parseDouble(dsZ[d]);
                    p[d].addData(metaZ, v);
                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        g.draw(output);


    }
}
