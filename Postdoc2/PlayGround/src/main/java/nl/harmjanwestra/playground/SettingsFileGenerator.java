package nl.harmjanwestra.playground;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class SettingsFileGenerator {

    public static void main(String[] args) {

        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-08-30-Titration\\2019-07-08-metaqtlsettings-EURandAFR-trans.xml";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-08-30-Titration\\SETTINGS-output\\trans-";
        String order = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-08-30-Titration\\datasetToSampleSize.txt";


        SettingsFileGenerator g = new SettingsFileGenerator();
        try {
//			g.run(input, order, output);

            String[] tissues = new String[]{
                    "cortex", "cortex-mixup", "cerebellum", "hippocampus", "basalganglia"
            };
            for (String tissue : tissues) {
                String datasetDefinition = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-14-mixup\\" + tissue + "-combos.txt";
                String platform = "gencode.v32";
                String annotation = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz";
                output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-02-14-mixup\\output\\" + tissue + "-combos.xml";
                g.generateDatasetRows(datasetDefinition, platform, annotation, output);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void generateDatasetRows(String datasetDefinition, String platform, String annotation, String output) throws IOException {
        TextFile outf = new TextFile(output, TextFile.W);
        outf.writeln("\t<datasets>");
        outf.writeln();

        TextFile in = new TextFile(datasetDefinition, TextFile.R);
        in.readLine();
        String[] elems = in.readLineElems(TextFile.tab);
        while (elems != null) {


            outf.writeln("\t\t<dataset>");

            outf.writeln("\t\t\t<name>" + elems[0] + "</name>");
            outf.writeln("\t\t\t<location>" + elems[2] + "</location>");
            outf.writeln("\t\t\t<genometoexpressioncoupling>" + elems[3] + "</genometoexpressioncoupling>");
            outf.writeln("\t\t\t<expressiondata>" + elems[1] + "</expressiondata>");
            outf.writeln("\t\t\t<probeannotation>" + annotation + "</probeannotation>");
            outf.writeln("\t\t\t<expressionplatform>" + platform + "</expressionplatform>");
            outf.writeln("\t\t\t<covariates/>");
            outf.writeln("\t\t\t<quantilenormalize>false</quantilenormalize>");
            outf.writeln("\t\t\t<logtransform>false</logtransform>");

            outf.writeln("\t\t</dataset>");
            outf.writeln();

            elems = in.readLineElems(TextFile.tab);
        }
        in.close();
        outf.writeln("\t</datasets>");
        outf.writeln("</settings>");
        outf.close();

    }


    class Dataset {
        String name;
        String location;
        String gte;
        String exp;
        String probe;
        String platform;
    }

    public void run(String input, String order, String outdir) throws IOException {

        TextFile tf = new TextFile(order, TextFile.R);
        ArrayList<String> list = tf.readAsArrayList();
        ArrayList<String> listtmp = tf.readAsArrayList();
        for (String s : list) {
            listtmp.add(s.split("\t")[0]);
        }
        tf.close();
        list = listtmp;

        ArrayList<String> lines = new ArrayList<>();
        ArrayList<Dataset> datasets = new ArrayList<>();

        TextFile in = new TextFile(input, TextFile.R);
        String ln = in.readLine();
        boolean dsblock = false;
        boolean dsopen = false;
        Dataset currentDataset = null;
        while (ln != null) {
            String lntmp = ln;
            while (lntmp.startsWith(" ")) {
                lntmp = lntmp.substring(1);
            }
            while (lntmp.startsWith("\t")) {
                lntmp = lntmp.substring(1);
            }

            if (lntmp.startsWith("<datasets>")) {
                dsblock = true;
            }


            if (!dsblock) {
                lines.add(ln);
            } else {
                if (lntmp.startsWith("<dataset>")) {
                    dsopen = true;
                    currentDataset = new Dataset();
                } else if (lntmp.startsWith("</dataset>")) {
                    datasets.add(currentDataset);
                    dsopen = false;
                } else if (dsopen) {
                    if (lntmp.startsWith("<name>")) {
                        currentDataset.name = getValue(lntmp);
                    } else if (lntmp.startsWith("<location>")) {
                        currentDataset.location = getValue(lntmp);
                    } else if (lntmp.startsWith("<genometoexpressioncoupling>")) {
                        currentDataset.gte = getValue(lntmp);
                    } else if (lntmp.startsWith("<expressiondata>")) {
                        currentDataset.exp = getValue(lntmp);
                    } else if (lntmp.startsWith("<probeannotation>")) {
                        currentDataset.probe = getValue(lntmp);
                    } else if (lntmp.startsWith("<expressionplatform>")) {
                        currentDataset.platform = getValue(lntmp);
                    }
                }

            }

            ln = in.readLine();
        }


        // order datasets according to datasetnames in "list"
        HashMap<String, Integer> dsIdx = new HashMap<String, Integer>();
        for (int i = 0; i < list.size(); i++) {
            dsIdx.put(list.get(i), i);
        }

        Dataset[] dsordered = new Dataset[datasets.size()];
        for (Dataset d : datasets) {
            Integer id = dsIdx.get(d.name);
            if (id == null) {
                System.out.println("Could not find dataset: " + d.name);
            } else {
                dsordered[id] = d;
            }

        }
        System.out.println(datasets.size() + " datasets :D");

        for (int i = 1; i < datasets.size(); i++) {

            int numToselect = i + 1;
            TextFile outf = new TextFile(outdir + numToselect + "-settings.xml", TextFile.W);

            for (String line : lines) {

                String lntmp = line;
                while (lntmp.startsWith(" ")) {
                    lntmp = lntmp.substring(1);
                }
                while (lntmp.startsWith("\t")) {
                    lntmp = lntmp.substring(1);
                }

                if (lntmp.startsWith("<outputdirectory>")) {
                    // replace outputdir, if needed
                    String outdirvar = getValue(lntmp);
                    outdirvar = outdirvar + numToselect + "/";
                    outf.writeln("<outputdirectory>" + outdirvar + "</outputdirectory>");
                } else {
                    outf.writeln(line);
                }
            }

            outf.writeln("\t<datasets>");
            outf.writeln();
            for (int j = 0; j < numToselect; j++) {
                outf.writeln("\t\t<dataset>");

                outf.writeln("\t\t\t<name>" + dsordered[j].name + "</name>");
                outf.writeln("\t\t\t<location>" + dsordered[j].location + "</location>");
                outf.writeln("\t\t\t<genometoexpressioncoupling>" + dsordered[j].gte + "</genometoexpressioncoupling>");
                outf.writeln("\t\t\t<expressiondata>" + dsordered[j].exp + "</expressiondata>");
                outf.writeln("\t\t\t<probeannotation>" + dsordered[j].probe + "</probeannotation>");
                outf.writeln("\t\t\t<expressionplatform>" + dsordered[j].platform + "</expressionplatform>");
                outf.writeln("\t\t\t<covariates/>");
                outf.writeln("\t\t\t<quantilenormalize>false</quantilenormalize>");
                outf.writeln("\t\t\t<logtransform>false</logtransform>");

                outf.writeln("\t\t</dataset>");
                outf.writeln();
            }
            outf.writeln("\t</datasets>");
            outf.writeln("</settings>");
            outf.close();
        }

    }

    private String getValue(String lntmp) {
        int open = lntmp.indexOf(">") + 1;
        int close = lntmp.lastIndexOf("<");


        if (close < open) {
            System.out.println("??");
        }
        String val = lntmp.substring(open, close);
        return val;
    }

}
