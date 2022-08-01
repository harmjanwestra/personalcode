package nl.harmjanwestra.playground.biogen.freeze2dot1.vcf;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import nl.harmjanwestra.playground.legacy.vcf.*;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.broad.tribble.util.popgen.HardyWeinbergCalculation;
import org.checkerframework.checker.units.qual.A;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.HWE;
import umcg.genetica.util.Primitives;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;

public class VariantQual {

    public static void main(String[] args) {
//nl.harmjanwestra.playground.biogen.freeze2dot1.vcf.VariantQual
        VariantQual q = new VariantQual();
        q.run(args[0], args[1]);

    }

    //    public void run(String vcf, String output, String sampleToDataset) {
    public void run(String vcf, String output) {

        try {

//            HashMap<String, String> sampleMap = new HashMap<>();
//            TextFile tf1 = new TextFile(sampleToDataset, TextFile.R);
//            String[] selems = tf1.readLineElems(TextFile.tab);
//            while (selems != null) {
//                sampleMap.put(selems[0], selems[1]);
//                selems = tf1.readLineElems(TextFile.tab);
//            }
//            tf1.close();

            TextFile tf = new TextFile(vcf, TextFile.R);

            int[] sampleIndex = null;
            int dsctr = 0;
            HashMap<String, Integer> dsMap = new HashMap<>();
            String line = tf.readLine();
            TextFile out = new TextFile(output, TextFile.W);
            out.writeln("ID\tDosageRsq\tN\tCallrate\tMAF\tHWE\tMeanDelta\tSDDelta\tCVDelta");
            int ctr = 0;
            while (line != null) {
                if (line.startsWith("#CHROM")) {
                    String[] elems = line.split("\t");
                    sampleIndex = new int[elems.length];
//                    for (int i = 9; i < elems.length; i++) {
//                        String ds = sampleMap.get(elems[i]);
//                        Integer id = dsMap.get(ds);
//                        if (id == null) {
//                            id = dsctr;
//                            dsMap.put(ds, id);
//                            dsctr++;
//                        }
//                        sampleIndex[i] = id;
//                    }
                } else if (!line.startsWith("#")) {
                    String[] elems = line.split("\t");
                    String id = elems[2];
                    ArrayList<Double> genotypes = new ArrayList<>();
                    ArrayList<Double> deltas = new ArrayList<>();
                    ArrayList<Double> dosages = new ArrayList<>();
                    double callrate = 0;


//                    ArrayList<ArrayList<Double>> gtsPerDs = new ArrayList<>();
//                    ArrayList<ArrayList<Double>> dosagesPerDs = new ArrayList<>();
//                    for (int i = 0; i < dsMap.size(); i++) {
//                        gtsPerDs.add(new ArrayList<>());
//                        dosagesPerDs.add(new ArrayList<>());
//                    }
                    for (int i = 9; i < elems.length; i++) {
                        String val = elems[i];
                        String[] gtelems = val.split(":");
                        String gt = gtelems[0];
                        String ds = gtelems[1];
                        if (!gt.equals("./.")) {
                            String[] gtelems2 = gt.split("/");
                            double gti = Integer.parseInt(gtelems2[0]) + Integer.parseInt(gtelems2[1]);

                            double dosage = Double.parseDouble(ds);
                            genotypes.add(gti);
                            dosages.add(dosage);
                            double delta = gti - dosage;
                            deltas.add(delta);
                            int dsid = sampleIndex[i];
//                            gtsPerDs.get(dsid).add(gti);
//                            dosagesPerDs.get(dsid).add(dosage);
                            callrate++;
                        }
                    }

                    callrate /= (elems.length - 9);

                    PearsonsCorrelation correlation = new PearsonsCorrelation();
                    double r = correlation.correlation(Primitives.toPrimitiveArr(genotypes), Primitives.toPrimitiveArr(dosages));


                    double rsq = r * r;
                    double maf = maf(genotypes);
                    double hwe = hwe(genotypes);
                    double meanDelta = JSci.maths.ArrayMath.mean(Primitives.toPrimitiveArr(deltas));
                    double sdDelta = JSci.maths.ArrayMath.standardDeviation(Primitives.toPrimitiveArr(deltas));
                    double cv = (meanDelta / sdDelta);
                    if (maf == 0 || Double.isNaN(rsq)) {
                        rsq = 0;
                        cv = 0;
                    }

//                    ArrayList<Double> mafs = new ArrayList<>();
//                    ArrayList<Double> hwes = new ArrayList<>();
//                    for (int i = 0; i < dsMap.size(); i++) {
//                        mafs.add(maf(gtsPerDs.get(i)));
//                        mafs.add(hwe(gtsPerDs.get(i)));

                    if (Double.isNaN(rsq)) {
                        for (int i = 0; i < genotypes.size(); i++) {
                            System.out.println(genotypes.get(i) + "\t" + dosages.get(i));
                        }
                        System.out.println();
                        System.out.println(id + "\t" + rsq + "\t" + genotypes.size() + "\t" + callrate + "\t" + maf + "\t" + hwe + "\t" + meanDelta + "\t" + sdDelta + "\t" + cv);


                        System.exit(0);
                    }

                    out.writeln(id + "\t" + genotypes.size() + "\t" + rsq + "\t" + callrate + "\t" + maf + "\t" + hwe + "\t" + meanDelta + "\t" + sdDelta + "\t" + cv);
                }


                ctr++;
                if (ctr % 1000 == 0) {
                    System.out.print(ctr + " lines parsed\r");
                }
                line = tf.readLine();
            }

            System.out.print(ctr + " lines parsed\n");
            tf.close();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    private double hwe(ArrayList<Double> gts) {
        int obsHomA = 0;
        int obsHet = 0;
        int obsHomB = 0;
        for (int d = 0; d < gts.size(); d++) {
            double gti = gts.get(d);
            if (gti == 0) {
                obsHomA++;
            } else if (gti == 1) {
                obsHet++;
            } else {
                obsHomB++;
            }
        }
        return HWE.calculateExactHWEPValue(obsHet, obsHomA, obsHomB);
    }

    private double maf(ArrayList<Double> gts) {

        double[] freqs = new double[2];
        for (int d = 0; d < gts.size(); d++) {
            int gt = (int) Math.round(gts.get(d));
            if (gt == 0) {
                freqs[0] += 2;
            } else if (gt == 1) {
                freqs[0] += 1;
                freqs[1] += 1;
            } else {
                freqs[1] += 2;
            }
        }
        freqs[0] /= (gts.size() * 2);
        freqs[1] /= (gts.size() * 2);
        if (freqs[0] < freqs[1]) {
            return freqs[0];
        } else {
            return freqs[1];
        }
    }


}
