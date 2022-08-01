package nl.harmjanwestra.playground.biogen;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;

import java.io.IOException;

public class Bryois {

    public static void main(String[] args) {

        try {
            TextFile tf = new TextFile("D:\\transfer\\julienbryois2021_Astrocytes_replication_data.txt", TextFile.R);
            TextFile tfo = new TextFile("D:\\transfer\\julienbryois2021_Astrocytes_replication_data-out.txt", TextFile.W);
            String header = tf.readLine() + "\tZHJ\tbetaHJ\tseHJ\tseFromZ";
            tfo.writeln(header);
            String[] elems = tf.readLineElems(TextFile.tab);
            while (elems != null) {
                double p = Double.parseDouble(elems[10]);
                double maf = Double.parseDouble(elems[6]);
                double beta = Double.parseDouble(elems[11]);
                int n = Integer.parseInt(elems[13]);
                double z = ZScores.pToZTwoTailed(p);
                if (beta >= 0) {
                    z *= -1;
                }
                double[] betaAndSEEstimate = ZScores.zToBeta(z, maf, n);

                // z=beta/se
                double seFromZ = beta / z;

//                double rightside = ((z / beta) * (z / beta)) / ((n + (z * z)) / 2);
                // maf * (1.0D - maf) = ((z/beta)^2 / ((double)n + chi)/2)
                // x(1-x)= y
                // x-x^2 = y
                // - x^2 + x - y = 0
                // x = 1 + sqrt(1 - 4 * y) / 2

//                double derivedMAF = z / Math.sqrt()
//                System.out.println(rightside);
//                System.out.println(z + "\t" + maf + "\t" + mafderivedA + "\t" + mafderivedB);

                // maf * (1.0D - maf) = ((z/beta)^2 / ((double)n + chi)/2)

                tfo.writeln(Strings.concat(elems, Strings.tab) + "\t" + z + "\t" + betaAndSEEstimate[0] + "\t" + betaAndSEEstimate[1] + "\t" + seFromZ);
                System.out.println(p + "\t" + beta + "\t" + z + "\t" + betaAndSEEstimate[0] + "\t" + betaAndSEEstimate[1] + "\t" + seFromZ);
                System.exit(0);
                elems = tf.readLineElems(TextFile.tab);
            }
            tfo.close();
            tf.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static double[] zToBeta(double z, double maf, int n) {
        double chi = z * z;
        double a = 2.0D * maf * (1.0D - maf) * ((double) n + chi);
        double beta = z / Math.sqrt(a);
        double se = 1.0D / Math.sqrt(a);
        return new double[]{beta, se};
    }

//    public static double[] zToMAF(double z, double beta, int n) {
//        double chi = z * z;
//        //  (2.0D * maf * (1.0D - maf)) = a / ((double)n + chi)
//        //           maf * (1.0D - maf) = (a / ((double)n + chi)/2)
//        //           maf * (1.0D - maf) = ((z/beta)^2 / ((double)n + chi)/2)
//
//        double a = 2.0D * maf * (1.0D - maf) * ((double) n + chi);
//        double beta = z / Math.sqrt(a);
//        // Math.sqrt(a) = z / beta
//        // a = (z/beta)^2
//        double se = 1.0D / Math.sqrt(a);
//        return new double[]{beta, se};
//    }
}
