package nl.harmjanwestra.playground.biogen.freeze2dot1;

public class MedianAgreementieqtl {
    public static void main(String[] args) {

        // cortex-AFR
        double[] pi = new double[]{0.00, 0.00, 0.06, 0.00, 0.00, 0.00, 0.00};
        double[] rb = new double[]{0.98, 0.74, 0.91, -0.17, 0.77, 0.81, 0.75};
        double[] ac = new double[]{0.72, 0.64, 0.72, 0.43, 0.67, 0.71, 0.65};

        double mp1 = JSci.maths.ArrayMath.median(pi);
        double mrb = JSci.maths.ArrayMath.median(rb);
        double mac = JSci.maths.ArrayMath.median(ac);

        System.out.println("Cortex AFR: mrb " + mrb + "\tmp1 " + mp1 + "\tmac " + mac);

        // ROSMAP
        ac = new double[]{0.72, 0.57, 0.81, 0.48, 0.55, 0.77};
        pi = new double[]{0.16, 0.34, 0.63, 0.00, 0.38, 0.48};
        rb = new double[]{0.81, 0.21, 0.82, 0.13, 0.54, 0.81};


        mp1 = JSci.maths.ArrayMath.median(pi);
        mrb = JSci.maths.ArrayMath.median(rb);
        mac = JSci.maths.ArrayMath.median(ac);

        System.out.println("ROSMAP: mrb " + mrb + "\tmp1 " + mp1 + "\tmac " + mac);

        // bryois
        ac = new double[]{0.90, 0.86, 0.90, 0.63, 0.81, 0.90};
        pi = new double[]{0.83, 0.43, 0.78, 0.71, 0.69, 0.65};
        rb = new double[]{0.86, 0.79, 0.84, 0.26, 0.78, 0.85};

        ac = new double[]{0.90, 0.86, 0.90, 0.81, 0.90};
        pi = new double[]{0.83, 0.43, 0.78, 0.69, 0.65};
        rb = new double[]{0.86, 0.79, 0.84, 0.78, 0.85};

        mp1 = JSci.maths.ArrayMath.median(pi);
        mrb = JSci.maths.ArrayMath.median(rb);
        mac = JSci.maths.ArrayMath.median(ac);

        System.out.println("Bryois: mrb " + mrb + "\tmp1 " + mp1 + "\tmac " + mac);
        System.out.println(JSci.maths.ArrayMath.min(rb) + "\t" + JSci.maths.ArrayMath.max(rb));
        System.out.println(JSci.maths.ArrayMath.min(pi) + "\t" + JSci.maths.ArrayMath.max(pi));
        System.out.println(JSci.maths.ArrayMath.min(ac) + "\t" + JSci.maths.ArrayMath.max(ac));
    }
}
