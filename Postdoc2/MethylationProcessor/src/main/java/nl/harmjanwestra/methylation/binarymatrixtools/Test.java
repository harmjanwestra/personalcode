package nl.harmjanwestra.methylation.binarymatrixtools;

import java.io.IOException;

public class Test {

    public static void main(String[] args){

        String file = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test.txt";
        String out = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test.dat";
        String out2 = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test2.txt";

        MatrixConverter cv = new MatrixConverter();
        TransposeMatrix t = new TransposeMatrix();

        try {
//            cv.toBinary(file,out);
//            cv.toText(out,out2);

            String outtranspose = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test-transpose.dat";
            t.transpose(out,outtranspose);
            String outtransposetext = "S:\\projects\\2018-methylation\\GPL13534_450k_OUT_Merged2\\test-transpose.txt";
            cv.toText(outtranspose,outtransposetext);
        } catch (IOException e) {
            e.printStackTrace();
        }


    }
}
