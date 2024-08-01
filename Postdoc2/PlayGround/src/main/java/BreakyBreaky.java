import umcg.genetica.math.stats.ZScores;

public class BreakyBreaky {


    public static void main(String[] args) {
        int nextChr = 1;
        boolean vcfFileExists = true;
        String origfile = "chr";
        String newfile = origfile + "-chr" + nextChr;
        while (nextChr < 25 && vcfFileExists) {

            newfile = origfile + "-chr" + nextChr;
            if (nextChr == 2) {
                vcfFileExists = false;
            }

            if (nextChr == 21) {
                vcfFileExists = false;
            }
            System.out.println(nextChr + "\t" + vcfFileExists);
            nextChr += 1;
        }

        System.out.println("pick:" + nextChr + "\t" + vcfFileExists + "\t" + newfile);
    }
}
