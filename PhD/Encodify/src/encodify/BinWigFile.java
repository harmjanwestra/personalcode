/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package encodify;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class BinWigFile extends BinaryFile {

    public BinWigFile(String name, boolean mode) throws IOException {
        super(name, mode);
    }

    public void convertWigFile(String file) throws IOException {
        TextFile wigFile = new TextFile(file, TextFile.R);
        // get the header
        // check wiggle type

        boolean fixedStep = false;
        int span = 1;
        byte chr = -1;
        int start = -1;
        int step = 1;
        int nextPos = start;
        int nrHeaders = 0;
        int nrElems = 0;


        HashMap<Integer, Integer> nrElemsPerHeader = new HashMap<Integer, Integer>();
        String[] elems = wigFile.readLineElems(Strings.whitespace);
        while (elems != null) {
            // check whether this line defines a fixed/variableStep region
            if (elems.length > 1) {
                for (int i = 0; i < elems.length; i++) {
                    if (elems[i].toLowerCase().contains("variablestep")) {
                        if (nrHeaders > 0) {
                            System.out.println(fixedStep + "\t" + nrHeaders + "\t" + nrElems);
                            nrElemsPerHeader.put(nrHeaders - 1, nrElems);
                            nrElems = 0;
                        }
                        fixedStep = false;
                        nrHeaders++;
                    } else if (elems[i].toLowerCase().contains("fixedstep")) {
                        if (nrHeaders > 0) {
                            System.out.println(fixedStep + "\t" + nrHeaders + "\t" + nrElems);
                            nrElemsPerHeader.put(nrHeaders - 1, nrElems);
                            nrElems = 0;
                        }
                        fixedStep = true;
                        nrHeaders++;
                    } else if (elems[i].toLowerCase().contains("chrom")) {
                        chr = ChrAnnotation.parseChr(elems[i].substring(9));
                    } else if (elems[i].toLowerCase().contains("span")) {
                        span = Integer.parseInt(elems[i].substring(5));
                    } else if (elems[i].toLowerCase().contains("start")) {
                        start = Integer.parseInt(elems[i].substring(6));
                    } else if (elems[i].toLowerCase().contains("step")) {
                        step = Integer.parseInt(elems[i].substring(5));
                    }
                }
            }

            String swtch = elems[0].toLowerCase();
            if (swtch.equals("variablestep")) {
                fixedStep = false;
            } else if (swtch.equals("fixedstep")) {
                fixedStep = true;
                nextPos = start;
            } else {
                nrElems++;
            }

            elems = wigFile.readLineElems(Strings.whitespace);
        }
        System.out.println(nrHeaders + "\t" + nrElems);
        nrElemsPerHeader.put(nrHeaders - 1, nrElems);
        wigFile.close();

        wigFile.open();
        elems = wigFile.readLineElems(Strings.whitespace);
        nrHeaders = 0;
        while (elems != null) {
            // check whether this line defines a fixed/variableStep region
            if (elems.length > 1) {
                boolean lineIsHeader = false;
                for (int i = 0; i < elems.length; i++) {
                    if (elems[i].toLowerCase().contains("variablestep")) {
                        lineIsHeader = true;
                        fixedStep = false;
                    } else if (elems[i].toLowerCase().contains("fixedstep")) {
                        lineIsHeader = true;
                        fixedStep = true;
                    } else if (elems[i].toLowerCase().contains("chrom")) {
                        chr = ChrAnnotation.parseChr(elems[i].substring(9));
                    } else if (elems[i].toLowerCase().contains("span")) {
                        span = Integer.parseInt(elems[i].substring(5));
                    } else if (elems[i].toLowerCase().contains("start")) {
                        start = Integer.parseInt(elems[i].substring(6));
                    } else if (elems[i].toLowerCase().contains("step")) {
                        step = Integer.parseInt(elems[i].substring(5));
                    }
                }
                if (lineIsHeader) {
                    nrElems = nrElemsPerHeader.get(nrHeaders);
                    os.writeBoolean(fixedStep);
                    os.writeByte(chr);
                    os.writeInt(start);
                    os.writeInt(step);
                    os.writeInt(span);
                    os.writeInt(nrElems);
                    nrHeaders++;
                }
            }

            String swtch = elems[0].toLowerCase();

            if (swtch.equals("variablestep")) {
                fixedStep = false;
            } else if (swtch.equals("fixedstep")) {
                fixedStep = true;
                nextPos = start;
            } else {
                if (fixedStep) {
                    double value = Double.parseDouble(elems[0]);
                    os.writeDouble(value);
                } else {
                    int chrPos = Integer.parseInt(elems[0]);
                    double value = Double.parseDouble(elems[1]);
                    os.writeInt(chrPos);
                    os.writeDouble(value);

                }
                nrElems++;
            }



            elems = wigFile.readLineElems(Strings.whitespace);
        }
        System.out.println(nrHeaders + "\t" + nrElems);
        wigFile.close();
    }
}
