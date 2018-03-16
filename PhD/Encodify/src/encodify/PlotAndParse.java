/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package encodify;

import com.lowagie.text.DocumentException;
import com.lowagie.text.Rectangle;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;
import umcg.genetica.io.ucsc.PeakFile;
import umcg.genetica.io.ucsc.UCSCDataObject;
import umcg.genetica.io.ucsc.UCSCTrack;
import umcg.genetica.io.ucsc.WigFile;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class PlotAndParse {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            PlotAndParse p = new PlotAndParse();
//            String query = "rs4910714";
//            String proxyFile = "d:\\work\\proxies.txt"; // /Users/harmjan/proxies.txt
//            String snpAnnotationFile = "d:\\work\\SNPMappings.txt";
//            String encodeDir = "d:\\work\\wgEncodeBroadChipSeq\\";
//            String outDir = "d:\\work\\";

            String query = "rs653178";
            String proxyFile = "/Users/harmjan/proxies-rs653178.txt";// "/Users/harmjan/proxies.txt"; // 
            String snpAnnotationFile = "/Data/1000g/SNPMappings.txt";
            String encodeDir = "/Volumes/BackupDisk/Encode/hg18/";
            String outDir = "/Volumes/iSnackHD/EncodeOutput/";
            int margin = 50;

//            Gpio.createDir(outDir);
//            p.processFiles(query, proxyFile, snpAnnotationFile, encodeDir, margin, outDir);
//            
            query = "rs653178";
            proxyFile = "/Users/harmjan/proxies-rs653178.txt";// "/Users/harmjan/proxies.txt"; // 
            snpAnnotationFile = "/Data/1000g/SNPMappings.txt";
            encodeDir = "/Volumes/BackupDisk/Encode/hg18/";
            outDir = "/Volumes/iSnackHD/EncodeOutputPeakFiles/";
            Gpio.createDir(outDir);
            p.processPeakFiles(query, proxyFile, snpAnnotationFile, encodeDir, margin, outDir);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void processPeakFiles(String query, String proxyfile, String snpAnnotationFile, String encodeDir, int margin, String outDir) throws IOException {

        TextFile tf = new TextFile(proxyfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> querySNPs = new HashSet<String>();
        querySNPs.add(query);
        while (elems != null) {
            if (elems[0].equals(query)) {
                querySNPs.add(elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(snpAnnotationFile, TextFile.R);
        elems = tf2.readLineElems(TextFile.tab);
        ArrayList<SortableSNP> snps = new ArrayList<SortableSNP>();
        while (elems != null) {
            if (elems.length > 2 && querySNPs.contains(elems[2])) {
                System.out.println(Strings.concat(elems, Strings.tab));
                byte chr = ChrAnnotation.parseChr(elems[0]);
                int chrPos = Integer.parseInt(elems[1]);
                SortableSNP snp = new SortableSNP(elems[2], 0, chr, chrPos, SortableSNP.SORTBY.CHRANDCHRPOS);
                snps.add(snp);
                System.out.println("Added: " + snp);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf.close();

        if (snps.size() > 0) {
            Collections.sort(snps);


            byte chr = snps.get(0).chr;

            int chrPosStart = snps.get(0).chrpos;
            int chrPosEnd = snps.get(snps.size() - 1).chrpos;
            int size = (chrPosEnd) - (chrPosStart);
            int reqWidth = 500;
            int remainder = size % reqWidth;
            System.out.println(size + "\t" + remainder);
            int toAdd = 500 - remainder;
            if (toAdd % 2 > 0) {
                chrPosStart -= 1;
                toAdd--;
            }
            toAdd /= 2;
            chrPosStart -= toAdd;
            chrPosEnd += toAdd;
            System.out.println(toAdd);



            System.out.println("Range: " + chr + ": " + chrPosStart + " - " + chrPosEnd);

            size = (chrPosEnd) - (chrPosStart);
            System.out.println("Final range: " + chr + ": " + (chrPosStart - margin) + " - " + (chrPosEnd + margin) + " (" + size + "bp)");
            HashSet<String> filesInDir = getListOfFiles(encodeDir);

            int fctr = 0;
            ArrayList<UCSCTrack> tracks = new ArrayList<UCSCTrack>();
//            double min = Double.MAX_VALUE;
//            double max = Double.MIN_VALUE;
            String prevDir = null;
            ArrayList<String> allFiles = new ArrayList<String>();
            allFiles.addAll(filesInDir);
            Collections.sort(allFiles);
            System.out.println(allFiles.size() + " files to parse");
            int fileCtr = 0;
            for (String f : allFiles) {

                if ((f.contains("broadPeak") || f.contains("gappedPeak")) || f.contains("narrowPeak")) {
                    String origFileName = f;
                    double min = Double.MAX_VALUE;
                    double max = Double.MIN_VALUE;

                    String[] filenameElems = f.split("/");
                    if (filenameElems.length == 1) {
                        filenameElems = f.split("\\\\");
                    }
//                    System.out.println(Strings.concat(filenameElems, Strings.tab));
                    f = f.replace(encodeDir, "");
                    if (filenameElems.length > 3 && f.contains(filenameElems[filenameElems.length - 2])) {
                        f = f.replace(filenameElems[filenameElems.length - 2], "");
                        if (prevDir == null) {
                            prevDir = filenameElems[filenameElems.length - 2];
                            System.out.println("Iterating dir: " + prevDir);
                        } else {
                            // start new plot
                            if (!prevDir.equals(filenameElems[filenameElems.length - 2])) {
                                // new dir. 

                                // plot stuff
                                Collections.sort(tracks);
                                try {
                                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-" + prevDir + ".pdf");
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }

                                // new databin
                                prevDir = filenameElems[filenameElems.length - 2];
                                System.out.println("Iterating dir: " + prevDir);
                                tracks = new ArrayList<UCSCTrack>();
                            }
                        }
                    }



                    String trackName = f.replace(".wig.gz", "");
                    UCSCTrack t = new UCSCTrack(trackName);


                    PeakFile.PEAKFORMAT format = null;
                    if (f.toLowerCase().contains("broadpeak")) {
                        format = PeakFile.PEAKFORMAT.BROADPEAK;
                    } else if (f.toLowerCase().contains("narrowpeak")) {
                        format = PeakFile.PEAKFORMAT.NARROWPEAK;
                    } else if (f.toLowerCase().contains("gappedpeak")) {
                        format = PeakFile.PEAKFORMAT.GAPPEDPEAK;
                    }
                    PeakFile wf = new PeakFile(origFileName, PeakFile.R, format);
                    System.out.println("Parsing file: " + f + " - " + Encodify.convertSize(wf.size()));

                    ArrayList<Double> vals = new ArrayList<Double>();
                    UCSCDataObject d = wf.parseLn();
                    while (d != null) {
//                        System.out.println(d.toString());
                        if (d.getChr() == chr && d.getPositionEnd() >= chrPosStart && d.getPositionStart() <= chrPosEnd) {
//                            System.out.println(d.toString());
                            t.addDataObject(d);
                        }
                        if (d.getChr() == chr) {
                            vals.add(d.getValue());
                            if (d.getValue() > max) {
                                max = d.getValue();

                            } else if (d.getValue() < min) {
                                min = d.getValue();
                            }
                        }
                        d = wf.parseLn();
                    }
                    if (vals.size() > 0) {
                        double[] vals2 = new double[vals.size()];
                        for (int v = 0; v < vals.size(); v++) {
                            vals2[v] = vals.get(v);
                        }
                        double mean = JSci.maths.ArrayMath.mean(vals2);
                        double sd = JSci.maths.ArrayMath.standardDeviation(vals2);
                        System.out.println(fileCtr + "/" + allFiles.size() + " - " + t.getName() + ":\t" + t.getData().size() + " objects loaded for region. " + vals.size() + " total for chr.  Chr Min: " + min + " || Max: " + max + " || Avg: " + mean + "\t" + sd);
                        wf.close();
                        if (t.getData().size() > 0) {
                            tracks.add(t);
                            t.setMax(max);
                            t.setMin(min);
                            t.setSD(sd);
                            t.setMean(mean);
                            fctr++;
//                        if (fctr == 2) {
//                            break;
//                        }
                        }
                    } else {
                        System.out.println("No values for chr!");
                    }

                }
                fileCtr++;
            }

            // loaded tracks. now plot the fuckers..
            // get max filename:
            // nrOf Tracks: 

            try {
                Collections.sort(tracks);
                if (prevDir != null) {
                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-" + prevDir + ".pdf");
                } else {
                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-output.pdf");
                }

            } catch (Exception e) {
                e.printStackTrace();
            }


        }

    }

    public void processFiles(String query, String proxyfile, String snpAnnotationFile, String encodeDir, int margin, String outDir) throws IOException {

        TextFile tf = new TextFile(proxyfile, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> querySNPs = new HashSet<String>();
        querySNPs.add(query);
        while (elems != null) {
            if (elems[0].equals(query)) {
                querySNPs.add(elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile tf2 = new TextFile(snpAnnotationFile, TextFile.R);
        elems = tf2.readLineElems(TextFile.tab);
        ArrayList<SortableSNP> snps = new ArrayList<SortableSNP>();
        while (elems != null) {
            if (elems.length > 2 && querySNPs.contains(elems[2])) {
                System.out.println(Strings.concat(elems, Strings.tab));
                byte chr = ChrAnnotation.parseChr(elems[0]);
                int chrPos = Integer.parseInt(elems[1]);
                SortableSNP snp = new SortableSNP(elems[2], 0, chr, chrPos, SortableSNP.SORTBY.CHRANDCHRPOS);
                snps.add(snp);
                System.out.println("Added: " + snp);
            }
            elems = tf2.readLineElems(TextFile.tab);
        }
        tf.close();

        if (snps.size() > 0) {
            Collections.sort(snps);


            byte chr = snps.get(0).chr;

            int chrPosStart = snps.get(0).chrpos;
            int chrPosEnd = snps.get(snps.size() - 1).chrpos;
            int size = (chrPosEnd) - (chrPosStart);
            int reqWidth = 500;
            int remainder = size % reqWidth;
            System.out.println(size + "\t" + remainder);
            int toAdd = 500 - remainder;
            if (toAdd % 2 > 0) {
                chrPosStart -= 1;
                toAdd--;
            }
            toAdd /= 2;
            chrPosStart -= toAdd;
            chrPosEnd += toAdd;
            System.out.println(toAdd);



            System.out.println("Range: " + chr + ": " + chrPosStart + " - " + chrPosEnd);

            size = (chrPosEnd) - (chrPosStart);
            System.out.println("Final range: " + chr + ": " + (chrPosStart - margin) + " - " + (chrPosEnd + margin) + " (" + size + "bp)");
            HashSet<String> filesInDir = getListOfFiles(encodeDir);

            int fctr = 0;
            ArrayList<UCSCTrack> tracks = new ArrayList<UCSCTrack>();
//            double min = Double.MAX_VALUE;
//            double max = Double.MIN_VALUE;
            String prevDir = null;
            ArrayList<String> allFiles = new ArrayList<String>();
            allFiles.addAll(filesInDir);
            Collections.sort(allFiles);
            System.out.println(allFiles.size() + " files to parse");
            int fileCtr = 0;
            for (String f : allFiles) {

                if (f.endsWith(".wig.gz") && f.contains("wgEncodeBroadChipSeq") && !f.contains("ChromatinMap")) {
                    String origFileName = f;
                    double min = Double.MAX_VALUE;
                    double max = Double.MIN_VALUE;

                    String[] filenameElems = f.split("/");
                    if (filenameElems.length == 1) {
                        filenameElems = f.split("\\\\");
                    }
//                    System.out.println(Strings.concat(filenameElems, Strings.tab));
                    f = f.replace(encodeDir, "");
                    if (filenameElems.length > 3 && f.contains(filenameElems[filenameElems.length - 2])) {
                        f = f.replace(filenameElems[filenameElems.length - 2], "");
                        if (prevDir == null) {
                            prevDir = filenameElems[filenameElems.length - 2];
                        } else {
                            // start new plot
                            if (!prevDir.equals(filenameElems[filenameElems.length - 2])) {
                                // new dir. 

                                // plot stuff
                                Collections.sort(tracks);
                                try {
                                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-" + prevDir + ".pdf");
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }

                                // new databin
                                prevDir = filenameElems[filenameElems.length - 2];
                                tracks = new ArrayList<UCSCTrack>();
                            }
                        }
                    }



                    String trackName = f.replace(".wig.gz", "");
                    UCSCTrack t = new UCSCTrack(trackName);


                    WigFile wf = new WigFile(origFileName, WigFile.R);
                    System.out.println("Parsing file: " + f + " - " + Encodify.convertSize(wf.size()));

                    ArrayList<Double> vals = new ArrayList<Double>();
                    UCSCDataObject d = wf.parseLn();
                    while (d != null) {
                        if (d.getChr() == chr && d.getPositionEnd() >= chrPosStart && d.getPositionStart() <= chrPosEnd) {
//                            System.out.println(d.toString());
                            t.addDataObject(d);
                        }
                        if (d.getChr() == chr) {
                            vals.add(d.getValue());
                            if (d.getValue() > max) {
                                max = d.getValue();

                            } else if (d.getValue() < min) {
                                min = d.getValue();
                            }
                        }
                        d = wf.parseLn();
                    }
                    double[] vals2 = new double[vals.size()];
                    for (int v = 0; v < vals.size(); v++) {
                        vals2[v] = vals.get(v);
                    }
                    double mean = JSci.maths.ArrayMath.mean(vals2);
                    double sd = JSci.maths.ArrayMath.standardDeviation(vals2);
                    System.out.println(fileCtr + "/" + allFiles.size() + " - " + t.getName() + ":\t" + t.getData().size() + " objects loaded for region. " + vals.size() + " total for chr.  Chr Min: " + min + " || Max: " + max + " || Avg: " + mean + "\t" + sd);
                    wf.close();
                    if (t.getData().size() > 0) {
                        tracks.add(t);
                        t.setMax(max);
                        t.setMin(min);
                        t.setSD(sd);
                        t.setMean(mean);
                        fctr++;
//                        if (fctr == 2) {
//                            break;
//                        }
                    }
                }
                fileCtr++;
            }

            // loaded tracks. now plot the fuckers..
            // get max filename:
            // nrOf Tracks: 

            try {
                Collections.sort(tracks);
                if (prevDir != null) {
                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-" + prevDir + ".pdf");
                } else {
                    plot(reqWidth, 5, chrPosStart, chrPosEnd, snps, tracks, outDir + query + "-output.pdf");
                }

            } catch (Exception e) {
                e.printStackTrace();
            }


        }

    }
    private static final Font COL_FONT = new Font("Verdana", Font.PLAIN, 6);

    public void plot(int requestedWidth, int boxWidth, int start, int end, ArrayList<SortableSNP> snps, ArrayList<UCSCTrack> tracks, String filename) throws DocumentException, FileNotFoundException, IOException {
        System.out.println("Plotting: " + filename);
        String longestString = "";
        int longestStrSize = 0;
        for (UCSCTrack t : tracks) {
            if (t.getName().length() > longestStrSize) {
                longestString = t.getName();
                longestStrSize = t.getName().length();
            }
        }



        int topMargin = 50;
        int plotMargin = 10; // 10 px margin;
        int boxSpacer = 5;
        int boxHeight = boxWidth;
        int plotHeight = topMargin + plotMargin * 2 + tracks.size() * boxWidth;
        int plotWidth = plotMargin * 2 + requestedWidth;

        System.out.println("Plot will be: " + plotHeight + " x " + plotWidth);
        BufferedImage bi = new BufferedImage(plotWidth, plotHeight, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = bi.createGraphics();



        FontMetrics fontmetrics = g2d.getFontMetrics();
        int leftMargin = 0;
        int rightSpacer = fontmetrics.stringWidth(longestString);

        plotWidth = plotWidth + boxSpacer + rightSpacer + plotMargin;
        System.out.println("Plot will be: " + plotHeight + " x " + plotWidth);
        Rectangle rectangle = new Rectangle(plotWidth, plotHeight);
        com.lowagie.text.Document document = new com.lowagie.text.Document(rectangle);
        com.lowagie.text.pdf.PdfWriter writer = com.lowagie.text.pdf.PdfWriter.getInstance(document, new java.io.FileOutputStream(filename));
//        bi = new BufferedImage(plotWidth, plotHeight, BufferedImage.TYPE_INT_RGB);
//        g2d = bi.createGraphics();
        document.open();
        com.lowagie.text.pdf.PdfContentByte cb = writer.getDirectContent();
        com.lowagie.text.pdf.DefaultFontMapper fontMap = new com.lowagie.text.pdf.DefaultFontMapper();
        g2d = cb.createGraphics(plotWidth, plotHeight);

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, plotWidth, plotHeight);



        // draw the tracks
        int range = end - start;

        int nrBins = requestedWidth / boxWidth;
        double bpPerBin = (double) range / nrBins; // bp / bin

        System.out.println("Nr Bins: " + nrBins);
        System.out.println("BP per bin: " + bpPerBin);
        g2d.setFont(COL_FONT);

        // draw a line
        g2d.setColor(Color.gray);

        int lnXStart = plotMargin;
        int lnXEnd = plotMargin + requestedWidth;
        int lnY = topMargin;

        g2d.drawLine(lnXStart, lnY - 2, lnXStart, lnY + 2);
        g2d.drawLine(lnXEnd, lnY - 2, lnXEnd, lnY + 2);
        g2d.drawLine(lnXStart, lnY, lnXEnd, lnY);

        g2d.setColor(Color.black);

        // draw each SNP:
        for (SortableSNP s : snps) {
            int pos = s.chrpos - start;
            int posRange = end - start;
            double posPErc = (double) pos / posRange;
            posPErc *= requestedWidth;
            int posX = plotMargin + (int) Math.floor(posPErc);

//            System.out.println(s.chrpos + "\t" + start + "\t" + pos + "\t" + posX);


            g2d.drawLine(posX, lnY - 2, posX, lnY + 2);
            AffineTransform orig = g2d.getTransform();

            double angle = Math.toRadians(45);
            g2d.rotate(-angle, posX, lnY - 2 - 5);
            g2d.drawString(s.name, posX, lnY - 2 - 5);
            g2d.setTransform(orig);

        }




        for (int i = 0; i < tracks.size(); i++) {
            double[] sumPerBin = new double[nrBins];
            int[] nrValsInBin = new int[nrBins];
            int yPos = topMargin + plotMargin + (boxHeight * i) + boxHeight;
            // sort the values on the tracks
            UCSCTrack t = tracks.get(i);
            t.sort();

            // now get the values

            // correct values for start and end
            double maxZ = Double.MIN_VALUE;
            double minZ = Double.MAX_VALUE;
            for (UCSCDataObject d : t.getData()) {
                int sta = d.getPositionStart() - start;
                int sto = d.getPositionEnd() - start;



                // get the bin number
                double midpoint = (double) sto + sta / 2;
                int binNo = (int) (Math.floor(midpoint / bpPerBin));
                if (binNo >= nrBins) {
                    binNo = nrBins - 1;
                }

                nrValsInBin[binNo]++;
                double z = (d.getValue() - t.getMean() / t.getSD());
                sumPerBin[binNo] += z;
                if (z > maxZ) {
                    maxZ = z;
                } else if (z < minZ) {
                    minZ = z;
                }
//                System.out.println(d.positionStart + " - " + sta + "\t" + d.positionEnd + " - " + sto + "\t" + d.value + "\t" + d.value + " - z: " + (d.value - t.mean / t.sd) + "\t" + midpoint + "\t" + binNo);
            }

//            System.out.println("Max: " + maxZ + " || min: " + minZ);

            for (int v = 0; v < sumPerBin.length; v++) {
                double avg = sumPerBin[v] / nrValsInBin[v];
                Color c = getColor(minZ, 15, avg);
                g2d.setColor(c);
                int xPos = plotMargin + v * boxHeight;
                g2d.fillRect(xPos, yPos - boxHeight, boxHeight, boxHeight);
            }

            // start now: < x > range



            // plot name

            int xPos = plotMargin + sumPerBin.length * boxHeight + boxSpacer;
            DecimalFormat df = new DecimalFormat("#.##");
            String maxZStr = "" + df.format(maxZ);
            g2d.setColor(Color.black);
            g2d.drawString(t.getName() + " (" + maxZStr + ")", xPos, yPos);


        }

        g2d.dispose();
        bi.flush();
//        ImageIO.write(bi, "png", new File(filename));
        document.close();
        writer.close();
    }

    public HashSet<String> getListOfFiles(String dir) {
        System.out.println("Iterating: " + dir);
        HashSet<String> fileList = new HashSet<String>();
        String[] files = Gpio.getListOfFiles(dir);
        for (String file : files) {
            String loc = dir + file;
            if (Gpio.isDir(loc)) {
                fileList.addAll(getListOfFiles(loc + "/"));
            } else {
                fileList.add(loc);
            }
        }
        return fileList;
    }

    public Color getColor(double min, double max, double v) {
        int r = 0;
        int g = 0;
        int b = 0;
        double perc = 0;
        if (min < 0) {
            min = Math.abs(min);
            max += min;
            perc = v + min;
        } else {
            max -= min;
            perc = v - min;
        }

        perc /= max;
        int aPerc = (int) Math.floor(255 * perc); //  
        if (aPerc < 0) {
            aPerc = 0;
        } else if (aPerc > 255) {
            aPerc = 255;
        }
//        System.out.println(v + "\t" + aPerc);
        return new Color(r, g, b, aPerc);
    }
}
