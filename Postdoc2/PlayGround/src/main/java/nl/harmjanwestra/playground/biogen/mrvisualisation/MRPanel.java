package nl.harmjanwestra.playground.biogen.mrvisualisation;

import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.Panel;
import umcg.genetica.graphics.themes.DefaultTheme;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class MRPanel extends Panel {

    private final int marginBetweenDots;
    private final int marginBetweenBoxes;
    private final int marginWithinBoxes;
    private final int dotsize;
    private final int halfdotsize;
    private final int boxSize;
    private HashMap<String, ArrayList<MREvent>> data;
    private Color[] greyScale;
    private Color[] colors;
    private DefaultTheme deftheme;
    private ORIENTATION orientation = ORIENTATION.VERTICAL;
    private Range rangeBeta;
    private Range rangeWR;
    private boolean determineRangeOverAllEvents;
    private Graphics2D g2d;
    private double pp4thresh = 0.7;
    private boolean overlapBoxes;

    private enum ORIENTATION {
        HORIZONTAL,
        VERTICAL
    }

    public void determineRangeOverAllEvents() {
        this.determineRangeOverAllEvents = true;
    }

    public void overlapBoxes() {
        this.overlapBoxes = true;
    }

    public MRPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);

        marginBetweenDots = 5;
        marginBetweenBoxes = 20;
        marginWithinBoxes = 10;
        dotsize = 5;
        halfdotsize = dotsize / 2;
        boxSize = 150;

        // define colors
        // metabrain color palette
        colors = new Color[]{
                new Color(0, 0, 0),
                new Color(230, 159, 0, 75),
                new Color(86, 180, 233, 75),
                new Color(204, 121, 167),
                new Color(240, 228, 66),
                new Color(0, 114, 178),
                new Color(213, 94, 0),
                new Color(0, 158, 115)
        };
        greyScale = new Color[]{
                new Color(0, 0, 0),
                new Color(128, 128, 128),
                new Color(169, 169, 169),
                new Color(211, 211, 211),
                new Color(232, 232, 232)
        };
    }

    public void setOrientation(ORIENTATION orientation) {
        this.orientation = orientation;
    }

    public void setData(ArrayList<MREvent> data) {
        HashMap<String, ArrayList<MREvent>> tmp = new HashMap<>();
        tmp.put(data.get(0).outcome, data);
        this.data = tmp;
    }

    @Override
    public void draw(DefaultGraphics defaultGraphics) {
        g2d = defaultGraphics.getG2d();
        g2d.setColor(greyScale[0]);

        deftheme = (DefaultTheme) theme;
        deftheme.setColors(colors);

        ArrayList<String> keys = new ArrayList<>();
        keys.addAll(data.keySet());
        Collections.sort(keys);

        if (determineRangeOverAllEvents) {
            ArrayList<MREvent> tmp = new ArrayList<>();
            for (int i = 0; i < keys.size(); i++) {
                String outcome = keys.get(i);
                ArrayList<MREvent> events = data.get(outcome);
                tmp.addAll(events);
            }
            determineRanges(tmp);
        }

        //g2d.setFont(deftheme.getSmallFontBold());
        //            g2d.drawString("Gene", x0 - getStringWidth(g2d, "Gene") - marginBetweenDots, y1 - 3);

        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[1]);
            g2d.setFont(deftheme.getSmallFont());
            int strw = getStringWidth(g2d, "Beta / Wald ratio");
            g2d.drawString("Beta / Wald ratio", x0 + (boxSize / 2) - strw, y0 - 20);
        }

        int startX = x0;
        int startY = y0;
        int octr = 0;
        for (int i = 0; i < keys.size(); i++) {
            String outcome = keys.get(i);
            ArrayList<MREvent> events = data.get(outcome);
            if (!determineRangeOverAllEvents) {
                rangeBeta = null;
                rangeWR = null;
                determineRanges(events);
            }

            System.out.println("Range WR: " + rangeWR);
            System.out.println("Range Beta: " + rangeBeta);

            plotData(events, startX, startY, outcome, defaultGraphics);
            if (orientation == ORIENTATION.HORIZONTAL) {
                int boxWidth = events.size() * (dotsize + marginBetweenDots);
                startX += boxWidth + marginBetweenBoxes;
            } else {
                int boxHeight = events.size() * (dotsize + marginBetweenDots) + (2 * marginWithinBoxes);
                startY += boxHeight + marginBetweenBoxes;
            }
            octr++;

        }
    }

    private void drawBox(int x1, int y1, int boxHeight, int boxWidth) {

        // draw 0 line
        g2d.setColor(greyScale[2]);
        Stroke defStroke = g2d.getStroke();
        g2d.setStroke(deftheme.strokeDashed);
        if (orientation == ORIENTATION.HORIZONTAL) {
            g2d.drawLine(x1, y1 + getPos(0, boxHeight, rangeBeta), x1 + boxWidth, y1 + getPos(0, boxHeight, rangeBeta));
        } else {
            g2d.drawLine(x1 + getPos(0, boxWidth, rangeBeta), y1, x1 + getPos(0, boxWidth, rangeBeta), y1 + boxHeight);
        }
        g2d.setStroke(defStroke);

        g2d.setColor(greyScale[0]);
        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[2]);
        }
        g2d.drawLine(x1, y1, x1, y1 + boxHeight); // |--
        g2d.drawLine(x1 + boxWidth, y1, x1 + boxWidth, y1 + boxHeight); // --|

        g2d.setColor(greyScale[4]);
        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[0]);
        }
        g2d.drawLine(x1, y1, x1 + boxWidth, y1); // --
        g2d.drawLine(x1, y1 + boxHeight, x1 + boxWidth, y1 + boxHeight); // __


    }

    private void plotData(ArrayList<MREvent> events, int x1, int y1, String outcome, DefaultGraphics defaultGraphics) {
        int boxWidth = 0;
        int boxHeight = 0;
        if (orientation == ORIENTATION.HORIZONTAL) {
            boxHeight = boxSize;
            boxWidth = events.size() * (dotsize + marginBetweenDots);
            drawBox(x1, y1, boxHeight, boxWidth); // draw beta box

            int y2 = y1 + boxHeight + marginBetweenBoxes;
            drawBox(x1, y2, boxHeight, boxWidth); // draw WR box


            // plot beta values
            int ectr = 0;
            for (MREvent event : events) {
                g2d.setColor(colors[1]);
                // plot Beta
                int posX = ectr * (dotsize + marginBetweenDots);
                g2d.fillOval(x1 + posX, y1 + getPos(event.betaExposure, boxHeight, rangeBeta), dotsize, dotsize);

                // plot outcome
                g2d.setColor(colors[2]);
                g2d.fillOval(x1 + posX, y1 + getPos(event.betaOutcome, boxHeight, rangeBeta), dotsize, dotsize);

                // plot WR
                if (event.pp4bon > pp4thresh || event.pp4def > pp4thresh) {
                    g2d.setColor(colors[3]);
                } else {
                    g2d.setColor(greyScale[2]);
                }
                g2d.fillOval(x1 + posX, y2 + getPos(event.waldratio, boxHeight, rangeWR), dotsize, dotsize);

                // event snp
                // event gene
                // event allele
            }
        } else {
            boxWidth = boxSize;
            boxHeight = events.size() * (dotsize + marginBetweenDots) + (2 * marginWithinBoxes);

            // draw shading every 5 lines
            if (events.size() > 10) {
                boolean shade = true;
                for (int i = 0; i < events.size(); i++) {
                    if (i > 0 && i % 5 == 0) {
                        if (shade) {
                            // shade the next 5 genes
                            int stop = i + 5;
                            if (stop > events.size()) {
                                stop = events.size();
                            }
                            int posY = y1 + i * (dotsize + marginBetweenDots) + (dotsize / 2) + marginWithinBoxes;
                            int posY2 = y1 + stop * (dotsize + marginBetweenDots) + (dotsize / 2) + marginWithinBoxes;
                            g2d.setColor(greyScale[4]);
                            g2d.fillRect(0, posY, width, posY2 - posY);
                            shade = false;
                        } else {
                            shade = true;
                        }
                    }
                }
            }

            // plot the scales
            g2d.setColor(greyScale[1]);
            DecimalFormat df = new DecimalFormat("#.##");
            g2d.drawString(df.format(rangeBeta.getMinX()), x1 - (getStringWidth(g2d, df.format(rangeBeta.getMinX())) / 2), y1 - 5);
            g2d.drawString(df.format(0), x1 + getPos(0, boxWidth, rangeBeta) - (getStringWidth(g2d, "0") / 2), y1 - 5);
            g2d.drawString(df.format(rangeBeta.getMaxX()), x1 + boxWidth - (getStringWidth(g2d, df.format(rangeBeta.getMaxX())) / 2), y1 - 5);


            drawBox(x1, y1, boxHeight, boxWidth); // draw beta box
            int x2 = x1;
            if (!overlapBoxes) {
                x2 = x1 + boxWidth + marginBetweenBoxes;
                drawBox(x1 + boxWidth + marginBetweenBoxes, y1, boxHeight, boxWidth); // draw WR box
            }
            int x3 = x2 + boxWidth + marginBetweenDots;


            g2d.setFont(deftheme.getSmallFont());
            int smallFontHeight = g2d.getFontMetrics().getHeight();

            int maxSNPStrWidth = 0;
            int maxGeneStrWidth = 0;
            for (MREvent event : events) {
                int strWidth = super.getStringWidth(g2d, event.EA + " - " + event.snp);
                if (strWidth > maxSNPStrWidth) {
                    maxSNPStrWidth = strWidth;
                }
                int strWidthGene = super.getStringWidth(g2d, event.gene);
                if (strWidthGene > maxGeneStrWidth) {
                    maxGeneStrWidth = strWidthGene;
                }
            }

            // plot outcome name
            g2d.setColor(greyScale[0]);
            g2d.setFont(deftheme.getSmallFontBold());
            g2d.drawString(outcome, x3, y1 - 3 + marginWithinBoxes);


            g2d.setFont(deftheme.getSmallFont());
            // plot beta values
            int ectr = 0;


            for (MREvent event : events) {

                g2d.setColor(colors[1]);
                // plot Beta exposure
                int posY = ectr * (dotsize + marginBetweenDots) + (dotsize / 2) + marginWithinBoxes;
                g2d.fillOval(x1 + getPos(event.betaExposure, boxWidth, rangeBeta) - (dotsize / 2), y1 + posY, dotsize, dotsize);

                // plot beta exposure confidence
                double conf1 = event.betaExposure - 1.96 * event.seExposure;
                double conf2 = event.betaExposure + 1.96 * event.seExposure;
                g2d.drawLine(x1 + getPos(conf1, boxWidth, rangeBeta), y1 + posY + (dotsize / 2), x1 + getPos(conf2, boxWidth, rangeBeta), y1 + posY + (dotsize / 2));


                // plot outcome
                g2d.setColor(colors[2]);
                g2d.fillOval(x1 + getPos(event.betaOutcome, boxWidth, rangeBeta) - (dotsize / 2), y1 + posY, dotsize, dotsize);

                conf1 = event.betaOutcome - 1.96 * event.seOutcome;
                conf2 = event.betaOutcome + 1.96 * event.seOutcome;
                g2d.drawLine(x1 + getPos(conf1, boxWidth, rangeBeta), y1 + posY + (dotsize / 2), x1 + getPos(conf2, boxWidth, rangeBeta), y1 + posY + (dotsize / 2));


                // plot WR
                if (event.pp4bon > pp4thresh || event.pp4def > pp4thresh) {
                    g2d.setColor(colors[3]);
                } else {
                    g2d.setColor(greyScale[2]);
                }
                g2d.fillOval(x2 + getPos(event.waldratio, boxWidth, rangeWR) - (dotsize / 2), y1 + posY, dotsize, dotsize);

                conf1 = event.waldratio - 1.96 * event.seWaldratio;
                conf2 = event.waldratio + 1.96 * event.seWaldratio;
                g2d.drawLine(x1 + getPos(conf1, boxWidth, rangeBeta), y1 + posY + (dotsize / 2), x1 + getPos(conf2, boxWidth, rangeBeta), y1 + posY + (dotsize / 2));


                // event snp
                g2d.setColor(greyScale[1]);
                g2d.setFont(deftheme.getSmallFont());

                posY = ectr * (dotsize + marginBetweenDots) + dotsize + marginBetweenDots + marginWithinBoxes; // + (smallFontHeight);
                g2d.drawString(event.snp + " - " + event.EA, x3, y1 + posY);
                // event gene
                g2d.drawString(event.gene, x0 - getStringWidth(g2d, event.gene) - marginBetweenDots, y1 + posY);

                ectr++;

            }
        }
    }

    private int getPos(double val, int size, Range range) {
        double percB = range.getRelativePositionX(val);
        int pixels = (int) Math.ceil(percB * size);
        return pixels;
    }

    private void determineRanges(ArrayList<MREvent> events) {
        // determine min, max
        double minWR = Double.MAX_VALUE;
        double minBeta = Double.MAX_VALUE;
        double maxWR = -Double.MAX_VALUE;
        double maxBeta = -Double.MAX_VALUE;
        for (MREvent e : events) {
            if (e.waldratio < minWR) {
                minWR = e.waldratio;
            }
            if (e.betaOutcome < minBeta) {
                minBeta = e.betaOutcome;
            }
            if (e.betaExposure < minBeta) {
                minBeta = e.betaExposure;
            }

            if (e.waldratio > maxWR) {
                maxWR = e.waldratio;
            }
            if (e.betaOutcome > maxBeta) {
                maxBeta = e.betaOutcome;
            }
            if (e.betaExposure > maxBeta) {
                maxBeta = e.betaExposure;
            }
        }

        // make symmetrical
        minBeta = -Math.max(Math.abs(minBeta), Math.abs(maxBeta));
        maxBeta = -minBeta;
        minWR = -Math.max(Math.abs(minWR), Math.abs(maxWR));
        maxWR = -minWR;

        // determine ranges of values and round
        rangeBeta = new Range(minBeta, minBeta, maxBeta, maxBeta);
        rangeBeta.roundX();
        rangeBeta.roundY();
        rangeWR = new Range(minWR, minWR, maxWR, maxWR);
        rangeWR.roundX();
        rangeWR.roundY();

    }

    public void setData(HashMap<String, ArrayList<MREvent>> mrEvents) {
        this.data = mrEvents;
    }
}
