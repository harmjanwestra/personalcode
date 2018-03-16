/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package csvtott;

/**
 *
 * @author harmjan
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.regex.Pattern;
import umcg.genetica.containers.Triple;
import umcg.genetica.text.Strings;

/**
 *
 * @author harmjan
 */
public class ByteParseTask implements Callable<Triple<Integer, String, byte[]>> {

    String stringData;
    int offset = 0;
    int row = 0;
    private final int[] columnIndex;
    private final HashSet<String> hashRowsToInclude;
    private int nrColumnsToInclude;
    private Pattern separator;

    public ByteParseTask(String d, int offset, int row, int[] sampleIndex, HashSet<String> hashProbesToInclude) {
        this.stringData = d;
        this.offset = offset;
        this.row = row;
        this.columnIndex = sampleIndex; // converts column position in output to column position 
        this.hashRowsToInclude = hashProbesToInclude;
        if (columnIndex == null) {
            this.nrColumnsToInclude = -1;
        } else {
            this.nrColumnsToInclude = sampleIndex.length;
        }
    }
    
    public ByteParseTask(String d, int offset, int row, int[] sampleIndex, HashSet<String> hashProbesToInclude, Pattern separator) {
        this.stringData = d;
        this.offset = offset;
        this.row = row;
        this.columnIndex = sampleIndex; // converts column position in output to column position 
        this.hashRowsToInclude = hashProbesToInclude;
        if (columnIndex == null) {
            this.nrColumnsToInclude = -1;
        } else {
            this.nrColumnsToInclude = sampleIndex.length;
        }
        this.separator = separator;
    }

    @Override
    public Triple<Integer, String, byte[]> call() throws Exception {
        if (stringData == null) {
            return new Triple<Integer, String, byte[]>(-1, null, null);
        }

        if(separator == null){
            separator = Strings.tab;
        }
        
        String[] splitData = separator.split(stringData);
        String rowname = new String(splitData[0].getBytes());

        if (this.nrColumnsToInclude == -1) {
            this.nrColumnsToInclude = splitData.length - 1;
        }

        if (hashRowsToInclude == null || hashRowsToInclude.contains(rowname)) {
            byte[] output = new byte[nrColumnsToInclude];
            for (int s = 0; s < nrColumnsToInclude; s++) {
                int columnPositionInFile = s + offset;
                if (columnIndex != null) {
                    columnPositionInFile = columnIndex[s] + offset;
                }

                try {
                    double d = Double.parseDouble(splitData[columnPositionInFile]);
                    int dosageInt = (int) Math.round(d * 100d);
                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                    output[s] = dosageByte;
                } catch (NumberFormatException e) {
                    System.err.println("ERROR! Value is not a double: " + splitData[columnPositionInFile] + "\trow: " + row + "\tcol:" + columnPositionInFile);
                }
            }
            this.stringData = null;
            splitData = null;
            return new Triple<Integer, String, byte[]>(row, rowname, output);
        } else {
            this.stringData = null;
            splitData = null;
            return new Triple<Integer, String, byte[]>(-1, null, null);
        }

    }
}

