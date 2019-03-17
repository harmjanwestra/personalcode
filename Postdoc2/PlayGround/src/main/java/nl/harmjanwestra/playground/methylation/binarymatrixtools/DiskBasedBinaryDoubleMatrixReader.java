package nl.harmjanwestra.playground.methylation.binarymatrixtools;

import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.util.HashMap;

public class DiskBasedBinaryDoubleMatrixReader {

    private final File path;
    private final long headerLen;
    private DataInputStream is;
    private CountingInputStream counter;
    private final FileChannel channel;
    private final long filesize;

    public String[] rowIds;
    int nrRows;
    public String[] colIds;
    int nrCols;
    int bytesPerRow;
    HashMap<String, Integer> rowIndex;
    HashMap<String, Integer> colIndex;
    long currentPos;

    int buffersize = 32 * 1024;
    private ByteBuffer bytebuffer;

    public DiskBasedBinaryDoubleMatrixReader(File path) throws IOException {

        this.path = path;
        channel = new RandomAccessFile(path, "r").getChannel();
        filesize = channel.size();
        counter = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
        is = new DataInputStream(counter);

        rowIndex = new HashMap<>();
        nrRows = is.readInt();
        rowIds = new String[nrRows];
        for (int i = 0; i < nrRows; i++) {
            String rowId = is.readUTF();
            rowIndex.put(rowId, i);
            rowIds[i] = rowId;
        }

        colIndex = new HashMap<>();
        nrCols = is.readInt();
        colIds = new String[nrCols];
        for (int i = 0; i < nrCols; i++) {
            String colId = is.readUTF();
            colIndex.put(colId, i);
            colIds[i] = colId;
        }
        headerLen = counter.getByteCount();
        bytesPerRow = 8 * nrCols;

        currentPos = headerLen;
        System.out.println("Rows: " + nrRows + ", Cols: " + nrCols + ", bytes per row: " + bytesPerRow);
    }


    public double[] getNextRow(double[] buffer) throws IOException {
        for (int d = 0; d < nrCols; d++) {
            buffer[d] = is.readDouble();
        }
        currentPos += bytesPerRow;
        return buffer;
    }

    public double[] getNextRow() throws IOException {
        double[] output = new double[nrCols];
        return getNextRow(output);
    }


    public double get(int row, int col) throws IOException {
        return getRow(row)[col];
    }

    public double[] getRow(int row) throws IOException {
        return getRow(row, new double[nrCols]);
    }

    public double[] getRow(int row, double[] buffer) throws IOException {
        long seekLoc = ((long) row * bytesPerRow) + headerLen;

        if (seekLoc > filesize) {
            throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + filesize);
        }
        // if row is the next row, don't use random access..
        if (seekLoc - currentPos == 0) {

            return getNextRow(buffer);
        }

        // else use random access
        channel.position(seekLoc);
        if (bytebuffer == null) {
            bytebuffer = ByteBuffer.wrap(new byte[bytesPerRow]);

        }
        channel.read(bytebuffer);
        currentPos = seekLoc + bytesPerRow;

        // TODO: change the position of the datainputstream? this is probably very very slow...
        // question is whether this is needed when the position of the channel underlying a datainputstream is changed?
        // does DataInputStream have a counter of it's own?
        counter = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
        is = new DataInputStream(counter);
        return bytebuffer.asDoubleBuffer().array();

    }

    public void close() throws IOException {
        this.is.close();
        this.counter.close();
        this.channel.close();
    }

    public String[] getCols() {
        return colIds;
    }

    public int getNrCols() {
        return nrCols;
    }

    public int getNrRows() {
        return nrRows;
    }

    public String[] getRowIds() {
        return rowIds;
    }


}
