package nl.harmjanwestra.playground.methylation.binarymatrixtools;

import org.apache.commons.io.input.CountingInputStream;
import org.apache.commons.io.output.CountingOutputStream;
import umcg.genetica.io.bin.BinaryFile;

import java.io.*;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;

public class DiskBasedBinaryDoubleMatrixWriter {


    FileChannel channel;
    private int nrCols;
    private int nrRows;
    private int buffersize = 32 * 1024;
    private CountingOutputStream counter;
    private DataOutputStream os;
    private int bytesPerRow;
    private long headerLen;
    private long currentPos;
    private ByteBuffer bytebuffer;


    public void initializeFullMatrix(String[] rows, String[] cols, String out) throws IOException {
        initialize(rows, cols, out);
        initializeWithNaNs();
    }

    public void open(String loc) throws IOException {
        channel = new RandomAccessFile(loc, "rw").getChannel();
        CountingInputStream countertmp = new CountingInputStream(new BufferedInputStream(Channels.newInputStream(channel), buffersize));
        DataInputStream is = new DataInputStream(countertmp);
        nrRows = is.readInt();
        for (int i = 0; i < nrRows; i++) {
            is.readUTF();
        }
        nrCols = is.readInt();
        for (int i = 0; i < nrCols; i++) {
            is.readUTF();
        }
        headerLen = countertmp.getByteCount();
        currentPos = headerLen;
        bytesPerRow = 8 * nrCols;

        channel.position(currentPos);
        counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
        os = new DataOutputStream(counter);

        System.out.println("Read header. current pos: " + channel.position());
        System.out.println("Header: " + headerLen);
    }


    public void initialize(String[] rows, String[] cols, String out) throws IOException {
        channel = new RandomAccessFile(out, "rw").getChannel();
        counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
        os = new DataOutputStream(counter);

        writeHeader(rows);
        writeHeader(cols);
        headerLen = counter.getByteCount();
        currentPos = counter.getByteCount();
        this.nrCols = cols.length;
        this.nrRows = rows.length;
        bytesPerRow = cols.length * 8;
    }

    private void writeHeader(String[] rows) throws IOException {
        os.writeInt(rows.length);
        for (String s : rows) {
            os.writeUTF(s);
        }
    }

    private void initializeWithNaNs() throws IOException {
        for (int row = 0; row < nrRows; row++) {
            for (int col = 0; col < nrCols; col++) {
                os.writeDouble(Double.NaN);
            }
        }
    }


    public void writeRow(double[] cols) throws IOException {
        if (cols.length != nrCols) {
            throw new IllegalArgumentException("Length of cols not equal: " + cols.length + " found, " + nrCols + " expected.");
        }
        for (int c = 0; c < cols.length; c++) {
            os.writeDouble(cols[c]);
        }
        currentPos += bytesPerRow;
    }

    public void write(double d) throws IOException {
        os.writeDouble(d);
        currentPos += 8;
    }

    public void writeRow(int row, double[] cols) throws IOException {
        long seekLoc = ((long) row * bytesPerRow) + headerLen;

        if (seekLoc > channel.size()) {
            throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + channel.size());
        }

        // if row is the next row, just write.
        if (seekLoc - currentPos == 0) {

            writeRow(cols);
        } else {
            // else, random access to new location

            channel.position(seekLoc);
            if (bytebuffer == null) {
                bytebuffer = ByteBuffer.wrap(new byte[bytesPerRow]);

            }
            channel.write(bytebuffer);

            // this is probably extremely slow?
            counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
            os = new DataOutputStream(counter);
            currentPos = seekLoc + bytesPerRow;


        }
    }

    ByteBuffer singledouble;

    public void write(int row, int col, double val) throws IOException {
        long seekLoc = ((long) row * bytesPerRow) + headerLen + (col * 8);
//        System.out.println(row + "\t" + col + "\t" + seekLoc + "\t" + currentPos + "\t" + val);
        if (seekLoc - currentPos == 0) {

            os.writeDouble(val);
            os.flush();
            currentPos = seekLoc + 8;

        } else {
            if (seekLoc > channel.size()) {
                throw new IllegalArgumentException("Seek location for row: " + row + ", " + seekLoc + " is outside file size: " + channel.size());
            }

            if (singledouble == null) {
                singledouble = ByteBuffer.allocate(8);
            }

            singledouble.putDouble(val);
            singledouble.flip();

//            System.out.println("Seeking: " + seekLoc);
            channel.position(seekLoc);

            channel.write(singledouble);
            currentPos = seekLoc + 8;
            singledouble.compact();


            // this is probably extremely slow?
            counter = new CountingOutputStream(new BufferedOutputStream(Channels.newOutputStream(channel), buffersize));
            os = new DataOutputStream(counter);
        }

    }

    ByteBuffer blockbuffer;

    public void writeBlock(int startRow, int startCol, double[] vals) throws IOException {
        long seekLoc = ((long) startRow * bytesPerRow) + headerLen + (startCol * 8);
        System.out.println(startRow + "\t" + startCol + "\t" + seekLoc + "\t" + vals[0]);
        channel.position(seekLoc);

        if (blockbuffer == null || blockbuffer.limit() != vals.length * 8) {
            blockbuffer = ByteBuffer.allocate(vals.length * 8);
        }

        for (int b = 0; b < vals.length; b++) {
            blockbuffer.putDouble(vals[b]);
        }
        blockbuffer.flip();

        channel.write(blockbuffer);
        blockbuffer.compact();
    }

    public void close() throws IOException {
        os.close();
        counter.close();
        channel.close();
    }

}
