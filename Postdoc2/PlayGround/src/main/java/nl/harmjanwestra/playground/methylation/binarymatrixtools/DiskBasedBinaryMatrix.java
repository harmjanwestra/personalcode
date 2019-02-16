package nl.harmjanwestra.playground.methylation.binarymatrixtools;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.HashMap;

public class DiskBasedBinaryMatrix {
	
	public String[] rowIds;
	int nrRows;
	public String[] colIds;
	int nrCols;
	long bytesPerRow;
	HashMap<String, Integer> rowIndex;
	HashMap<String, Integer> colIndex;
	RandomAccessFile filehandle;
	long currentPos;
	
	
	public DiskBasedBinaryMatrix(File path) throws IOException {

//		BinaryFile b = new BinaryFile(path, BinaryFile.R);
		
//		filehandle = new RandomAccessFile(path, "r");
//		FileChannel gtChannel = filehandle.getChannel();
//
//		bytesPerRow = nrCols * Float.SIZE;
		
	}
	
	public void get(int row, int col) {
		long rowLoc = row * bytesPerRow;
	}
	
	
	public void getRow(int row) {
		
		
//		gtChannel.position(seekLoc);
//		if (mappedGenotypeHandle == null || mappedGenotypeHandle.capacity() != maplentouse) {
//			mappedGenotypeHandle = ByteBuffer.allocateDirect(maplentouse);
//		} else {
//			((Buffer) mappedGenotypeHandle).clear();
//		}
	}
	
}
