package nl.harmjanwestra.playground;

import net.jpountz.lz4.LZ4Compressor;
import net.jpountz.lz4.LZ4Factory;
import umcg.genetica.util.RunTimer;

import java.util.BitSet;

public class LZ4 {
	
	public static void main(String[] args) {
		
		
		RunTimer rt = new RunTimer();
		rt.start();
		
		long bytescomp = 0;
		long bytesdecomp = 0;
//		for (int perm = 0; perm < 1000; perm++) {
//
//			short[][] matrix = new short[10000][20000];
//			for (int i = 0; i < matrix.length; i++) {
//				for (int j = 0; j < matrix[i].length; j++) {
//					matrix[i][j] = (short) (Math.random() * 1000);
//				}
//			}
//
//			LZ4Factory factory = LZ4Factory.fastestInstance();
//
//			RunTimer rt2 = new RunTimer();
//			rt2.start();
//
//
//			byte[] concatenated = concat(matrix);
//
//			final int decompressedLength = concatenated.length;
//			bytesdecomp += decompressedLength;
//
//			// compress data
//			LZ4Compressor compressor = factory.fastCompressor();
//			int maxCompressedLength = compressor.maxCompressedLength(decompressedLength);
//			byte[] compressed = new byte[maxCompressedLength];
//			int compressedLength = compressor.compress(concatenated, 0, decompressedLength, compressed, 0, maxCompressedLength);
//			bytescomp += compressedLength;
//			rt2.stop();
//			double mbcomp = (double) compressedLength / 1048576;
//			double mbdecomp = (double) decompressedLength / 1048576;
//			System.out.println(perm + "\t" + rt2.getTimeDesc() + "\t" + mbcomp + " comp\t" + mbdecomp + " decomp\tratio: " + (mbcomp / mbdecomp));
//
//
//		}
		
		rt.stop();

//		System.out.println(rt.getTimeDesc());
//		System.out.println((double) bytescomp / 1048576 + " mb compressed");
//		System.out.println((double) bytesdecomp / 1048576 + " mb decompressed");
//		System.out.println((double) bytescomp / bytesdecomp + " avg ratio");
		
		
		Bits b = new Bits();
		for (int i = 0; i < 2001; i++) {
			BitSet bitset = b.convert(i);
			System.out.println(i + "\t" + bitset.toString());
		}
		
		
	}
	
	private static byte[] concat(short[][] matrix) {
		int len = matrix.length * matrix[0].length * 2; // 2 bytes per short
		byte[] concate = new byte[len];
		int pos = 0;
		for (int i = 0; i < matrix.length; i++) {
			byte[] conv = convertShortArr(matrix[i]);
			System.arraycopy(conv, 0, concate, pos, conv.length);
			pos += conv.length;
		}
		return concate;
	}
	
	public static byte[] convertShort(short x) {
		byte[] ret = new byte[2];
		ret[0] = (byte) x;
		ret[1] = (byte) (x >> 8);
		return ret;
	}
	
	public static byte[] convertShortArr(short[] x) {
		byte[] ret = new byte[x.length * 2];
		for (int i = 0; i < x.length; i++) {
			short s = x[i];
			int idx = i * 2;
			ret[idx] = (byte) (s >> 8);
			ret[idx + 1] = (byte) s;
		}
		return ret;
	}
	
	
	public static class Bits {
		
		public BitSet convert(long value) {
			BitSet bits = new BitSet();
			int index = 0;
			while (value != 0L) {
				if (value % 2L != 0) {
					bits.set(index);
				}
				++index;
				value = value >>> 1;
			}
			return bits;
		}
		
		public long convert(BitSet bits) {
			long value = 0L;
			for (int i = 0; i < bits.length(); ++i) {
				value += bits.get(i) ? (1L << i) : 0L;
			}
			return value;
		}
	}
}
