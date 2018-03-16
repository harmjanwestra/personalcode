///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package misc;
//
//import java.nio.ByteBuffer;
//import java.util.BitSet;
//import jsc.util.BitVector;
//
///**
// *
// * @author harmjan
// */
//public class BitSetTest {
//
//    /**
//     * @param args the command line arguments
//     */
//    public static void main(String[] args) {
//        // TODO code application logic here
//        short[] testSet = new short[10];
//        int requiredBits = 11 * 10;
//
//        int remainder = 8 - (requiredBits % 8);
//        requiredBits += remainder; // should use padding at the end...
//
//        System.out.println("Nr bits: " + requiredBits + "\t" + (requiredBits / 8) + " bytes");
//
//        BitSet finalBits = new BitSet(requiredBits);
//        int nextBit = 0;
//        for (int i = 0; i < testSet.length; i++) {
//            double val = Math.round(Math.random() * 2001d);
//            testSet[i] = (short) val;
//
//            ByteBuffer b = ByteBuffer.allocate(2);
//            b.putShort(testSet[i]);
//
//            boolean[] byte1 = BitSetTest.bits(b.array()[0]);
//            boolean[] byte2 = BitSetTest.bits(b.array()[1]);
//
//            for (int q = 0; q < byte1.length; q++) {
//                finalBits.set(nextBit, byte1[q]);
//                nextBit++;
//            }
//
//            for (int q = 0; q < byte2.length - 5; q++) {
//                finalBits.set(nextBit, byte2[q]);
//                nextBit++;
//            }
//        }
//
//        System.out.println(nextBit);
//        byte[] bytes = toByteArray(finalBits);
//
//        System.out.println(bytes.length + "\tbytes finally produced...");
//        // i know there are 10 values. that would leave me with requiredbits amount of bits
//
//        BitSet bits = fromByteArray(new byte[]{bytes[0], bytes[1]});
//
//        System.out.println(bits.length() + " bits read");
//        for (int i = 11; i < bytes.length * 8; i++) {
//            bits.clear(i);
//        }
//        byte[] origBytes = toByteArray(bits);
//
//        System.out.println(origBytes.length);
////        short newshort = (short) ((origBytes[1] << 8) + (origBytes[0] & 0xFF));
//
//
//
//        ByteBuffer b = ByteBuffer.wrap(bytes);
//        short output = b.getShort();
//        System.out.println(testSet[0] + "\t" + output);
//
//
//    }
//
//    static boolean[] bits(final byte b) {
//        return new boolean[]{
//                    (b & 1) != 0,
//                    (b & 2) != 0,
//                    (b & 4) != 0,
//                    (b & 8) != 0,
//                    (b & 0x10) != 0,
//                    (b & 0x20) != 0,
//                    (b & 0x40) != 0,
//                    (b & 0x80) != 0
//                };
//    }
//
//    public static byte[] toByteArray(BitSet bits) {
//        byte[] bytes = new byte[bits.length() / 8 + 1];
//        for (int i = 0; i < bits.length(); i++) {
//            if (bits.get(i)) {
//                bytes[bytes.length - i / 8 - 1] |= 1 << (i % 8);
//            }
//        }
//        return bytes;
//    }
//
//    public static BitSet fromByteArray(byte[] bytes) {
//        System.out.println((bytes.length * 8) + " required for byte -> bitconversion");
//        BitSet bits = new BitSet();
//        for (int i = 0; i < bytes.length * 8; i++) {
//            if ((bytes[bytes.length - i / 8 - 1] & (1 << (i % 8))) > 0) {
//                bits.set(i);
//            }
//        }
//        return bits;
//    }
//}
