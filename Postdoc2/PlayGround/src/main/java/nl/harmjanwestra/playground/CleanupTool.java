package nl.harmjanwestra.playground;

import org.apache.commons.io.FilenameUtils;
import umcg.genetica.io.bin.BinaryFile;
import umcg.genetica.io.text.TextFile;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

import org.apache.commons.codec.binary.Hex;

public class CleanupTool {

	public static void main(String[] args) {

		if (args.length < 1) {
			System.out.println("Usage: indir [commaseparatedlistofextensions]");
			System.out.println("Default extensions: txt tsv csv tab dat log bed gtf fa");
			System.exit(0);
		}

		String dir = args[0];
		String[] ext = null;
		if (args.length > 1) {
			ext = args[1].split(",");
		}
//		String dir = "D:\\asd\\";
		CleanupTool t = new CleanupTool();
		try {
			t.run(dir, ext);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private static String createMD5(final InputStream is) throws EOFException, IOException, NoSuchAlgorithmException {
		MessageDigest md = MessageDigest.getInstance("MD5");
		md.reset();
		int byteArraySize = 8096;
		byte[] bytes = new byte[byteArraySize];
		int numBytes;
		while ((numBytes = is.read(bytes)) != -1) {
			md.update(bytes, 0, numBytes);
		}
		byte[] digest = md.digest();
		String result = new String(Hex.encodeHex(digest));
		return result;
	}

	private static void debug(Path v, Path gv) throws IOException {
		System.out.println("Comparing " + v.toString() + " vs " + gv.toString());
		InputStream is = Files.newInputStream(v);
		GZIPInputStream gis = new GZIPInputStream(Files.newInputStream(gv));
		byte[] buffer0 = new byte[100];
		byte[] buffer1 = new byte[100];
		is.read(buffer0);
		gis.read(buffer1);
		System.out.println("Compare...");
	}

	private static boolean compare(Path v, Path gv) throws IOException, EOFException {
		System.out.println("Comparing " + v.toString() + " vs " + gv.toString());
		InputStream is = Files.newInputStream(v);
		String md5 = null;
		try {
			md5 = createMD5(is);
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
		is.close();

		GZIPInputStream gis = new GZIPInputStream(Files.newInputStream(gv));
		String md5gz = null;
		try {
			md5gz = createMD5(gis);
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
		gis.close();

		if (md5.equals(md5gz)) {
			System.out.println("EQUAL: " + v.toString() + " md5: " + md5 + " vs " + gv.toString() + " md5: " + md5gz);
		} else {
			System.out.println("UNEQUAL: " + v.toString() + " md5: " + md5 + " vs " + gv.toString() + " md5: " + md5gz);
		}
		return md5gz.equals(md5);
	}

	HashSet<Path> visited = new HashSet<>();

	public synchronized boolean visited(Path v) {
		if (!visited.contains(v)) {
			visited.add(v);
			return false;
		}
		return true;
	}

	private void run(String dir, String[] ext) throws IOException {
		TextFile errorlog = new TextFile(dir + "/errors.txt.gz", TextFile.W);

		HashSet<String> extset = new HashSet<>();
		if (ext == null) {
			extset.add("txt");
			extset.add("tsv");
			extset.add("tab");
			extset.add("csv");
			extset.add("dat");
			extset.add("log");
			extset.add("bed");
			extset.add("gtf");
			extset.add("fa");
		} else {
			for (String s : ext) {
				extset.add(s);
			}
		}

		try (Stream<Path> paths = Files.walk(Paths.get(dir))) {
			paths.parallel().forEach(v -> {
				boolean visited = visited(v);

				if (!Files.isDirectory(v) && !visited) {

					String fStr = v.toString();
					String gfStr = fStr + ".gz";
					Path gv = Paths.get(gfStr);

					try {


						if (!fStr.equals(dir + "/errors.txt.gz")) {
							String currentExt = FilenameUtils.getExtension(fStr);
							if (extset.contains(currentExt) && !currentExt.equals("gz") && !currentExt.equals("zip")) {
								boolean gzippedVersionPresent = Files.exists(Paths.get(gfStr));
								if (gzippedVersionPresent) {
									boolean ok = compare(v, gv);
									if (ok) {
										System.out.println("Deleting " + fStr);
										Files.delete(v);
									}
								} else {
									// gzip the file
									System.out.println("Gzipping " + fStr);
									BinaryFile bf = new BinaryFile(fStr, BinaryFile.R);
									BinaryFile bfo = new BinaryFile(fStr + ".gz", BinaryFile.W);
									int totalSize;
									byte[] buffer = new byte[8096];
									while ((totalSize = bf.read(buffer)) > 0) {
										bfo.write(buffer, 0, totalSize);
									}
									bfo.close();
									bf.close();

									// mark the file for deletion.
									boolean ok = compare(v, gv);
									if (ok) {
										System.out.println("Deleting " + fStr);
										Files.delete(v);
									}
								}
							}
						}

					} catch (EOFException e2) {
						System.out.println("Unexpected EOF: " + fStr + " or " + gfStr);
						e2.printStackTrace();
						try {
							errorlog.writelnsynced("EOF problem: " + fStr + " or " + gfStr);
						} catch (IOException e) {
							throw new RuntimeException(e);
						}
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
			});
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		errorlog.close();
	}

}
