package nl.harmjanwestra.playground;

import org.apache.commons.codec.binary.Hex;
import org.apache.commons.io.FilenameUtils;
import umcg.genetica.io.text.TextFile;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.AccessDeniedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

public class CleanupToolDirCompare {

	public static void main(String[] args) {

		if (args.length < 2) {
			System.out.println("Usage: olddir newdir [OK]");
			System.out.println("Remove files from olddir that are also in newdir by comparing md5sums. Defaults to a dry run. Add OK to actually delete files.");
			System.exit(0);
		}

		String olddir = args[0];
		String newdir = args[1];

		boolean dryrun = true;
		if (args.length > 2) {
			if (args[2].equals("OK")) {
				dryrun = false;
			}
		}


//		String olddir = "D:\\tmp\\dir1";
//		String newdir = "D:\\tmp\\dir2";

//		String dir = "D:\\asd\\";
		CleanupToolDirCompare t = new CleanupToolDirCompare();
		try {
			t.run(olddir, newdir, dryrun);
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

	private static boolean compare(Path oldFile, Path newFile) throws IOException, EOFException {
		System.out.println("Comparing " + oldFile.toString() + " vs " + newFile.toString());
		InputStream is = Files.newInputStream(oldFile);
		String oldMD5 = null;
		try {
			oldMD5 = createMD5(is);
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
		is.close();

		InputStream gis = Files.newInputStream(newFile);
		String newMD5 = null;
		try {
			newMD5 = createMD5(gis);
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
		gis.close();

		if (oldMD5.equals(newMD5)) {
			System.out.println("EQUAL: " + oldFile.toString() + " md5: " + oldMD5 + " vs " + newFile.toString() + " md5: " + newMD5);
		} else {
			System.out.println("UNEQUAL: " + oldFile.toString() + " md5: " + oldMD5 + " vs " + newFile.toString() + " md5: " + newMD5);
		}
		return newMD5.equals(oldMD5);
	}

	private static boolean compareWithGzippedVersion(Path v, Path gv) throws IOException, EOFException {
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

	public static final String ANSI_RESET = "\u001B[0m";
	public static final String ANSI_BLACK = "\u001B[30m";
	public static final String ANSI_RED = "\u001B[31m";
	public static final String ANSI_GREEN = "\u001B[32m";
	public static final String ANSI_YELLOW = "\u001B[33m";
	public static final String ANSI_BLUE = "\u001B[34m";
	public static final String ANSI_PURPLE = "\u001B[35m";
	public static final String ANSI_CYAN = "\u001B[36m";
	public static final String ANSI_WHITE = "\u001B[37m";

	private void run(String olddir, String newdir, boolean dryrun) throws IOException {
		System.out.println("Old dir: " + olddir);
		System.out.println("New dir: " + newdir);
		System.out.println("Performing dry run: " + dryrun);
		TextFile errorlog = new TextFile(newdir + "/errors.txt.gz", TextFile.W);

		HashSet<String> extset = new HashSet<>();
		extset.add("txt");
		extset.add("tsv");
		extset.add("tab");
		extset.add("csv");
		extset.add("dat");
		extset.add("log");
		extset.add("bed");
		extset.add("gtf");
		extset.add("fa");

		Path olddirPath = Paths.get(olddir);
		int olddirlen = olddirPath.getNameCount();
		System.out.println("Old dir path has " + olddirlen + " length");
		Path newdirPath = Paths.get(newdir);

		AtomicInteger nrDeleted = new AtomicInteger(0);
		AtomicInteger nrFiles = new AtomicInteger(0);
		AtomicInteger notFound = new AtomicInteger(0);

		try (Stream<Path> paths = Files.walk(olddirPath)) {
			paths.forEach(oldFile -> {
				boolean visited = visited(oldFile);
				int f = -1;
				int d = -1;
				if (!Files.isDirectory(oldFile) && !visited) {
					f = nrFiles.getAndIncrement();
					System.out.println("Old file: " + oldFile);
					// check if a directory with a similar name is present in newdir
					int plen = oldFile.getNameCount();
					Path relative = oldFile.subpath(olddirlen, plen);

//					System.out.println("Relative: " + relative);
					Path newFile = newdirPath.resolve(relative);

					try {

						System.out.println("Looking for: " + newFile);

						if (Files.exists(newFile) && !Files.isDirectory(newFile)) {
							System.out.println("Found: " + newFile);

							if (compare(oldFile, newFile)) {
								System.out.println(ANSI_GREEN + "Files are equal. Deleting old file: " + oldFile + ANSI_RESET);
								d = nrDeleted.getAndIncrement();
								if (!dryrun) {
									Files.delete(oldFile);
								} else {
									System.out.println("Dry run. Not actually removing file.");
								}
							} else {
								System.out.println(ANSI_PURPLE + "Files are NOT equal. Not deleting old file: " + oldFile + ANSI_RESET);
							}

						} else {
							System.out.println("Not found: " + newFile);
							String currentExt = FilenameUtils.getExtension(oldFile.toString());
							if (!Files.isDirectory(oldFile) && extset.contains(currentExt) && !currentExt.equals("gz") && !currentExt.equals("zip")) {
								// check if there is a gzipped version of the old file in the new directory
								System.out.println("Looking for GZipped version of " + oldFile);
								Path gzipNewFile = Paths.get(newFile.toString() + ".gz");
								if (Files.exists(gzipNewFile)) {

									System.out.println("Gzipped file exists: " + gzipNewFile);
									if (compareWithGzippedVersion(oldFile, gzipNewFile)) {
										System.out.println(ANSI_GREEN + "Gzipped file is equal. Deleting old file: " + oldFile + ANSI_RESET);
										d = nrDeleted.getAndIncrement();
										if (!dryrun) {
											Files.delete(oldFile);
										} else {
											System.out.println("Dry run. Not actually removing file.");
										}
									}

								}
							} else {
								notFound.getAndIncrement();
							}
						}


						System.out.println();

						if (f % 10 == 0) {
							d = nrDeleted.get();
							int nf = notFound.get();
							int tot = d + nf;
							System.out.println(ANSI_RED + f + " files processed, " + d + " deleted, " + nf + " not found, " + tot + " total." + ANSI_RESET);
							System.out.println();
						}
					} catch (AccessDeniedException e4){
						try {
							System.out.println("Error locating file: " + oldFile + " or " + newFile);
							errorlog.writelnsynced("AccessDeniedException: " + oldFile + " or " + newFile);
						} catch (IOException e) {
							throw new RuntimeException(e);
						}
					} catch (java.nio.file.NoSuchFileException e3) {
						try {
							System.out.println("Error locating file: " + oldFile + " or " + newFile);
							errorlog.writelnsynced("NoSuchFileException: " + oldFile + " or " + newFile);
						} catch (IOException e) {
							throw new RuntimeException(e);
						}
					} catch (EOFException e2) {
						System.out.println("Unexpected EOF: " + oldFile + " or " + newFile);
						e2.printStackTrace();
						try {
							errorlog.writelnsynced("EOF problem: " + oldFile + " or " + newFile);
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
