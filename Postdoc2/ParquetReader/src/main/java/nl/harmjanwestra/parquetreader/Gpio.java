package nl.harmjanwestra.parquetreader;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class Gpio {
    public Gpio() {
    }

    public static String[] getListOfFiles(String dir, String extension, String containsThis, String doesntContainThis) {
        File loc = new File(dir);
        String[] fileList = loc.list();
        ArrayList<String> al = new ArrayList();
        String[] output = fileList;
        int var8 = fileList.length;

        for(int var9 = 0; var9 < var8; ++var9) {
            String fileloc = output[var9];
            File file = new File(fileloc);
            if (!file.isDirectory()) {
                String fileName = file.getName();
                int mid = fileName.lastIndexOf(".");
                if (mid >= 0) {
                    String fname = fileName.substring(0, mid);
                    if (fname.contains(containsThis) && (doesntContainThis.isEmpty() || !fname.contains(doesntContainThis))) {
                        String ext = fileName.substring(mid + 1, fileName.length());
                        if (ext.toLowerCase().equals(extension)) {
                            al.add(loc.getAbsolutePath() + getFileSeparator() + fileName);
                        }
                    }
                }
            }
        }

        output = (String[])al.toArray(new String[0]);
        return output;
    }

    public static String[] getListOfFiles(String dir, String extension, String containsThis) {
        return getListOfFiles(dir, extension, containsThis, "");
    }

    public static String[] getListOfFiles(String dir, String extension) {
        return getListOfFiles(dir, extension, "", "");
    }

    public static void createDir(String dirName) throws IOException {
        if (!dirName.endsWith(getFileSeparator())) {
            dirName = dirName + getFileSeparator();
        }

        if (!exists(dirName)) {
            boolean success = (new File(dirName)).mkdirs();
            if (success) {
                System.out.println("Directory: " + dirName + " created");
            }
        }

    }

    public static boolean isDir(String dir) {
        File loc = new File(dir);
        return loc.isDirectory();
    }

    public static boolean exists(String dir) {
        return existsAndReadable(new File(dir));
    }

    public static boolean existsAndReadable(File file) {
        return file.exists() && file.canRead();
    }

    public static boolean createOuputDir(File dir) throws IOException {
        if (!dir.exists() && !dir.mkdirs()) {
            throw new IOException("Failed to create output dir at: " + dir.getAbsolutePath());
        } else {
            return dir.canRead() && dir.canWrite();
        }
    }

    public static void copyFile(File sourceFile, File destFile) throws IOException {
        if (!destFile.exists()) {
            destFile.createNewFile();
        }

        FileChannel source = null;
        FileChannel destination = null;

        try {
            source = (new FileInputStream(sourceFile)).getChannel();
            destination = (new FileOutputStream(destFile)).getChannel();
            destination.transferFrom(source, 0L, source.size());
        } finally {
            if (source != null) {
                source.close();
            }

            if (destination != null) {
                destination.close();
            }

        }

    }

    public static String[] getListOfFiles(String dir) {
        File loc = new File(dir);
        String[] fileList = loc.list();
        return fileList;
    }

    public static long getFileSize(String fileName) {
        File f = new File(fileName);
        long size = f.length();
        return size;
    }

    public static String humanizeFileSize(long s) {
        String output = "";
        DecimalFormat df = new DecimalFormat("##.##");
        if (s == 0L) {
            output = "0b";
        } else {
            double nrtb;
            if (s > 1099511627776L) {
                nrtb = (double)s / 1.099511627776E12;
                output = df.format(nrtb) + " TB";
            } else if (s > 1073741824L) {
                nrtb = (double)s / 1.073741824E9;
                output = df.format(nrtb) + " GB";
            } else if (s > 1048576L) {
                nrtb = (double)s / 1048576.0;
                output = df.format(nrtb) + " MB";
            } else if (s > 1024L) {
                nrtb = (double)s / 1024.0;
                output = df.format(nrtb) + " KB";
            }
        }

        return output;
    }

    public static void moveFile(String fileIn, String fileOut) throws IOException {
        File f1 = new File(fileIn);
        File f2 = new File(fileOut);
        copyFile(f1, f2);
        f1.delete();
    }

    public static void delete(String file) {
        File f = new File(file);
        f.delete();
    }

    public static boolean canRead(String fileName) {
        File f = new File(fileName);
        return f.canRead();
    }

    public static boolean canWrite(String fileName) {
        File f = new File(fileName);
        return f.canWrite();
    }

    public static boolean canExecute(String fileName) {
        File f = new File(fileName);
        return f.canExecute();
    }

    public static String getParentDir(String file) {
        File f = new File(file);
        return f.getParent();
    }

    public static String getFileName(String file) {
        File f = new File(file);
        return f.getName();
    }

    public static String getFileSeparator() {
        return System.getProperty("file.separator");
    }

    public static String formatAsDirectory(String loc) {
        if (!loc.endsWith(getFileSeparator())) {
            loc = loc + getFileSeparator();
        }

        return loc;
    }

    public static void copyFile(String filein, String fileout) throws IOException {
        copyFile(new File(filein), new File(fileout));
    }
}
