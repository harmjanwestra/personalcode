/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package encodify;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.SocketException;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.commons.net.ftp.FTPReply;
import umcg.genetica.io.Gpio;

/**
 *
 * @author harmjan
 */
public class Encodify {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            Encodify f = new Encodify();
            f.run();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
//    public void run() {
    FTPClient ftp = null;
    int reply = -1;
    String username = "anonymous";
    String passwd = "marco@live.nl";
    String address = "hgdownload.cse.ucsc.edu";
    String startWD = "/apache/htdocs/goldenPath/hg19/encodeDCC/";
    String startLocalWD = "/Volumes/BackupDisk/Encode/hg19/";
    long totalBytes = 0;
    int nrLongsUsed = 0;

    public Encodify() {
        ftp = new FTPClient();
    }

    public void run() throws IOException {

        // first check whether a previous run has collected all files

        boolean connected = connect(username, passwd, address);

        System.out.println("Buffer: " + ftp.getBufferSize());
        ftp.setFileTransferMode(FTP.STREAM_TRANSFER_MODE);
        ftp.setFileType(FTP.BINARY_FILE_TYPE);
        ftp.setControlKeepAliveTimeout(300);
        // set transfer mode to binary
        ftp.setBufferSize(8192);
        System.out.println("Buffer: " + ftp.getBufferSize());

        if (connected) {
            // change directory 
            iterateDir(startWD, startLocalWD);
            ftp.disconnect();
        }
        // list directories

        System.out.println("");
        System.out.println("--------------------------------------");
        System.out.println("Nr Longs used: " + nrLongsUsed);
        System.out.println("Max long: " + convertSize(Long.MAX_VALUE));
        System.out.println("Remainder: " + convertSize(totalBytes));
        System.out.println("--------------------------------------");
        System.out.println("");


        // iterate each directory, only go one directory deep..
        // get wig/bigwig/broadpeak/narrowpeak data
        // binarify the data

    }

    private void iterateDir(String inputDir, String localDir) throws IOException {
        System.out.println("");
        System.out.println("--------------------------------------");
        System.out.println("Nr Longs used: " + nrLongsUsed);
        System.out.println("Max long: " + convertSize(Long.MAX_VALUE));
        System.out.println("Remainder: " + convertSize(totalBytes));
        System.out.println("--------------------------------------");
        System.out.println("");
        boolean success = ftp.changeWorkingDirectory(inputDir);
        System.out.println("");
        System.out.println("Current WD: \t" + ftp.printWorkingDirectory());
        System.out.println("Local Dir: \t" + localDir);
        if (success) {

            FTPFile[] files = ftp.listFiles();
            int fCtr = 0;
            int peakfiles = 0;
            int wigfiles = 0;
            int bedfiles = 0;
            for (FTPFile file : files) {
                if (!file.isDirectory()) {
                    String fName = file.getName();
                    String fNameLc = fName.toLowerCase();
                    if (fNameLc.contains("raw")) {
                        // dont use raw files
                    } else {
                        if (fNameLc.contains("wig")) {
                            if (!Gpio.exists(localDir)) {
                                Gpio.createDir(localDir);
                            }
                            wigfiles++;
//                            String localFileName = localDir + fName;
//                            String convertedFileName = fName + ".wigbin";
                            // check whether the parsed file is already there. 
//                            if (Gpio.exists(localDir + convertedFileName)) {
//
//                            } else {
                            // retrieve file
                            retrieve(file, localDir);

//                            System.out.println(fName);
                        } else if (fNameLc.contains("peak")) {
                            if (!Gpio.exists(localDir)) {
                                Gpio.createDir(localDir);
                            }
                            retrieve(file, localDir);
                            peakfiles++;
                        } else if (fNameLc.contains("bed")) {
                            if (!Gpio.exists(localDir)) {
                                Gpio.createDir(localDir);
                            }
//                            retrieve(file, localDir);
                            bedfiles++;
//                            System.out.println(fName);
                        } else {
                        }
                    }

//                    String filename = file.getName();
//                    System.out.println("File: " + inputDir + filename);
                    fCtr++;
                }
            }
            System.out.println("Nr Files: " + fCtr + "\tWIG: " + wigfiles + "\tPeak: " + peakfiles + "\tBED: " + bedfiles);


            FTPFile[] dirs = ftp.listDirectories();
            for (FTPFile dir : dirs) {
                String dirName = inputDir + dir.getName() + "/";
                if (dir.getName().equals(".") || dir.getName().equals("..")) {
                    // skip dir
                } else {
                    String nextlocalDir = localDir + dir.getName() + "/";
                    iterateDir(dirName, nextlocalDir);
                }
            }
        }
    }

    private void retrieve(FTPFile f, String localPath) throws IOException {
        System.out.println("About to retrieve file: " + f.getName() + "\t to location: " + localPath + f.getName());
        String fullLocalPath = localPath + f.getName();
        File file = new File(fullLocalPath);
        long fszise = f.getSize();
        if (totalBytes + fszise > Long.MAX_VALUE) {
            long remainder = Long.MAX_VALUE - totalBytes;
            fszise -= remainder;
            totalBytes = fszise;
            nrLongsUsed++;
        } else {
            totalBytes += fszise;
        }

        boolean skipfile = false;
        if (file.exists()) {
            System.out.print("\t- " + fullLocalPath + "\talready exists. Checking file size. ");
            long localSize = Gpio.getFileSize(fullLocalPath);
            long remoteSize = f.getSize();

            if (localSize == remoteSize) {
                skipfile = true;
                System.out.print("Sizes are equal. Skipping file.\n");
            } else {
                System.out.print("Sizes differ (" + convertSize(Math.abs(localSize - remoteSize)) + "). Redownloading file.\n");
            }

        }
        if (!skipfile) {
            boolean success = false;
            while (!success) {
                System.out.println("\t- Retrieving file: " + f.getName() + " (" + convertSize(f.getSize()) + ")");
                OutputStream outputStream = new FileOutputStream(file);


                success = ftp.retrieveFile(f.getName(), outputStream);
                outputStream.close();
                long localSize = Gpio.getFileSize(fullLocalPath);
                long remoteSize = f.getSize();
                if (Math.abs(localSize - remoteSize) > 0) {
                    System.out.println("\t- File downloaded succesfully but not complete: " + convertSize(Math.abs(localSize - remoteSize)) + " bytes difference");
//                    success = false;
                }
            }
        }
        ftp.noop();
    }

    public static String convertSize(long size) {
        if (size < 1024) {
            return "" + size + " b";
        } else if (size >= 1024 && size < 1048576) {
            return "" + (size / 1024) + " kb";
        } else if (size >= 1048576 && size < (1048576 * 1024)) {
            return "" + (size / 1048576) + " mb";
        } else {
            return "" + (size / (1048576 * 1024)) + " gb";
        }
    }

    private boolean connect(String username, String password, String serverAdd) {
        try {

            ftp.connect(serverAdd);
            ftp.login(username, password);
            reply = ftp.getReplyCode();

            if (FTPReply.isPositiveCompletion(reply)) {
                System.out.println("Connected Success");
                return true;
            } else {
                System.out.println("Connection Failed");
                ftp.disconnect();
                return false;
            }

        } catch (SocketException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        return false;

    }
//}
}
