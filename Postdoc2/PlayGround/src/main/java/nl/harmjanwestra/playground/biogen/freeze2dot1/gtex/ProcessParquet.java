package nl.harmjanwestra.playground.biogen.freeze2dot1.gtex;


import ch.qos.logback.classic.Level;

import org.apache.parquet.column.page.PageReadStore;
import org.apache.parquet.example.data.simple.SimpleGroup;
import org.apache.parquet.example.data.simple.convert.GroupRecordConverter;
import org.apache.parquet.hadoop.ParquetFileReader;
import org.apache.parquet.hadoop.util.HadoopInputFile;
import org.apache.parquet.io.ColumnIOFactory;
import org.apache.parquet.io.MessageColumnIO;
import org.apache.parquet.io.RecordReader;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.Type;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.stream.IntStream;


public class ProcessParquet {
    public static void main(String[] args) {

        String folder = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_EUR_eQTL_all_associations\\";
        String outfolder = "U:\\2020-GTExV8\\GTEx_Analysis_v8_QTLs\\GTEx_Analysis_v8_EUR_eQTL_all_associations-flat\\";
        ProcessParquet p = new ProcessParquet();
        try {
            p.run(folder, outfolder);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run(String infolder, String outfolder) throws IOException {

        ch.qos.logback.classic.Logger root = (ch.qos.logback.classic.Logger) org.slf4j.LoggerFactory.getLogger(ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
        root.setLevel(Level.INFO);

        String[] files = Gpio.getListOfFiles(infolder);
        Arrays.sort(files);

        IntStream.range(0, files.length).parallel().forEach(v -> {
            String filename = files[v];
            try {
                processFile(infolder + filename, outfolder + filename + ".txt.gz");
            } catch (IOException e) {
                e.printStackTrace();
            }

        });

    }

    private void processFile(String infile, String outfile) throws IOException {


        Logger log = LoggerFactory.getLogger(ProcessParquet.class);
        log.info("Converting: " + infile + " --> " + outfile);
        TextFile output = new TextFile(outfile, TextFile.W);
        output.writeln("gene\tsnp\tmaf\tz\tp");

        org.apache.hadoop.fs.Path path = new org.apache.hadoop.fs.Path(infile);
        org.apache.hadoop.conf.Configuration config = new org.apache.hadoop.conf.Configuration(true);
        ParquetFileReader reader = ParquetFileReader.open(HadoopInputFile.fromPath(path, config));

        MessageType schema = reader.getFooter().getFileMetaData().getSchema();
        List<Type> fields = schema.getFields();

        int records = 0;
        PageReadStore pages;
        while ((pages = reader.readNextRowGroup()) != null) {
            long rows = pages.getRowCount();
            MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);

            RecordReader recordReader = columnIO.getRecordReader(pages, new GroupRecordConverter(schema));


            for (int i = 0; i < rows; i++) {
                SimpleGroup simpleGroup = (SimpleGroup) recordReader.read();


                String gene = simpleGroup.getString(0, 0);
                String snp = simpleGroup.getString(1, 0);
                float maf = simpleGroup.getFloat(3, 0);

                if (maf > 0) {
                    double p = simpleGroup.getDouble(6, 0);
                    double z = 0;
                    if (p < 1) {
                        z = (double) simpleGroup.getFloat(7, 0) / simpleGroup.getDouble(8, 0);
                    }
                    output.writeln(gene + "\t" + snp + "\t" + maf + "\t" + z + "\t" + p);
                } else {
                    output.writeln(gene + "\t" + snp + "\t" + 0 + "\t" + 0 + "\t" + 1.0);
                }

                records++;
                if (records % 1000000 == 0) {
                    log.info(infile + "\t" + records);
                }
            }

        }
        reader.close();
        output.close();
        log.info("Converting: " + infile + " --> " + outfile + "\tDONE.");

    }


}
