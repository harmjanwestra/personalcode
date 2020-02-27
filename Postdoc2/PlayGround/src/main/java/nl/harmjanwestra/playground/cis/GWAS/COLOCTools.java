package nl.harmjanwestra.playground.cis.GWAS;

import org.apache.commons.cli.*;

import java.io.IOException;

public class COLOCTools {

	private static final Options OPTIONS;

	static {
		OPTIONS = new Options();

		Option option = Option.builder()
				.longOpt("mergegwaswitheqtls")
				.desc("Merge GWAS with eQTLs")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("locusplot")
				.desc("Make locusplot")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("correlategenes")
				.desc("Correlate genes between traits")
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("colocparser")
				.desc("Merge coloc outputs, provide summary stats per gene")
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.longOpt("eqtl")
				.desc("eQTL mapping pipeline formatted eqtl file location")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("eqtlsummarystats")
				.desc("eQTL summary stats (mainly MAF)")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("gwas")
				.desc("PLINK GWAS file location")
				.hasArg()
				.build();
		OPTIONS.addOption(option);


		option = Option.builder()
				.longOpt("input")
				.desc("Input File/Dir")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("output")
				.desc("Output File/Dir")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("colocsummaryout")
				.desc("COLOC summary output")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("colocgeneoutput")
				.desc("COLOC gene summary output")
				.hasArg()
				.build();
		OPTIONS.addOption(option);

		option = Option.builder()
				.longOpt("ppthreshold")
				.desc("COLOC posterior probability threshold for calling a gene colocalized (default: 0.2)")
				.hasArg()
				.build();
		OPTIONS.addOption(option);
	}

	public static void main(String[] args) {
		COLOCTools c = new COLOCTools(args);
	}

	public COLOCTools(String[] args) {
		try {
			CommandLineParser parser = new DefaultParser();
			final CommandLine cmd = parser.parse(OPTIONS, args, false);


			if (cmd.hasOption("mergegwaswitheqtls")) {
				if (cmd.hasOption("eqtl") && cmd.hasOption("gwas") && cmd.hasOption("eqtlsummarystats") && cmd.hasOption("output")) {
					MergeWithGWAS v = new MergeWithGWAS();

					v.merge(cmd.getOptionValue("eqtl"),
							cmd.getOptionValue("gwas"),
							null,
							cmd.getOptionValue("eqtlsummarystats"),
							cmd.getOptionValue("output"));

				} else {
					System.out.println("mergegwaswitheqtls needs --eqtl --gwasfile --eqtlsummarystats and --output");
				}
			} else if (cmd.hasOption("colocparser")) {
				if (cmd.hasOption("input") && cmd.hasOption("colocsummaryout") && cmd.hasOption("colocgeneoutput")) {
					COLOCParser p = new COLOCParser();

					p.parseSummariesInDir(cmd.getOptionValue("input"), cmd.getOptionValue("colocsummaryout"));

					double PPthreshold = 0.2;
					if (cmd.hasOption("ppthreshold")) {
						PPthreshold = Double.parseDouble(cmd.getOptionValue("ppthreshold"));
					}
					p.identifyLociWithColocalizingSNPs(cmd.getOptionValue("input"), cmd.getOptionValue("colocgeneoutput"), PPthreshold);

				} else {
					System.out.println("colocparser needs --input --colocgeneoutput and --colocsummaryout");
				}


			} else {
				printHelp();
			}

		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(" ", OPTIONS);

		System.exit(-1);
	}

}
