package nl.harmjanwestra.playground.biogen.freeze2;

import umcg.genetica.io.text.TextFile;

import java.io.IOException;

public class SettingsFileWriter {

	public static void main(String[] args) {
		String filein = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\settingsfiles\\combos-EUR.txt";
		String fileout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\settingsfiles\\combos-EUR.xml";

		try {
			TextFile tf = new TextFile(filein, TextFile.R);

			tf.readLine();

			TextFile tfo = new TextFile(fileout, TextFile.W);
			String[] elems = tf.readLineElems(TextFile.tab);
			while (elems != null) {

				String name = elems[0];
				String gt = elems[1];
				String gte = elems[2];
				String exp = elems[3];
				String annot = elems[4];
				String platform = elems[5];


				tfo.writeln("<dataset>");
				tfo.writeln("\t<name>" + name + "</name>");
				tfo.writeln("\t<location>" + gt + "</location>");
				tfo.writeln("\t<genometoexpressioncoupling>" + gte + "</genometoexpressioncoupling>");
				tfo.writeln("\t<expressiondata>" + exp + "</expressiondata>");
				tfo.writeln("\t<probeannotation>" + annot + "</probeannotation>");
				tfo.writeln("\t<expressionplatform>" + platform + "</expressionplatform>");
				tfo.writeln("\t<covariates></covariates>");
				tfo.writeln("\t<quantilenormalize>false</quantilenormalize>");
				tfo.writeln("\t<logtranform>false</logtranform>");
				tfo.writeln("</dataset>");

				elems = tf.readLineElems(TextFile.tab);
			}


			tfo.close();
			tf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}
}
