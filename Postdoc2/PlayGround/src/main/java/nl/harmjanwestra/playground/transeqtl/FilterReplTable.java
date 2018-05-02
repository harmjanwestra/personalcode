package nl.harmjanwestra.playground.transeqtl;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.HashSet;

public class FilterReplTable {
	
	
	public static void main(String[] args) {
		FilterReplTable r = new FilterReplTable();
		try {
			r.run();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void run() throws IOException {
		String replfile = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-03-ReplicationTables\\2018-04-03-transEQTL-Replication-PBMCAndWholeBlood.txt.gz";
		String out = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-03-ReplicationTables\\2018-04-03-transEQTL-Replication-PBMCAndWholeBlood-WithoutReplicatingInPurifiedCellTypes.txt";
		String out2 = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-03-ReplicationTables\\2018-04-03-transEQTL-Replication-PBMCAndWholeBlood-OnlyReplicatingInPurifiedCellTypes.txt";
		String tpfile = "D:\\Sync\\SyncThing\\Postdoc2\\2017-11-eQTLMeta\\data\\2018-04-03-ReplicationTables\\2016-04-03-transEQTL-ReplicatingInPurifiedCellTypes.txt";
		HashSet<String> tp = new HashSet<String>();
		TextFile tf = new TextFile(tpfile, TextFile.R);
		String[] elems = tf.readLineElems(TextFile.tab);
		while (elems != null) {
			String q = elems[0] + "_" + elems[1];
			tp.add(q);
			elems = tf.readLineElems(TextFile.tab);
		}
		tf.close();
		
		TextFile tout = new TextFile(out, TextFile.W);
		TextFile tout2 = new TextFile(out2, TextFile.W);
		TextFile tin = new TextFile(replfile, TextFile.R);
		String head1 = tin.readLine();
		tout.writeln(head1);
		tout2.writeln(head1);
		String head2 = tin.readLine();
		tout.writeln(head2);
		tout2.writeln(head2);
		elems = tin.readLineElems(TextFile.tab);
		while (elems != null) {
			String q = elems[0] + "_" + elems[3];
			if (!tp.contains(q)) {
				tout.writeln(Strings.concat(elems, Strings.tab));
			} else {
				tout2.writeln(Strings.concat(elems, Strings.tab));
			}
			elems = tin.readLineElems(TextFile.tab);
		}
		tin.close();
		tout.close();
		tout2.close();
	}
	
	
}
