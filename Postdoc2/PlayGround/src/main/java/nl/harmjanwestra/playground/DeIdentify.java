package nl.harmjanwestra.playground;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import umcg.genetica.io.trityper.util.BaseAnnot;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

public class DeIdentify {
	
	
	public static void main(String[] args) {
		
		
		if (args.length < 3) {
			System.out.println("Usage: referenceFasta bamin bamout");
		} else {
			DeIdentify d = new DeIdentify();
			try {
				d.run(args[1], args[2], args[0]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	public void run(String bamFile, String outputSamOrBamFile, String fasta) throws IOException {
		
		
		// read the reference genome
		LinkedHashMap<String, byte[]> referenceGenome = new LinkedHashMap<>();
		htsjdk.samtools.reference.FastaSequenceFile fastaFile = new htsjdk.samtools.reference.FastaSequenceFile(new File(fasta), false);
		ReferenceSequence seq = fastaFile.nextSequence();
		while (seq != null) {
			String name = seq.getName();
			byte[] bases = seq.getBases();
			referenceGenome.put(name, bases);
			System.out.println("Loaded " + name + " - len: " + bases.length);
			seq = fastaFile.nextSequence();
		}

//        System.out.println("Checking sequence library.");
//        {
//            SamReader reader = SamReaderFactory.make().enable(new SamReaderFactory.Option[]{SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES}).validationStringency(ValidationStringency.DEFAULT_STRINGENCY).open(new File(bamFile));
//
//            Integer lastseqindex = null;
//            for (SAMRecord samRecord : reader) {
//                Integer index = samRecord.getReferenceIndex();
//                if (index != null && index > -1) {
//                    if (lastseqindex == null || !lastseqindex.equals(index)) {
//                        SAMSequenceRecord sequence = reader.getFileHeader().getSequence(index);
//                        String seqname = sequence.getSequenceName();
//                        System.out.println("Looking for chr " + seqname + "\t" + index);
//                        byte[] refbases = referenceGenome.get(seqname);
//                        if (refbases == null) {
//                            System.err.println("Could not find sequence in reference: " + seqname);
//                            System.exit(-1);
//                        }
//                        lastseqindex = index;
//                    }
//                }
//
//            }
//            reader.close();
//        }
//
//        System.out.println("Done checking sequence library.");
		
		SamReader reader = SamReaderFactory.make().enable(new SamReaderFactory.Option[]{SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES}).validationStringency(ValidationStringency.DEFAULT_STRINGENCY).open(new File(bamFile));
		SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(),
				true, new File(outputSamOrBamFile));
		
		Integer lastseqindex = null;
		
		long baseschanged = 0;
		long basesN = 0;
		long totalbases = 0;
		long ctr = 0;
		long written = 0;
		byte N = 78;
		byte[] refbases = null;
		for (SAMRecord samRecord : reader) {
			
			Integer index = samRecord.getReferenceIndex();
			if (index != null && index > -1) {
				if (refbases == null || lastseqindex == null || !lastseqindex.equals(index)) {
					SAMSequenceRecord sequence = reader.getFileHeader().getSequence(index);
					String seqname = sequence.getSequenceName();
					System.out.println("Looking for " + seqname);
					refbases = referenceGenome.get(seqname);
					if (refbases != null) {
						System.out.println(seqname + " has " + refbases.length + " bases.");
					}
					lastseqindex = index;
				}
				
				if (refbases == null) {
					SAMSequenceRecord sequence = reader.getFileHeader().getSequence(index);
					String seqname = sequence.getSequenceName();
					System.err.println("Could not find sequence in reference: " + seqname);
					System.exit(-1);
				} else {
					List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();
					byte[] bases = samRecord.getReadBases();
					if (bases.length > 0 && blocks.size() > 0) {
						byte[] newbases = new byte[bases.length];
						boolean[] baseset = new boolean[bases.length];
						int delta = 0;
						for (AlignmentBlock block : blocks) {
							int refst = block.getReferenceStart() - 1;
							int len = block.getLength();
							int readst = block.getReadStart() - 1;
							int readsto = readst + len;
							int b = 0;
							try {
								for (int i = readst; i < readsto; i++) {
									int refpos = refst + b;
									byte refb = refbases[refpos];
									if (refb > 96) {
										refb = capitalize(refb);
									}
									newbases[i] = refb;
									// capitalize
									if (newbases[i] != bases[i]) {
										baseschanged++;
										delta++;
									}
									baseset[i] = true;
									b++;
								}
							} catch (NullPointerException e) {
								e.printStackTrace();
								
								System.out.println("Watte?");
							}
						}
						
						// fill the remainder with N
						for (int b = 0; b < newbases.length; b++) {
							if (!baseset[b]) {
								newbases[b] = N;
								basesN++;
							}
						}
						
						
						if (delta > 10) {
							String basestr = "";
							String basestr2 = "";
							for (int b = 0; b < bases.length; b++) {
								basestr += BaseAnnot.toString(bases[b]);
								basestr2 += BaseAnnot.toString(newbases[b]);
//                System.out.println(BaseAnnot.toString(bases[b]) + " --> " + BaseAnnot.toString(newbases[b]));
							}
							System.out.println();
							System.out.println(delta + " bases changed " + samRecord.getMappingQuality() + "\t" + samRecord.getCigarString() + "\t" + basestr + "\t" + basestr2);
							System.out.println();
						}
						samRecord.setReadBases(newbases);
						outputSam.addAlignment(samRecord);
						totalbases += bases.length;
						written++;
					}
				}
				
				
			}
			ctr++;
			if (ctr % 100000 == 0) {
				System.out.print("\r" + ctr + " reads parsed. " + written + " (" + ((double) written / ctr) + ") written.\t" + totalbases + "\t" + baseschanged + " (" + ((double) baseschanged / totalbases) + ") changed\t" + basesN + " bases Ned. Last index: " + lastseqindex);
			}
		}
		
		outputSam.close();
		reader.close();
		
	}
	
	private byte capitalize(byte b) {
		switch (b) {
			case 97: // a
				return 65;
			case 99: // c
				return 67;
			case 116: // t
				return 84;
			case 103: // g
				return 71;
		}
		return b;
	}
}
