package nl.harmjanwestra.playground.biogen.gwasfilter;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static java.nio.charset.StandardCharsets.*;

public class GroupGWASTraits {

    public static void main(String[] args) {
        GroupGWASTraits t = new GroupGWASTraits();
        String input = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-11-15-GWASUpdateForTransFX\\2021-11-15-2021-11-15-gwaslist-wgwascatalog-wALS-wMetaBrain.txt.gz";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-11-15-GWASUpdateForTransFX\\2021-11-15-2021-11-15-gwaslist-wgwascatalog-wALS-wMetaBrain-braintraits.txt.gz";
        String assoc = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-11-15-GWASUpdateForTransFX\\2021-11-15-2021-11-15-allTopAssociations-wgwascatalog-wALS-wMetaBrain-MetaBrain2dot1IDs.txt.gz";
        String assocout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-11-15-GWASUpdateForTransFX\\2021-11-15-2021-11-15-allTopAssociations-wgwascatalog-wALS-wMetaBrain-MetaBrain2dot1IDs-braintraits.txt.gz";
        String tophitsout = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2021-06-04-Rebuttal\\2021-11-15-GWASUpdateForTransFX\\2021-11-15-2021-11-15-topHitsUniqueRSids-wgwascatalog-wALS-wMetaBrain-MetaBrain2dot1IDs-braintraits.txt.gz";
        try {
            t.run(input, output, assoc, assocout, tophitsout);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static String removeNonAlphanumeric(String str) {
        // replace the given string
        // with empty string
        // except the pattern "[^a-zA-Z0-9]"
//        str = str.replaceAll("[^a-zA-Z0-9\\s]", "");
        // return string
        return str;
    }

    public void run(String input, String output, String assoc, String assocout, String tophitsout) throws IOException {
        TextFile tf = new TextFile(input, TextFile.R);
        String header = tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        HashSet<String> traits = new HashSet<String>();
        HashMap<String, String> traitIdToTrait = new HashMap<>();
        while (elems != null) {
            String traitid = elems[0];
            traits.add(elems[3]);
            traitIdToTrait.put(traitid, elems[3]);
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        System.out.println(traits.size() + " traits loaded");
        ArrayList<String> allTraits = new ArrayList<>();
        allTraits.addAll(traits);
        Collections.sort(allTraits);
        HashMap<String, String> traitToAnnotation = new HashMap<>();
        int ctrnotincluded = 0;
        int braineqtltraits = 0;
        int braintraits = 0;
        int eqtltraits = 0;
        HashMap<String, Integer> wordcount = new HashMap<>();
        for (String trait : allTraits) {
            trait = trait.toLowerCase();

            if (trait.startsWith("metabrain")) {
                traitToAnnotation.put(trait, "metabrain");
                braineqtltraits++;
            } else if (trait.contains("brain")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("cerebral") || trait.contains("cerebell") || trait.contains("putamen")
                    || trait.contains("amygdala") || trait.contains("pituitary") || trait.contains("basal ganglia")
                    || trait.contains("hippocamp")
            ) {
                traitToAnnotation.put(trait, "brain");
//                System.err.println(trait);
                braintraits++;
            } else if (trait.contains("epilepsy")) {
                traitToAnnotation.put(trait, "brain");

                braintraits++;
            } else if (trait.contains("schooling") || trait.contains("education") || trait.contains("cognition") || trait.contains("cognitive")
                    || trait.contains("mathematical ability") || trait.contains("dyslexia") || trait.contains("intelligence")) {
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("alzheimer") || trait.contains("parkinson") || trait.contains("multiple sclerosis") ||
                    trait.contains("amyotropic")) {
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("dementia")) {
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("psycho")) {
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("irritability") || trait.contains("depression") || trait.contains("schizo") || trait.contains("depressive")
                    || trait.contains("bipolar") || trait.contains("autism") || trait.contains("anorexia") || trait.contains("adhd")
                    || trait.contains("attention deficit hyperactivity disorder")
                    || trait.contains("worry") || trait.contains("neuroticism") || trait.contains("neurotic") || trait.contains("tic disorders") || trait.contains("anxiety")
                    || trait.contains("social disorders") || trait.contains("sensation seeking")
                    || trait.contains("sensitivity / hurt feelings") || trait.contains("sensitivity to environmental stress and adversity")
                    || trait.contains("psychiatric") || trait.contains("personality disorders")
                    || trait.contains("persistent mood disorders") || trait.contains("persistent delusional disorders")
                    || trait.contains("obsessive compulsive") || trait.contains("obsessive-compulsive") || trait.contains("mood")
                    || trait.startsWith("mixed disorders of conduct and emotions") || trait.startsWith("mixed and other personality disorders")
                    || trait.contains("manic episode") || trait.contains("manic/hyper")
                    || trait.startsWith("mental") || trait.contains("behavior") || trait.contains("behavioural")
                    || trait.startsWith("general happiness") || trait.contains("feeling") || trait.contains("emotional")
                    || trait.contains("depressed affect") || trait.contains("dependent personality disorder")
                    || trait.contains("conscientiousness") || trait.contains("cognitive") || trait.contains("anxious")
                    || trait.contains("bulimia nervosa")
            ) { // psych phenotypes
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("net100") || trait.startsWith("net25") || trait.startsWith("a2009s") || trait.startsWith("nodeamps25")) { // weird codified phenotypes
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("volume ") || trait.startsWith("idp dmri") || trait.startsWith("idp t2") || trait.startsWith("idp swi")
                    || trait.startsWith("dkatlas")) { // MRI brain volume
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.equals("tau")) { // MRI brain volume
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("nerve growth factor")) { // MRI brain volume
//                System.err.println(trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("body mass index") || trait.contains("bmi") || trait.contains("waist-to-hip ratio") || trait.contains("waist-hip ratio")) {
//                System.err.println("BMI: " + trait);
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("verbal")) {
//                System.err.println("memory: " + trait);
            } else if (trait.contains("white matter")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("viral meningitis") || trait.contains("viral encephalitis")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("subarachnoid")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("stroke")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("speech")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("sleeplessness") || trait.contains("insomnia") || trait.contains("sleep")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("progressive supranuclear palsy")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("nervous system")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("suicide") || trait.contains("post-traumatic stress disorder")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("opioid") || trait.contains("addiction") || trait.contains("nicotine dependence")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("intracranial")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("neuro") || trait.contains("neura") || trait.contains("neure") || trait.contains("neuri")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("neuropeptide") || trait.contains("neuropilin") || trait.contains("neurotensin")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("neuropeptide")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("narcolepsy")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("neo-")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("motor disorders") || trait.startsWith("motor neuron")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("migraine")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("mental retardation")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("huntington's disease")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("handedness")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.startsWith("glioma") || trait.contains("glioblastoma")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("neuropathy")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("creutzfeldt-jakob disease")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("headache")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("cerebrovascular")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            } else if (trait.contains("cerebrospinal")) {
                traitToAnnotation.put(trait, "brain");
                braintraits++;
            }

            // other traits
            if (!traitToAnnotation.containsKey(trait)) {
                if (Character.isDigit(trait.charAt(0)) || trait.startsWith("gluc") || trait.startsWith("glut") || trait.startsWith("gly")) {
                    traitToAnnotation.put(trait, "protein or metabolite");
                } else if (trait.startsWith("trna") || trait.contains(" trna ") || trait.startsWith("rna") || trait.contains(" rna ")
                        || trait.startsWith("rrna") || trait.contains(" rrna ") || trait.startsWith("mrna") || trait.contains(" mrna ")) {
//                    System.err.println(trait);
                    traitToAnnotation.put(trait, "rna");
                    eqtltraits++;
                } else if (trait.startsWith("protein") || trait.contains("globulin")) {
                    traitToAnnotation.put(trait, "protein");
                } else if (trait.startsWith("ensg")) {
                    traitToAnnotation.put(trait, "eqtl");
                    eqtltraits++;
                } else if (trait.contains("serum")) {
                    traitToAnnotation.put(trait, "serum");
                    eqtltraits++;
                } else if (trait.contains("arthritis")) {
                    traitToAnnotation.put(trait, "arthritis");
                    eqtltraits++;
                } else if (trait.contains("serine") || trait.contains("threonine")) {
                    traitToAnnotation.put(trait, "aminoacids");
                    eqtltraits++;
                } else if (trait.contains("thyroid")) {
                    traitToAnnotation.put(trait, "thyroid");
                    eqtltraits++;
                } else if (trait.startsWith("interleukin") || trait.contains("cd4") || trait.contains("cd8")
                        || trait.contains("hla") || trait.contains("interferon") || trait.contains("white blood cell") ||
                        trait.contains("leukocyte") || trait.contains("platelet") || trait.contains("plasma")) {
                    traitToAnnotation.put(trait, "immune");
                    braineqtltraits++;
                } else if (trait.contains("lipids")) {
                    traitToAnnotation.put(trait, "lipids");
                    braineqtltraits++;
                } else if (trait.startsWith("kallikrein") || trait.startsWith("zinc finger")) {
                    traitToAnnotation.put(trait, "pqtl");
                    braineqtltraits++;
                } else if (trait.contains("job soc coding") || trait.contains("job coding") ||
                        trait.contains("leisure/social activities") || trait.contains("gap coding")) {
                    traitToAnnotation.put(trait, "questionnaire");
                } else if (trait.startsWith("x-")) {
                    traitToAnnotation.put(trait, "pqtl");
                } else if (trait.contains("atherosclerosis") || trait.contains("cardio") || trait.contains("heart") ||
                        trait.contains("pulse") || trait.contains("blood pressure") || trait.contains("hypertension")
                        || trait.contains("pericard") || trait.contains("vascular") || trait.contains("myocard")) {
                    traitToAnnotation.put(trait, "cardiovascular");
                } else if (trait.contains("pulmonary")) {
                    traitToAnnotation.put(trait, "pulmonary");
                } else if (trait.contains("neoplasm") || trait.contains("carcinoma") || trait.contains("cancer")) {
                    traitToAnnotation.put(trait, "cancer");
                } else if (trait.contains("vitamin")) {
                    traitToAnnotation.put(trait, "vitamin");
                } else if (trait.contains("viral")) {
                    traitToAnnotation.put(trait, "viral");
                } else if (trait.startsWith("pct")) {
                    traitToAnnotation.put(trait, "questionnaire");
                } else if (trait.startsWith("other") || trait.startsWith("operative procedures") || trait.startsWith("operation")) {
                    traitToAnnotation.put(trait, "other");
                } else if (trait.contains("osteo")) {
                    traitToAnnotation.put(trait, "bone");
                }
            }
            if (!traitToAnnotation.containsKey(trait)) {

//                System.out.println(trait);

//                String[] telems = trait.split(" ");
//                for (String word : telems) {
//                    if (word.length() > 5) {
//                        Integer ctr = wordcount.get(word);
//                        if (ctr == null) {
//                            ctr = 0;
//                        }
//                        ctr++;
//
//                        wordcount.put(word, ctr);
//                    }
//                }
                ctrnotincluded++;
            }
        }

//        ArrayList<Pair<String, Integer>> pairs = new ArrayList<Pair<String, Integer>>();
//        for (String word : wordcount.keySet()) {
//            pairs.add(new Pair<>(word, wordcount.get(word), Pair.SORTBY.RIGHT));
//        }
//        Collections.sort(pairs);
//
//        for (Pair<String, Integer> word : pairs) {
//            System.out.println(word.getLeft() + "\t" + word.getRight());
//        }

        System.out.println(ctrnotincluded + " traits without annotation");
        System.out.println(braintraits + " brain traits");
        System.out.println(braineqtltraits + " brain eqtl traits");
        System.out.println(eqtltraits + " eqtl traits");
        System.out.println();
        System.out.println("---------------------------------------");
        for (String trait : allTraits) {
            trait = trait.toLowerCase();
            String annot = traitToAnnotation.get(trait);
            if (annot != null && annot.equals("brain")) {
                System.out.println(trait);
            }
        }

        tf.open();
        TextFile tfo = new TextFile(output, TextFile.W);
        tfo.writeln(tf.readLine());
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[0];
            String trait = elems[3];
            trait = trait.toLowerCase();
            String annot = traitToAnnotation.get(trait);
            if (annot != null && annot.equals("brain")) {
                tfo.writeln(Strings.concat(elems, Strings.tab));
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        tfo.close();

        tf = new TextFile(assoc, TextFile.R);
        tfo = new TextFile(assocout, TextFile.W);
        TextFile tfo2 = new TextFile(tophitsout, TextFile.W);
        tfo.writeln(tf.readLine());
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String id = elems[0];
            String trait = traitIdToTrait.get(id);
            trait = trait.toLowerCase();
            String annot = traitToAnnotation.get(trait);
            if (annot != null && annot.equals("brain")) {
                tfo.writeln(Strings.concat(elems, Strings.tab));
                tfo2.writeln(elems[1]);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();
        tfo.close();
        tfo2.close();
    }
}
