package nl.harmjanwestra.playground.biogen.freeze2dot1.convergence;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;
import org.graphstream.stream.file.FileSink;
import org.graphstream.stream.file.FileSinkFactory;
import org.graphstream.stream.file.FileSinkSVG2;
import org.graphstream.ui.layout.Layout;
import org.graphstream.ui.layout.Layouts;
import org.graphstream.ui.layout.springbox.implementations.SpringBox;
import org.graphstream.ui.view.Viewer;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Objects;

public class DetermineConvergentEffects {


    public static void main(String[] args) {

        String cisFile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-15-assoc\\cortex-cis-EUR\\eQTLsFDR0.05-diseaseSNPs-iter1-4.txt.gz";
        String transFile = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-15-assoc\\cortex-trans-EUR\\eQTLsFDR0.05-diseaseAnnotation.txt.gz";

        String gwaslist = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-15-assoc\\2020-05-03-gwaslist-wgwascatalog-wALS.txt.gz";
        String query = "schizo";

//        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-15-assoc\\cortex-EUR-convergence\\schizophrenia.svg";
        String output = "D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2020-01-Freeze2dot1\\2020-05-15-assoc\\cortex-EUR-convergence\\allphenotypes.txt";
        DetermineConvergentEffects d = new DetermineConvergentEffects();
        try {
//            d.run(cisFile, transFile, gwaslist, query, false, output);
            d.toText(cisFile, transFile, gwaslist, null, false, output);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*

    snp --> cis
    othersnp --> trans

     */

    class EQTLSNP {
        HashSet<String> cisGenes = new HashSet<>();
        HashSet<String> transGenes = new HashSet<>();
    }

    class EQTLGene {
        HashSet<String> cisSNPs = new HashSet<>();
        HashSet<String> transSNPs = new HashSet<>();
    }

    class NODE {
        String name = null;
        NODETYPE type = NODETYPE.GENE;
        HashSet<EDGE> edges = new HashSet<>();

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            NODE node = (NODE) o;
            return Objects.equals(name, node.name);
        }

        @Override
        public int hashCode() {
            return Objects.hash(name);
        }

        public void addEdge(NODE genenode, NODE snpnode, EDGETYPE edgetype) {
            EDGE edge = new EDGE();
            edge.node1 = genenode;
            edge.node2 = snpnode;
            edge.type = edgetype;
            edges.add(edge);
        }

        public int countSharedEdges(NODE cisnode2) {
            int shared = 0;
            for (EDGE e : edges) {
                for (EDGE e2 : cisnode2.edges) {
                    if ((e.node1.equals(e2.node1) && e.node2.equals(e2.node2)) || (e.node1.equals(e2.node2) && e.node2.equals(e2.node1))) {
                        shared++;
                    }
                }
            }
            return shared;
        }

        public NODE edgelessCopy() {
            NODE copy = new NODE();
            copy.name = name;
            copy.type = type;
            return copy;
        }

        public void addEdge(EDGE e) {
            edges.add(e);
        }
    }

    enum EDGETYPE {
        CIS,
        TRANS
    }

    enum NODETYPE {
        SNP,
        GENE
    }

    class EDGE {
        EDGETYPE type = EDGETYPE.CIS;
        NODE node1 = null;
        NODE node2 = null;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            EDGE edge = (EDGE) o;


            return (Objects.equals(node1, edge.node1) &&
                    Objects.equals(node2, edge.node2)) || (Objects.equals(node1, edge.node2) &&
                    Objects.equals(node2, edge.node1));
        }

        @Override
        public int hashCode() {
            return Objects.hash(node1, node2);
        }

        public String getName() {
            return node1.name + "-" + node2.name;
        }
    }

    class NETWORK {
        ArrayList<NODE> nodes = new ArrayList<>();

        public void addNode(NODE node) {
            nodes.add(node);
        }
    }


    public void toText(String cisfile, String transfile, String gwaslistfile, String query, boolean prune, String output) throws IOException {
        if (query != null) {
            HashSet<String> allowedGWASIds = null;
            if (gwaslistfile != null) {
                allowedGWASIds = new HashSet<>();
                TextFile t = new TextFile(gwaslistfile, TextFile.R);
                t.readLine();
                String[] elems = t.readLineElems(TextFile.tab);

                while (elems != null) {
                    if (query == null || elems[2].contains(query)) {
                        allowedGWASIds.add(elems[0]);
                    }
                    elems = t.readLineElems(TextFile.tab);
                }
                t.close();
            }
        } else {
            HashSet<String> allGWASIds = new HashSet<String>();
            TextFile t = new TextFile(gwaslistfile, TextFile.R);
            t.readLine();
            String[] elems = t.readLineElems(TextFile.tab);

            HashMap<String, String> idToTrait = new HashMap<String, String>();
            while (elems != null) {
                allGWASIds.add(elems[0]);
                String trait = elems[3];
                if (elems[0].startsWith("eqtl-a")) {
                    trait = "eqtlgen-" + elems[3];
                }
                idToTrait.put(elems[0], trait);
                elems = t.readLineElems(TextFile.tab);
            }
            t.close();

            TextFile out = new TextFile(output, TextFile.W);
            out.write("GWAS\tTrait\tnCis\tnTrans\tCisGenes\tTransGenes\tCisSNPs\tTransSNPs");
            for (String id : allGWASIds) {
                String trait = idToTrait.get(id);
                System.out.println(id + "\t" + trait);
                HashSet<String> allowedGWASIDs = new HashSet<>();
                allowedGWASIDs.add(id);
                ArrayList<NETWORK> subnets = getSubnets(allowedGWASIDs, cisfile, transfile, prune);

                for (NETWORK net : subnets) {
                    ArrayList<String> cisgenes = new ArrayList<>();
                    ArrayList<String> transgenes = new ArrayList<>();
                    ArrayList<String> cissnps = new ArrayList<>();
                    ArrayList<String> transsnps = new ArrayList<>();

                    for (NODE n : net.nodes) {
                        if (n.type.equals(NODETYPE.GENE)) {
                            for (EDGE e : n.edges) {
                                if (e.type.equals(EDGETYPE.CIS)) {
                                    cisgenes.add(n.name);
                                } else {
                                    transgenes.add(n.name);
                                }
                            }
                        } else {
                            for (EDGE e : n.edges) {
                                if (e.type.equals(EDGETYPE.CIS)) {
                                    cissnps.add(n.name);
                                } else {
                                    transsnps.add(n.name);
                                }
                            }
                        }
                    }

                    out.writeln(id + "\t" + trait + "\t" +
                            cisgenes.size() + "\t" + transgenes.size() +
                            Strings.concat(cisgenes, Strings.semicolon) + "\t" +
                            Strings.concat(transgenes, Strings.semicolon) + "\t" +
                            Strings.concat(cissnps, Strings.semicolon) + "\t" +
                            Strings.concat(transsnps, Strings.semicolon) + "\t"
                    );

                }
            }
            out.close();
        }
    }

    public ArrayList<NETWORK> getSubnets(HashSet<String> allowedGWASIds, String cisfile, String transfile, boolean prune) throws IOException {
        HashMap<String, EQTLGene> genes = new HashMap<>();
        HashMap<String, EQTLSNP> snps = new HashMap<>();

        TextFile tf = new TextFile(cisfile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String gene = elems[2];
            String gwasidstr = elems[5];
            String[] gwasids = gwasidstr.split(";");
            boolean ok = false;
            for (String gwasid : gwasids) {
                if (allowedGWASIds == null || allowedGWASIds.contains(gwasid)) {
                    ok = true;
                }
            }
            if (ok) {
                EQTLGene g = genes.get(gene);
                if (g == null) {
                    g = new EQTLGene();
                }
                g.cisSNPs.add(snp);
                genes.put(gene, g);
                EQTLSNP s = snps.get(snp);
                if (s == null) {
                    s = new EQTLSNP();
                }
                s.cisGenes.add(gene);
                snps.put(snp, s);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        tf = new

                TextFile(transfile, TextFile.R);
        tf.readLine();
        elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            String snp = elems[1];
            String gene = elems[2];
            String gwasidstr = elems[5];
            String[] gwasids = gwasidstr.split(";");
            boolean ok = false;
            for (String gwasid : gwasids) {
                if (allowedGWASIds == null || allowedGWASIds.contains(gwasid)) {
                    ok = true;
                }
            }
            if (ok) {
                EQTLGene g = genes.get(gene);
                if (g == null) {
                    g = new EQTLGene();
                }
                g.transSNPs.add(snp);
                genes.put(gene, g);
                EQTLSNP s = snps.get(snp);
                if (s == null) {
                    s = new EQTLSNP();
                }
                s.transGenes.add(gene);
                snps.put(snp, s);
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        // build network

        HashMap<String, NODE> snpToNode = new HashMap<>();
        HashMap<String, NODE> geneToNode = new HashMap<>();

        NETWORK network = new NETWORK();
        for (
                String gene : genes.keySet()) {
            EQTLGene g = genes.get(gene);

            NODE genenode = geneToNode.get(gene);

            if (genenode == null) {
                genenode = new NODE();
                genenode.name = gene;
                genenode.type = NODETYPE.GENE;
                geneToNode.put(gene, genenode);
                network.addNode(genenode);
            }

            for (String snp : g.cisSNPs) {
                NODE snpnode = snpToNode.get(snp);
                if (snpnode == null) {
                    snpnode = new NODE();
                    snpnode.name = snp;
                    snpnode.type = NODETYPE.SNP;
                    snpToNode.put(snp, snpnode);
                    network.addNode(snpnode);
                }

                genenode.addEdge(genenode, snpnode, EDGETYPE.CIS);
                snpnode.addEdge(snpnode, genenode, EDGETYPE.CIS);
            }

            for (String snp : g.transSNPs) {
                NODE snpnode = snpToNode.get(snp);
                if (snpnode == null) {
                    snpnode = new NODE();
                    snpnode.name = snp;
                    snpnode.type = NODETYPE.SNP;
                    snpToNode.put(snp, snpnode);
                    network.addNode(snpnode);
                }

                genenode.addEdge(genenode, snpnode, EDGETYPE.TRANS);
                snpnode.addEdge(snpnode, genenode, EDGETYPE.TRANS);
            }
        }

        System.out.println(geneToNode.size() + " gene nodes");
        System.out.println(snpToNode.size() + " SNP nodes");
        System.out.println("Full network has: " + network.nodes.size() + " nodes.");

        // make subnets
        HashSet<NODE> visitedNodes = new HashSet<NODE>();
        ArrayList<NETWORK> subnets = new ArrayList<>();

        int nrVisited = 0;
        for (NODE n : network.nodes) {
            if (!visitedNodes.contains(n)) {
                // walk the initial set of edges
                HashSet<NODE> connectedNodes = new HashSet<>();
                connectedNodes.add(n); // add the 'seed' node
                HashSet<NODE> newlyAddedNodes = new HashSet<>();
                for (EDGE e : n.edges) {
                    if (!connectedNodes.contains(e.node1)) {
                        connectedNodes.add(e.node1);
                        newlyAddedNodes.add(e.node1);
                    }
                    if (!connectedNodes.contains(e.node2)) {
                        connectedNodes.add(e.node2);
                        newlyAddedNodes.add(e.node2);
                    }
                }

                // and iterate until no new nodes are added
                while (newlyAddedNodes.size() > 0) {
                    HashSet<NODE> curNewlyAddedNodes = new HashSet<>();

                    // iterate only newly added nodes
                    for (NODE n2 : newlyAddedNodes) {
                        for (EDGE e : n2.edges) {
                            // don't consider nodes already connected
                            if (!connectedNodes.contains(e.node1)) {
                                connectedNodes.add(e.node1);
                                curNewlyAddedNodes.add(e.node1);
                            }
                            if (!connectedNodes.contains(e.node2)) {
                                connectedNodes.add(e.node2);
                                curNewlyAddedNodes.add(e.node2);
                            }
                        }
                    }
                    newlyAddedNodes = curNewlyAddedNodes;
                }
                NETWORK subnet = new NETWORK();
                subnet.nodes.addAll(connectedNodes);
                subnets.add(subnet);
                visitedNodes.addAll(connectedNodes);
            }
            System.out.println(visitedNodes.size() + " out of " + network.nodes.size() + " visited.");
        }

        System.out.println("There are " + subnets.size() + " unconnected subnets.");


        // prune the subnets to exclude multiple cisfx
        ArrayList<NETWORK> prunedsubnets = new ArrayList<>();
        if (prune) {
            for (
                    NETWORK subnet : subnets) {


                ArrayList<NODE> geneNodes = new ArrayList<>();
                for (NODE n : subnet.nodes) {
                    if (n.type.equals(NODETYPE.GENE)) {
                        geneNodes.add(n);
                    }
                }

                HashMap<NODE, HashSet<NODE>> synonyms = new HashMap<>();
                HashSet<NODE> nodesThatAreSynonyms = new HashSet<>();

                for (NODE n : geneNodes) {
                    // get a list of cis nodes
                    HashSet<NODE> cisnodes = new HashSet<>();
                    for (EDGE e : n.edges) {
                        if (e.type.equals(EDGETYPE.CIS)) {
                            if (e.node1.type.equals(NODETYPE.SNP)) {
                                cisnodes.add(e.node1);
                            }
                            if (e.node2.type.equals(NODETYPE.SNP)) {
                                cisnodes.add(e.node2);
                            }
                        }
                    }

                    ArrayList<NODE> cisNodeList = new ArrayList<>();
                    cisNodeList.addAll(cisnodes);

                    // determine nodes with equal edges
                    for (int i = 0; i < cisNodeList.size(); i++) {
                        NODE cisnode1 = cisNodeList.get(i);
                        if (!nodesThatAreSynonyms.contains(cisnode1)) {
                            HashSet<NODE> syn = new HashSet<>();
                            for (int j = i + 1; j < cisNodeList.size(); j++) {
                                NODE cisnode2 = cisNodeList.get(i);
                                if (!nodesThatAreSynonyms.contains(cisnode2)) {
                                    int nrshared = cisnode1.countSharedEdges(cisnode2);
                                    if (nrshared == cisnode1.edges.size()) {
                                        // identical list of edges for this node
                                        syn.add(cisnode2);
                                    }
                                }
                            }
                            nodesThatAreSynonyms.addAll(syn);
                            synonyms.put(cisnode1, syn);
                        }
                    }
                }
                // make new network using pruned nodes and edges
                NETWORK prunedSubnet = new NETWORK();
                for (NODE n : subnet.nodes) {
                    // exclude synonymous nodes
                    if (!nodesThatAreSynonyms.contains(n)) {
                        NODE prunedCopy = n.edgelessCopy();
                        // add edges that still remain
                        for (EDGE e : n.edges) {
                            if (!nodesThatAreSynonyms.contains(e.node1) && !nodesThatAreSynonyms.contains(e.node2)) {
                                prunedCopy.addEdge(e);
                            }
                        }
                        prunedSubnet.addNode(prunedCopy);
                    }
                }
                prunedsubnets.add(prunedSubnet);
            }
        } else {
            prunedsubnets = subnets;
        }
        return prunedsubnets;
    }

//    public void run(String cisfile, String transfile, String gwaslistfile, String query, boolean prune, String output) throws IOException {
//
//
//
//
//
//        // find graph with maximal number of gene nodes
//        int max = 0;
//        NETWORK chosen = null;
//        for (NETWORK subnet : prunedsubnets) {
//            int ctr = 0;
//            for (NODE n : subnet.nodes) {
//                if (n.type.equals(NODETYPE.GENE)) {
//                    ctr++;
//                }
//            }
//            if (ctr > max) {
//                max = ctr;
//                chosen = subnet;
//            }
//        }
//
//        System.out.println("Max subnet has " + max + " nodes.");
//
//        // draw subnets
////        for (NETWORK subnet : subnets) {
////        System.setProperty("gs.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
//        Layout layout = Layouts.newLayoutAlgorithm();
//        Graph graph = Graphs.synchronizedGraph(new MultiGraph("Max subnet"));
//
//        graph.addSink(layout);
//        layout.addAttributeSink(graph);
//        graph.addAttribute("ui.quality", 0);
//        graph.addAttribute("ui.antialias");
//        graph.addAttribute("ui.stylesheet",
//                "node.gene { fill-color: #0476D9; text-alignment: center; text-background-mode: rounded-box; text-background-color: #0476D9; text-mode: normal; text-color: white; } " +
//                        "node.snp { fill-color: #93bf34; text-alignment: center; size-mode: fit; text-background-mode: rounded-box; text-background-color: #93bf34; text-mode: normal; text-color: white; }" +
//                        "edge.cis { fill-color: grey; }" +
//                        "edge.trans {fill-color: #f2b705;}");
//
////        for(NETWORK net: prunedsubnets){
////
////        }
//        for (NODE node : chosen.nodes) {
//            Node n = graph.addNode(node.name);
//            if (node.type.equals(NODETYPE.GENE)) {
//                n.setAttribute("ui.class", "gene");
//                n.setAttribute("layout.weight", 1);
//            } else {
//                n.setAttribute("ui.class", "snp");
//                n.setAttribute("layout.weight", 1);
//            }
//            n.setAttribute("ui.label", node.name);
//        }
//        for (NODE node : chosen.nodes) {
//            for (EDGE e : node.edges) {
//                try {
//
//                    Edge edge = graph.addEdge(e.getName(), e.node1.name, e.node2.name);
//                    if (e.type.equals(EDGETYPE.TRANS)) {
//                        edge.setAttribute("layout.weight", 1);
//                        edge.setAttribute("ui.class", "trans");
//                    } else {
//                        edge.setAttribute("layout.weight", 1);
//                        edge.setAttribute("ui.class", "cis");
//                    }
//
//                } catch (EdgeRejectedException ex) {
//
//                }
//            }
//        }
//
//
////
//
//
//        /*
//        gene: blue
//        snp: green
//        cis: grey
//        trans: yellow
//         */
//
//
//        // iterate the compute() method a number of times
////        while (layout.getStabilization() < 1) {
////            layout.compute();
////        }
//
//
//        layout.shake();
//        while (layout.getStabilization() < 1) {
//            layout.compute();
//        }
//        FileSink out = new FileSinkSVG2();
//
//        out.writeAll(graph, output);
//        out.flush();
//
//        Viewer viewer = graph.display();
//        viewer.disableAutoLayout();
//
//
////        System.out.println("Graph has " + graph.getNodeCount() + " nodes");
////
////        // Let the layout work ...
////
////// Do some work ...
////
////        }
//
//    }

}
