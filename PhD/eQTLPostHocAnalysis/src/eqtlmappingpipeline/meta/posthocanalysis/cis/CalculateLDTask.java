/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package eqtlmappingpipeline.meta.posthocanalysis.cis;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import umcg.genetica.containers.Pair;
import umcg.genetica.containers.SortableSNP;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author harmjan
 */
public class CalculateLDTask implements Callable<Pair<String, ArrayList<String>>> {

    private String query = null;
    private SortableSNP[] sortedSNPs = null;
    private TriTyperGenotypeData d = null;

    public CalculateLDTask(TriTyperGenotypeData d, String querySNP, SortableSNP[] sortedSNPs) {
        this.d = d;
        this.query = querySNP;
        this.sortedSNPs = sortedSNPs;
    }

    @Override
    public Pair<String, ArrayList<String>> call() {

        try {
            int querySNPPos = -1;
            for (int refsnp = 0; refsnp < sortedSNPs.length; refsnp++) {
                if (sortedSNPs[refsnp].name.equals(query)) {
                    // found it :)
                    querySNPPos = refsnp;
                    break;
                }
            }



            if (querySNPPos > -1) {
                SNPLoader loader = d.createSNPLoader();
                DetermineLD ldcalc = new DetermineLD();
                // we have a starting point.. 
                ArrayList<String> proxies = new ArrayList<String>();
                byte querySNPChr = sortedSNPs[querySNPPos].chr;
                int querySNPChrPos = sortedSNPs[querySNPPos].chrpos;

                SNP qSNPObj = d.getSNPObject(sortedSNPs[querySNPPos].id);
                loader.loadGenotypes(qSNPObj);


                // go 250kb downstream and upstream
                for (int refsnp = querySNPPos - 1; refsnp > -1; refsnp--) {
                    SortableSNP refSNPObj = sortedSNPs[refsnp];
                    if (refSNPObj.chr == querySNPChr) {
                        int distance = Math.abs(refSNPObj.chrpos - querySNPChrPos);
                        if (distance <= 250000) {
                            // calculate LD
                            SNP rSNPObj = d.getSNPObject(sortedSNPs[refsnp].id);
                            loader.loadGenotypes(rSNPObj);
                            double r2 = ldcalc.getRSquared(qSNPObj, rSNPObj, d, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                            if (r2 >= 1d) {
                                proxies.add(rSNPObj.getName());
                            }
                            rSNPObj.clearGenotypes();
                        }
                    }
                }

                for (int refsnp = querySNPPos + 1; refsnp < sortedSNPs.length; refsnp++) {
                    SortableSNP refSNPObj = sortedSNPs[refsnp];
                    if (refSNPObj.chr == querySNPChr) {
                        int distance = Math.abs(refSNPObj.chrpos - querySNPChrPos);
                        if (distance <= 250000) {
                            // calculate LD
                            SNP rSNPObj = d.getSNPObject(sortedSNPs[refsnp].id);
                            loader.loadGenotypes(rSNPObj);
                            double r2 = ldcalc.getRSquared(qSNPObj, rSNPObj, d, DetermineLD.RETURN_R_SQUARED, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
                            if (r2 >= 1d) {
                                proxies.add(rSNPObj.getName());
                            }
                            rSNPObj.clearGenotypes();
                        }
                    }
                }
                qSNPObj.clearGenotypes();

                

                loader.close();
                // done finding proxies, now print the f**kers..
//            tfout.writeln(query);
                if (proxies.size() > 0) {
//                snpsWithProxies++;
//                for (String proxy : proxies) {
//                    tfout.writeln(proxy);
//                }
                    return new Pair<String, ArrayList<String>>(query, proxies);
                } else {
                    return new Pair<String, ArrayList<String>>(query, null);
                }



                


            } else {
                // no proxies because the SNP is not in the reference.. print an error though
                return new Pair<String, ArrayList<String>>(query, null);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return new Pair<String, ArrayList<String>>(query, null);
    }
}
