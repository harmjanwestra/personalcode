/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 *
 * @author harmjan
 */
public class TestHash {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        Integer a = 1;
        Integer b = 2;
        Integer c = 3;
        Integer d = 4200;
        Integer e = 4200;

        Set<Integer> h1 = new TreeSet<Integer>();
        Set<Integer> h2 = new TreeSet<Integer>();

        System.out.println(h1.hashCode());
        System.out.println(h2.hashCode());

        h1.add(a);
        h2.add(a);
        h2.add(0);

        for (int i = -4; i < 6; i++) {
            h1.add(i);
        }
        
//        String bla = h1.toString();

        System.out.println(h1);
        System.out.println(h1.hashCode());
        System.out.println(h2.hashCode());

        h1.add(b);
        h1.add(c);
        h2.add(a);
        h2.add(c);
        h1.add(d);
        h2.add(e);

        System.out.println(h1.hashCode());
        System.out.println(h2.hashCode());


    }
}
