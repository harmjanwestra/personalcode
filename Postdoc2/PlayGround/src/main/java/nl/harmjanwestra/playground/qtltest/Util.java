package nl.harmjanwestra.playground.qtltest;

import java.util.ArrayList;
import java.util.HashMap;

public class Util {

    public static HashMap<String, Integer> hash(ArrayList<String> list) {
        HashMap<String, Integer> hash = new HashMap<>();
        for (int i = 0; i < list.size(); i++) {
            hash.put(list.get(i), i);
        }
        return hash;
    }
}
