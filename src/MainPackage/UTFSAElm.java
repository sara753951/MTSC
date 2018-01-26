package MainPackage;

import javafx.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class UTFSAElm {
    HashSet<Integer> objSet = new HashSet<Integer>();
    HashMap<Integer, Integer> objCount = new HashMap<Integer, Integer>();
    ArrayList<Pair<Integer, Integer>> tfsaID = new ArrayList<Pair<Integer, Integer>>();
    int lenOfSubseqInData;

    public UTFSAElm() {
    }

    public UTFSAElm(HashSet<Integer> objSet, HashMap<Integer, Integer> objCount, ArrayList<Pair<Integer, Integer>> tfsaID, int lenOfSubseqInData) {
        this.objSet = objSet;
        this.objCount = objCount;
        this.tfsaID = tfsaID;
        this.lenOfSubseqInData = lenOfSubseqInData;
    }

    public HashSet<Integer> getObjSet() {
        return objSet;
    }

    public void setObjSet(HashSet<Integer> objSet) {
        this.objSet = objSet;
    }

    public HashMap<Integer, Integer> getObjCount() {
        return objCount;
    }

    public void setObjCount(HashMap<Integer, Integer> objCount) {
        this.objCount = objCount;
    }

    public ArrayList<Pair<Integer, Integer>> getTfsaID() {
        return tfsaID;
    }

    public void setTfsaID(ArrayList<Pair<Integer, Integer>> tfsaID) {
        this.tfsaID = tfsaID;
    }

    public int getLenOfSubseqInData() {
        return lenOfSubseqInData;
    }

    public void setLenOfSubseqInData(int lenOfSubseqInData) {
        this.lenOfSubseqInData = lenOfSubseqInData;
    }
}
