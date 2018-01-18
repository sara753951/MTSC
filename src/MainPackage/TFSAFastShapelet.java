package MainPackage;

import java.util.ArrayList;
import java.util.HashSet;

public class TFSAFastShapelet {

    public Shapelet ShapeletDiscovery (ArrayList<ArrayList<Double>> data) {

        Shapelet sh = new Shapelet();
        HashSet<Integer> KeyPoints = new HashSet<Integer>();

        for(int i=0;i<data.size();i++)
        {
            ArrayList<TFSA_word> TFSAListofT;
            KeyPoints.clear();
            int len=data.get(0).size();
            ArrayList<Double> T = data.get(i);
            AS3(T,Math.tan(30.0),len/min_len);
            TFSAListofT=Symbolization(T);
            DataofTFSA.push_back(TFSAListofT);

        }

        return sh;
    }

    public AS(ArrayList<Double> ){

    }
}
