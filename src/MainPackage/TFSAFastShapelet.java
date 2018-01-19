package MainPackage;

import java.util.ArrayList;

public class TFSAFastShapelet {

    private static final double PI = 3.14;
    ArrayList<Integer> KeyPoints = new ArrayList<Integer>();

    public Shapelet ShapeletDiscovery(ArrayList<ArrayList<Double>> data) {

        int min_len = 1;
        Shapelet sh = new Shapelet();

        ArrayList<ArrayList<TFSAWord>> DataofTFSA = new ArrayList<ArrayList<TFSAWord>>();

        for (int i = 0; i < data.size(); i++) {
            ArrayList<TFSAWord> TFSAListofT;
            KeyPoints.clear();
            int len = data.get(0).size();
            ArrayList<Double> T = data.get(i);
            AS(T, Math.tan(30.0), len / min_len);
            TFSAListofT = Symbolization(T);
            DataofTFSA.add(TFSAListofT);
        }



        return sh;
    }

    public void AS(ArrayList<Double> S, double angle_tol, int k) {
        int l = (int) Math.ceil(S.size() / k);
        //ubseq_type S=T;
        ArrayList<Double> S1 = (ArrayList<Double>) S.subList(0, 0 + l);
        ArrayList<Double> S2;
        double K1 = Slope(S1);
        double K2;
        int i = 1;
        while (i < (S.size() - l + 1)) {
            S2 = (ArrayList<Double>) S.subList(i, i + l);
            K2 = Slope(S2);
            if (Math.abs(K1 - K2) > angle_tol) {
                KeyPoints.add(l + i - 2);
                i = l + i - 2;
                S1 = (ArrayList<Double>) S.subList(i, i + l);
                K1 = Slope(S1);
                i++;
            } else {
                K1 = (K1 + K1) / 2;
                i++;
            }
        }
    }

    public double Slope(ArrayList<Double> S) {
        double vt = 0, v = 0, t = 0, v2 = 0;
        //int len=S.size();

        for (int i = 0; i < S.size(); i++) {
            vt += S.get(i) * i;
            v += S.get(i);
            t += i;
            v2 += Math.pow(S.get(i), 2);
        }
        double a = (S.size() * vt - v * t) / (S.size() * v2 - 2 * v);

        //double a=(S[len-1]-S[0])/(len-1);
        return a;
    }

    public ArrayList<TFSAWord> Symbolization(ArrayList<Double> T) {
        ArrayList<TFSAWord> TFSAList = new ArrayList<TFSAWord>();

        ArrayList<Double> subseq;
        int it_previous;
        if (KeyPoints.size() == 0) {
            KeyPoints.add((int) Math.floor((double) T.size() / 12));
            KeyPoints.add((int) Math.floor((double) (2 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (3 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (4 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (5 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (6 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (7 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (8 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (9 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (10 * T.size() / 12)));
            KeyPoints.add((int) Math.floor((double) (11 * T.size() / 12)));
        }
        subseq = (ArrayList<Double>) T.subList(0, KeyPoints.get(0));
        if (subseq.size() > 1)
            TFSAList.add(Symbolize(subseq, 0));

        int i = 1;
        for (; i < KeyPoints.size(); i++) {

            subseq = (ArrayList<Double>) T.subList(KeyPoints.get(i - 1), KeyPoints.get(i));
            TFSAList.add(Symbolize(subseq, KeyPoints.get(i - 1)));
        }
        subseq = (ArrayList<Double>) T.subList(KeyPoints.get(i - 1), T.size() - 1);
        TFSAList.add(Symbolize(subseq, KeyPoints.get(i - 1)));

        return TFSAList;
    }

    public TFSAWord Symbolize(ArrayList<Double> T, int posOfBegin) {
        TFSAWord tfsaWord = new TFSAWord();

        double a = Slope(T);

        if (a >= 0) {
            tfsaWord.setDegree((int) Math.ceil(Math.atan(a) * 180 / (10 * PI)));
        } else {
            tfsaWord.setDegree((int) Math.floor(Math.atan(a) * 180 / (10 * PI)));
        }
        double lastpoint = T.get(T.size() - 1);
        if (lastpoint < -0.67)
            tfsaWord.setLastPoint('A');
        else if (lastpoint < 0)
            tfsaWord.setLastPoint('B');
        else if (lastpoint < 0.67)
            tfsaWord.setLastPoint('C');
        else
            tfsaWord.setLastPoint('D');
        tfsaWord.setLen(T.size());
        tfsaWord.setPosOfBeginInData(posOfBegin);

        return tfsaWord;
    }
}
