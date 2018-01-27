package MainPackage;

import javafx.util.Pair;
import org.omg.CORBA.INTERNAL;

import java.util.*;

public class TFSAFastShapelet {


    private static final double PI = 3.14;
    private static final int SH_MIN_LEN = 5;
    ArrayList<Integer> KeyPoints = new ArrayList<Integer>();
    ArrayList<ArrayList<TFSAWord>> DataofTFSA = new ArrayList<ArrayList<TFSAWord>>();
    ArrayList<ArrayList<TFSAWord>> TFSASet = new ArrayList<ArrayList<TFSAWord>>();
    ArrayList<ArrayList<TFSAWord>> TFSANewSet = new ArrayList<ArrayList<TFSAWord>>();
    LinkedHashMap<Integer, UTFSAElm> UTFSAMap = new LinkedHashMap<Integer, UTFSAElm>();
    ArrayList<Integer> Label = new ArrayList<Integer>();
    int num_class;
    int num_obj;
    int subseq_len;
    int[] Class_Freq = new int[num_class];
    ArrayList<Pair<Integer, Double>> TFSAScoreList = new ArrayList<Pair<Integer, Double>>();
    double[][] Data;
    double class_entropy;


    public Shapelet ShapeletDiscovery(ArrayList<ArrayList<Double>> data) {

        this.Data = readData(data);
        int min_len = 1, max_len = 1;
        int step = 1;
        Shapelet sh = new Shapelet();
        int subseq_len = 0;
        int sax_max_len = 15;
        int card = 4;
        int R = 10;
        double percent_mask = 0.25;
        int top_k = 10;
        ArrayList<Shapelet> allLengthBestShapelet = new ArrayList<Shapelet>();
        ArrayList<Shapelet> oneLengthBestShapelet = new ArrayList<Shapelet>();


        for (int i = 0; i < data.size(); i++) {
            ArrayList<TFSAWord> TFSAListofT;
            KeyPoints.clear();
            int len = data.get(0).size();
            ArrayList<Double> T = data.get(i);
            AS(T, Math.tan(30.0), len / min_len);
            TFSAListofT = Symbolization(T);
            DataofTFSA.add(TFSAListofT);
        }

        for (subseq_len = min_len; subseq_len <= max_len; subseq_len += step) {
            /// Shapelet cannot be too short, e.g. len=1.
            if (subseq_len < SH_MIN_LEN) continue;

            int sax_len = sax_max_len;
            int w = (int) Math.ceil((1.0 * subseq_len / sax_len));//ceil():上取整
            sax_len = (int) Math.ceil(1.0 * subseq_len / w);

            CreateTFSAList3(sax_len);

            //printf("R=%d,percent_mask=%e,sax_len=%d\n\n",R,percent_mask,sax_len);
            //RandomProjection(R, percent_mask, sax_len);
            //printf("b");
            TFSARandomProjection(R, percent_mask, sax_len);

            //ScoreAllSAX(R);
            //printf("c");
            ScoreAllTFSA(R);

            //result  = FindBestSAX(top_k);
            //printf("d");
            oneLengthBestShapelet = FindBestTFSA(top_k);
            //printf("e\n");


            for (int num = 0; num < oneLengthBestShapelet.size(); num++) {
                allLengthBestShapelet.add(oneLengthBestShapelet.get(num));
            }


            UTFSAMap.clear();
            TFSASet.clear();
            TFSANewSet.clear();
            TFSAScoreList.clear();
        }


        return sh;
    }

    private double[][] readData(ArrayList<ArrayList<Double>> data) {

        double[][] d = new double[data.size()][data.get(0).size()];
        for(int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data.get(i).size(); j++) {
                d[i][j] = data.get(i).get(j);
            }
        }
        return d;
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

    public void CreateTFSAList3(int tfsa_len) {
        int id, j;
        int k = tfsa_len;
        int l = 10;
        double angle_tol = Math.tan(30.0);
        UTFSAElm ptr = new UTFSAElm();
        ArrayList<TFSAWord> TFSAList;

        for (id = 0; id < DataofTFSA.size(); id++) {
            for (j = 0; j < Math.max(DataofTFSA.get(id).size() - tfsa_len, 0); j++) {
                TFSAList = SubTFSAListofT(DataofTFSA.get(id), j, tfsa_len);
                int len_of_subseq = 0;
                //printf("DataofTFSA[id].size()=%d tfsa_len=%d j=%d",DataofTFSA[id].size(),tfsa_len,j);
                len_of_subseq = DataofTFSA.get(id).get(j + tfsa_len - 1).getPosOfBeginInData()
                        - DataofTFSA.get(id).get(j).getPosOfBeginInData()
                        + DataofTFSA.get(id).get(j + tfsa_len - 1).getLen();
                int x = isinclude(TFSAList, TFSASet);
                if (x == -1) {
                    TFSASet.add(TFSAList);
                    ptr = UTFSAMap.get(TFSASet.size() - 1);

                    HashSet tempObjSet = ptr.objSet;
                    tempObjSet.add(id);
                    ptr.setObjSet(tempObjSet);

                    ArrayList tempTfsaID = ptr.getTfsaID();
                    tempTfsaID.add(new Pair(id, DataofTFSA.get(id).get(j).getPosOfBeginInData()));
                    ptr.setTfsaID(tempTfsaID);

                    ptr.setLenOfSubseqInData(len_of_subseq);

                    UTFSAMap.put(TFSASet.size() - 1, ptr);
                } else {
                    ptr = UTFSAMap.get(x);

                    HashSet tempObjSet = ptr.objSet;
                    tempObjSet.add(id);
                    ptr.setObjSet(tempObjSet);

                    ArrayList tempTfsaID = ptr.getTfsaID();
                    tempTfsaID.add(new Pair(id, DataofTFSA.get(id).get(j).getPosOfBeginInData()));
                    ptr.setTfsaID(tempTfsaID);

                    ptr.setLenOfSubseqInData(len_of_subseq);

                    UTFSAMap.put(x, ptr);
                }
                KeyPoints.clear();
            }
        }
    }

    public ArrayList<TFSAWord> SubTFSAListofT(ArrayList<TFSAWord> T, int pos, int len) {
        ArrayList<TFSAWord> subTfsa = new ArrayList<TFSAWord>();
        int i;
        for (i = pos; i < Math.min(pos + len, T.size()); i++) {
            subTfsa.add(T.get(i));
        }
        return subTfsa;
    }

    public int isinclude(ArrayList<TFSAWord> TFSAList, ArrayList<ArrayList<TFSAWord>> TFSASet) {
        int pos = -1;
        ArrayList<TFSAWord> TFSAList_temp;
        if (TFSASet.size() == 0) {
            return pos;
        }
        for (int i = 0; i < TFSASet.size(); i++) {
            TFSAList_temp = TFSASet.get(i);
            if (Isequal(TFSAList_temp, TFSAList)) {
                return i;
            }
        }
        return pos;
    }

    public boolean Isequal(ArrayList<TFSAWord> TFSAList1, ArrayList<TFSAWord> TFSAList2) {
        if (TFSAList1.size() == TFSAList2.size()) {
            for (int i = 0; i < TFSAList1.size(); i++) {
                if ((TFSAList1.get(i).getDegree() == TFSAList2.get(i).getDegree())
                        && (TFSAList1.get(i).getLastPoint() == TFSAList2.get(i).getLastPoint()))
                //if((TFSAList1[i].degree==TFSAList2[i].degree)&&(TFSAList1[i].last_point==TFSAList2[i].last_point)&&(TFSAList1[i].trend==TFSAList2[i].trend))
                {

                } else {
                    return false;
                }
            }
        } else {
            return false;
        }
        return true;
    }

    public void TFSARandomProjection(int R, double percent_mask, int tfsa_len) {
        HashMap<Integer, HashSet<Integer>> TFSA_Hash_Mark = new HashMap<Integer, HashSet<Integer>>();

        ArrayList<TFSAWord> tfsa_list, new_list;
        int mask_word;
        HashSet<Integer> obj_set = new HashSet<Integer>();
        HashSet<Integer> ptr = new HashSet<Integer>();
        int num_mask = (int) Math.ceil(percent_mask * tfsa_len);

        for (int r = 0; r < R; r++) {
            mask_word = CreateMaskWord(num_mask, tfsa_len);
            /// random projection and mark non-duplicate object
            for (Map.Entry<Integer, UTFSAElm> entry : UTFSAMap.entrySet()) {
                tfsa_list = TFSASet.get(entry.getKey());
                obj_set = entry.getValue().getObjSet();
                new_list = OR(tfsa_list, mask_word);
                TFSANewSet.set(entry.getKey(), new_list);
                ptr = TFSA_Hash_Mark.get(entry.getKey());
                ptr.addAll(obj_set);
                TFSA_Hash_Mark.put(entry.getKey(), ptr);
            }
            /// hash again for keep the count
            for (Map.Entry<Integer, UTFSAElm> entry : UTFSAMap.entrySet()) {
                obj_set = TFSA_Hash_Mark.get(entry.getKey());
                Iterator iterator = obj_set.iterator();

                UTFSAElm temp = entry.getValue();
                HashMap<Integer, Integer> objCount = temp.getObjCount();
                while (iterator.hasNext()) {
                    int x = (Integer) iterator.next();
                    objCount.put(x, objCount.get(x) + 1);
                    temp.setObjCount(objCount);
                    UTFSAMap.put(entry.getKey(), temp);
                }

            }
            TFSA_Hash_Mark.clear();
        }
    }

    public int CreateMaskWord(int num_mask, int word_len) {
        int a, b;
        a = 0;
        Random r = new Random();
        for (int i = 0; i < num_mask; i++) {
            b = 1 << r.nextInt(word_len);
            a = a | b;
        }
        return a;
    }

    public ArrayList<TFSAWord> OR(ArrayList<TFSAWord> tfsa_list, int mask_word) {
        ArrayList<TFSAWord> new_list = new ArrayList<TFSAWord>();
        for (int i = 0; i < tfsa_list.size(); i++) {
            if ((mask_word & (1 << (tfsa_list.size() - i - 1))) != 0) {
                new_list.add(tfsa_list.get(i));
            } else {
                TFSAWord tfsa_word_null = new TFSAWord();
                tfsa_word_null.setLastPoint('N');
                new_list.add(tfsa_word_null);
            }
        }
        return new_list;
    }

    public void ScoreAllTFSA(int R) {
        ArrayList<TFSAWord> tfsa_list;
        double score;
        UTFSAElm utfsa;

        for (Map.Entry<Integer, UTFSAElm> entry : UTFSAMap.entrySet()) {
            tfsa_list = TFSASet.get(entry.getKey());
            utfsa = entry.getValue();
            score = CalScoreTFSA(utfsa, R);
            TFSAScoreList.add(new Pair(entry.getKey(), score));
        }
    }

    public double CalScoreTFSA(UTFSAElm utfsa, int R) {
        double score = -1;
        int cid, count;
        int[] c_in = new int[num_class];
        int[] c_out = new int[num_class];
        for (int i = 0; i < num_class; i++) {
            c_in[i] = 0;
            c_out[i] = 0;
        }

/// Note that if no c_in, then no c_out of that object
        HashMap<Integer, Integer> objCount = utfsa.getObjCount();
        for (Map.Entry<Integer, Integer> entry : objCount.entrySet()) {
            cid = Label.get(entry.getKey());
            count = entry.getValue();
            c_in[cid] += (count);
            c_out[cid] += (R - count);
        }
        score = CalScoreFromObjCount(c_in, c_out);
        return score;
    }

    public double CalScoreFromObjCount(int[] c_in, int[] c_out) {
        // /// 2 classes only
        //return abs((c_in[0]+c_out[1]) - (c_out[0]+c_in[1]));

        /// multi-class
        double diff, sum = 0, max_val = Double.MIN_VALUE, min_val = Double.MAX_VALUE;
        for (int i = 0; i < num_class; i++) {
            diff = (c_in[i] - c_out[i]);
            if (diff > max_val) max_val = diff;
            if (diff < min_val) min_val = diff;
            sum += Math.abs(diff);
        }
        return (sum - Math.abs(max_val) - Math.abs(min_val)) + Math.abs(max_val - min_val);
    }

    public ArrayList<Shapelet> FindBestTFSA(int top_k) {
        ArrayList<Shapelet> shapelets = new ArrayList<Shapelet>();

        ArrayList<Pair<Integer, Double>> Dist = new ArrayList<Pair<Integer, Double>>(num_obj);//建立长度为num_obj的数组
        ArrayList<TFSAWord> tfsa_list = new ArrayList<TFSAWord>();
        double gain, dist_th, gap;
        int[] bsf_c_in = new int[num_class];
        int[] bsf_c_out = new int[num_class];
        int q_obj, q_pos;
        UTFSAElm utfsa;
        int label, kk, total_c_in, num_diff;

        Shapelet sh = new Shapelet();
        Shapelet bsf_sh = new Shapelet();
        int[] c_in = new int[num_class];
        int[] c_out = new int[num_class];

        if (top_k > 0)
            Collections.sort(TFSAScoreList, new Comparator<Pair<Integer, Double>>() {
                @Override
                public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
                    if (o1.getValue() < o2.getValue())
                        return 1;
                    else if (o1.getValue() > o2.getValue())
                        return -1;
                    return 0;
                }
            });
        top_k = Math.abs(top_k);

        for (int k = 0; k < Math.min(top_k, (int) TFSAScoreList.size()); k++) {
            int TFSASet_id = TFSAScoreList.get(k).getKey();
            tfsa_list = TFSASet.get(TFSASet_id);
            utfsa = UTFSAMap.get(TFSASet_id);
            for (kk = 0; kk < Math.min((int) utfsa.getTfsaID().size(), 1); kk++) {
                q_obj = utfsa.getTfsaID().get(kk).getKey();
                q_pos = utfsa.getTfsaID().get(kk).getValue();
                int len_of_subseq = utfsa.getLenOfSubseqInData();
                double[] query = new double[len_of_subseq];

                for (int i = 0; i < num_class; i++) {
                    c_in[i] = 0;
                    c_out[i] = Class_Freq[i];
                }
                //printf("x");
                //printf("len_of_subseq=%d q_pos=%d Data[].size=%d ",len_of_subseq,q_pos,Data[q_obj].size());
                for (int i = 0; i < len_of_subseq; i++) {  // printf("i=%d ",i);
                    query[i] = Data[q_obj][q_pos + i];
                }
                //printf("y");
                double dist;
                int m = query.length;
                double[] Q = new double[m];
                int[] order = new int[m];
                for (int obj = 0; obj < num_obj; obj++) {
                    dist = NearestNeighborSearch(query, Data[obj], obj, Q, order);
                    Dist.set(obj, new Pair<Integer, Double>(obj, dist));//<时间序列序号，这条子序列到该条时间序列的距离>
                }

                Collections.sort(Dist, new Comparator<Pair<Integer, Double>>() {
                    @Override
                    public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
                        if (o1.getValue() > o2.getValue())
                            return 1;
                        else if (o1.getValue() < o2.getValue())
                            return -1;
                        return 0;
                    }
                });

                total_c_in = 0;
                for (int i = 0; i < (int) Dist.size() - 1; i++) {
                    dist_th = (Dist.get(i).getValue() + Dist.get(i + 1).getValue()) / 2;
                    //////gap = Dist[i+1].second - dist_th;
                    gap = (Dist.get(i + 1).getValue() - dist_th) / Math.sqrt((double) subseq_len);
                    label = Label.get(Dist.get(i).getKey());
                    c_in[label]++;
                    c_out[label]--;
                    total_c_in++;
                    num_diff = Math.abs(num_obj - 2 * total_c_in);
                    ///////gain = CalInfoGain1(c_in, c_out);
                    gain = CalInfoGain2(c_in, c_out, total_c_in, num_obj - total_c_in);

                    sh.setValueFew(gain, gap, dist_th);

                    if (bsf_sh.compareTo(sh) < 0) {
                        bsf_sh.setValueAll(gain, gap, dist_th, q_obj, q_pos, subseq_len, num_diff, c_in, c_out);
                    }


                }
            }
            shapelets.add(bsf_sh);
        }

        return shapelets;
    }

    public double NearestNeighborSearch(double[] query, double[] data, int obj_id, double[] Q, int[] order) {
        double d;
        long i;
        int j;
        double ex, ex2, mean, std;
        long loc = 0;

        int m = query.length;
        int M = data.length;

        double bsf = Double.MAX_VALUE;
        i = 0;
        j = 0;
        ex = ex2 = 0;

        if (obj_id == 0) {
            for (i = 0; i < m; i++) {
                d = query[(int) i];
                ex += d;
                ex2 += d * d;
                Q[(int) i] = d;
            }

            mean = ex / m;
            std = ex2 / m;
            std = Math.sqrt(std - mean * mean);

            for (i = 0; i < m; i++)
                Q[(int) i] = (Q[(int) i] - mean) / std;

            ArrayList<Pair<Double, Integer>> Q_tmp = new ArrayList<Pair<Double, Integer>>(m);

            for (i = 0; i < m; i++) {
                Q_tmp.add(new Pair(Q[(int) i], i));
            }

            Collections.sort(Q_tmp, new Comparator<Pair<Double, Integer>>() {
                @Override
                public int compare(Pair<Double, Integer> o1, Pair<Double, Integer> o2) {
                    if (Math.abs(o1.getKey()) - Math.abs(o2.getKey()) > 0)
                        return -1;
                    else if (Math.abs(o1.getKey()) - Math.abs(o2.getKey()) < 0)
                        return 1;
                    return 0;
                }
            });

            for (i = 0; i < m; i++) {
                Q[(int) i] = Q_tmp.get((int) i).getKey();
                order[(int) i] = Q_tmp.get((int) i).getValue();
            }
        }

        i = 0;
        j = 0;
        ex = ex2 = 0;


        double[] T = new double[2 * m];

        double dist = 0;
        while (i < M) {
            d = data[(int) i];
            ex += d;
            ex2 += d * d;
            T[(int) (i % m)] = d;
            T[(int) ((i % m) + m)] = d;


            if (i >= m - 1) {
                mean = ex / m;
                std = ex2 / m;
                std = Math.sqrt(std - mean * mean);

                j = (int) (i + 1) % m;
                dist = distance(Q, order, T, j, m, mean, std, bsf);
                if (dist < bsf) {
                    bsf = dist;
                    loc = i - m + 1;
                }
                ex -= T[j];
                ex2 -= T[j] * T[j];
            }
            i++;
        }
        return bsf;
    }

    public double CalInfoGain2(int[] c_in, int[] c_out, int total_c_in, int total_c_out) {
        return class_entropy - ((double) (total_c_in) / num_obj * EntropyArray(c_in, total_c_in) + (double) (total_c_out) / num_obj * EntropyArray(c_out, total_c_out));
    }

    /// new function still in doubt (as in Mueen's paper)
    public double EntropyArray(int[] A, int total) {
        double en = 0;
        double a;
        for (int i = 0; i < num_class; i++) {
            a = (double) A[i] / total;
            if (a > 0)
                en -= a * Math.log(a);
        }
        return en;
    }

    public double distance(double[] Q, int[] order, double[] T, int j, int m, double mean, double std, double best_so_far) {
        best_so_far = Double.MAX_VALUE;
        int i;
        double sum = 0;
        double bsf2 = best_so_far * best_so_far;
        for (i = 0; i < m && sum < bsf2; i++) {
            double x = (T[(order[i] + j)] - mean) / std;
            sum += (x - Q[i]) * (x - Q[i]);
        }
        return Math.sqrt(sum);
    }

}
