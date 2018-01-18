package MainPackage;

import java.util.ArrayList;
import java.util.Comparator;

public class Shapelet implements Comparable<Shapelet>{
    double  gain;
    double  gap;
    double  dist_th;
    int     obj;
    int     pos;
    int     len;
    int     num_diff;
    ArrayList<Integer> c_in;
    ArrayList<Integer> c_out;
    ArrayList<Double> ts;

    public Shapelet() {
        this.gain = Double.MIN_VALUE;
        this.gap = Double.MIN_VALUE;
        this.dist_th = Double.MAX_VALUE;
        this.obj = -1;
        this.pos = -1;
        this.len = -1;
        this.num_diff = -1;
    }

    public Shapelet(double gain, double gap, double dist_th, int obj, int pos, int len, int num_diff) {
        gain = -1;
        gap = 0;
        dist_th = 0;
        obj = 0;
        pos = 0;
        len = 0;
        num_diff = 0;
        this.gain = gain;
        this.gap = gap;
        this.dist_th = dist_th;
        this.obj = obj;
        this.pos = pos;
        this.len = len;
        this.num_diff = num_diff;
    }


    @Override
    public int compareTo(Shapelet s) {

        if((this.gain > s.gain) ||
                (this.gain == s.gain && this.num_diff < this.num_diff) ||
                (this.gain == s.gain && this.num_diff == s.num_diff && this.gap > s.gap))
            return 1;
        else if ((this.gain < s.gain) ||
                (this.gain == s.gain && this.num_diff > s.num_diff) ||
                (this.gain == s.gain && this.num_diff == s.num_diff && this.gap < s.gap)){
            return -1;
        }

        if(this.obj == s.obj && this.pos == s.pos && this.len == s.len)
            return 0;

        return 0;
    }


    public void setValueFew(double gain, double gap, double dist_th){
        this.gain = gain;
        this.gap = gap;
        this.dist_th = dist_th;
    }

    public void setValueAll(double gain, double gap, double dist_th, int obj, int pos, int len, int num_diff, ArrayList<Integer> c_in, ArrayList<Integer> c_out){
        this.gain = gain;
        this.gap = gap;
        this.dist_th = dist_th;
        this.obj = obj;
        this.pos = pos;
        this.len = len;
        this.num_diff = num_diff;
        this.c_in = c_in;
        this.c_out = c_out;
    }

    public void setValueTree( int obj, int pos, int len ,double dist_th){
        this.dist_th = dist_th;
        this.obj = obj;
        this.pos = pos;
        this.len = len;
    }

    public void setTS(ArrayList<Double> ts){
        this.ts = ts;
    }

}
