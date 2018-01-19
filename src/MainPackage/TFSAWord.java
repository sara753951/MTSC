package MainPackage;

public class TFSAWord {
    int degree;
    char lastPoint;
    int posOfBeginInData;
    int len;

    public TFSAWord() {
        this.degree = 0;
        this.lastPoint = 'A';
        this.posOfBeginInData = 0;
        this.len = 1;
    }

    public TFSAWord(int degree, char lastPoint, int posOfBeginInData, int len) {
        this.degree = degree;
        this.lastPoint = lastPoint;
        this.posOfBeginInData = posOfBeginInData;
        this.len = len;
    }

    public int getDegree() {
        return degree;
    }

    public char getLastPoint() {
        return lastPoint;
    }

    public int getPosOfBeginInData() {
        return posOfBeginInData;
    }

    public int getLen() {
        return len;
    }

    public void setDegree(int degree) {
        this.degree = degree;
    }

    public void setLastPoint(char lastPoint) {
        this.lastPoint = lastPoint;
    }

    public void setPosOfBeginInData(int posOfBeginInData) {
        this.posOfBeginInData = posOfBeginInData;
    }

    public void setLen(int len) {
        this.len = len;
    }
}
