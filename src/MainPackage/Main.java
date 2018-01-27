package MainPackage;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Main {
    public static void main(String[] args)
    {
        String path = args[0];//参数设置，训练集路径
        int dimension = Integer.parseInt(args[1]);//参数设置，维度数量
        int len = Integer.parseInt(args[2]);//参数设置，时间序列长度

        ArrayList<double[][]> datas = readTrainData(path, dimension, len);

    }

    private static ArrayList<double[][]> readTrainData(String path, int dimension, int len){

        ArrayList<double[][]> datas = new ArrayList<double[][]>();
        int number;
        for (int i = 0; i < dimension; i++) {
            double[][] data = new double[dimension][len];
            File directory = new File(path);
            if (directory.exists() && directory.isDirectory()) {
                String[] filesName = directory.list();
                number = filesName.length;
                for (int n = 0; n < number; n++) {
                    Workbook book = null;
                    try {
                        book = Workbook.getWorkbook(new File(path + "\\" + filesName[n]));
                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (BiffException e) {
                        e.printStackTrace();
                    }
                    Sheet sheet = book.getSheet(0);
                    for (int l = 0; l < len; l++) {
                        Cell cell = sheet.getCell(i, l);
                        data[i][l] = Double.parseDouble(cell.getContents());
                    }
                }
            } else {
                System.out.println("文件不存在 或 路径错误");
                return null;
            }
            datas.add(data);
        }

        return datas;
    }
}
