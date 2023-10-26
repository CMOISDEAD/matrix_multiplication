package matrix_multiplication;

public class App {
    public static void main(String[] args) {
        Utils utils = new Utils();

        // int[][] matrix1 = { { 1, 2 }, { 3, 4 } };
        // int[][] matrix2 = { { 5, 6 }, { 7, 8 } };
        // int[][] result = new Naiv().NaivKahan(matrix1, matrix2);
        // utils.pretty_print(result);

        double[][] matrix1 = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        double[][] matrix2 = { { 5.0, 6.0 }, { 7.0, 8.0 } };
        double[][] result = new Naiv().NaiveLoopUnrollingTwo(matrix1, matrix2);
        utils.pretty_print(result);
    }
}
