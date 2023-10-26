package matrix_multiplication;

public class Naiv {
  public int[][] NaivKahan(int[][] matrix1, int[][] matrix2) {
    int rows1 = matrix1.length;
    int cols1 = matrix1[0].length;
    int cols2 = matrix2[0].length;

    if (cols1 != matrix2.length) {
      throw new IllegalArgumentException("Las matrices no son multiplicables");
    }

    int[][] result = new int[rows1][cols2];

    for (int i = 0; i < rows1; i++) {
      for (int j = 0; j < cols2; j++) {
        for (int k = 0; k < cols1; k++) {
          result[i][j] += matrix1[i][k] * matrix2[k][j];
        }
      }
    }

    return result;
  }

  public double[][] NaiveLoopUnrollingTwo(double[][] A, double[][] B) {
    int n = A.length;
    int m = A[0].length;
    int p = B[0].length;

    if (m != B.length) {
      throw new IllegalArgumentException(
          "No se pueden multiplicar estas matrices.");
    }

    double[][] result = new double[n][p];

    int blockSize = 64; // Tamaño del bloque, ajusta según la caché de tu CPU

    for (int i = 0; i < n; i += blockSize) {
      for (int j = 0; j < p; j += blockSize) {
        for (int k = 0; k < m; k += blockSize) {
          for (int i1 = i; i1 < Math.min(i + blockSize, n); i1++) {
            for (int j1 = j; j1 < Math.min(j + blockSize, p); j1++) {
              for (int k1 = k; k1 < Math.min(k + blockSize, m); k1++) {
                result[i1][j1] += A[i1][k1] * B[k1][j1];
              }
            }
          }
        }
      }
    }

    return result;
  }
}
