package matrix_multiplication;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

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

  public double[][] loadMatrix(String filePath, int rows, int cols) {
    double[][] matrix = new double[rows][cols];
    try {
      BufferedReader br = new BufferedReader(new FileReader(filePath));
      for (int i = 0; i < rows; i++) {
        String[] row = br.readLine().split(" ");
        for (int j = 0; j < cols; j++) {
          matrix[i][j] = Double.parseDouble(row[j]);
        }
      }
      br.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
    return matrix;
  }

  public void NaiveLoopUnrollingThree(double[][] A, double[][] B,
      double[][] Result, int N, int P, int M) {

    int i, j, k;
    double aux;

    if (P % 3 == 0) {
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0.0;
          for (k = 0; k < P; k += 3) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j];
          }
          Result[i][j] = aux;
        }
      }
    } else if (P % 3 == 1) {
      int PP = P - 1;
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0.0;
          for (k = 0; k < PP; k += 3) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j];
          }
          Result[i][j] = aux + A[i][PP] * B[PP][j];
        }
      }
    } else {
      int PP = P - 2;
      int PPP = P - 1;
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0.0;
          for (k = 0; k < PP; k += 3) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j];
          }
          Result[i][j] = aux + A[i][PP] * B[PP][j] + A[i][PPP] * B[PPP][j];
        }
      }
    }
  }

  public void NaiveLoopUnrollingFour(int[][] A, int[][] B, int[][] Result,
      int N, int P, int M) {
    int i, j, k;
    int aux;
    if (P % 4 == 0) {
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0;
          for (k = 0; k < P; k += 4) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j] + A[i][k + 3] * B[k + 3][j];
          }
          Result[i][j] = aux;
        }
      }
    } else if (P % 4 == 1) {
      int PP = P - 1;
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0;
          for (k = 0; k < PP; k += 4) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j] + A[i][k + 3] * B[k + 3][j];
          }
          Result[i][j] = aux + A[i][PP] * B[PP][j];
        }
      }
    } else if (P % 4 == 2) {
      int PP = P - 2;
      int PPP = P - 1;
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0;
          for (k = 0; k < PP; k += 4) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j] + A[i][k + 3] * B[k + 3][j];
          }
          Result[i][j] = aux + A[i][PP] * B[PP][j] + A[i][PPP] * B[PPP][j];
        }
      }
    } else {
      int PP = P - 3;
      int PPP = P - 2;
      int PPPP = P - 1;
      for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
          aux = 0;
          for (k = 0; k < PP; k += 4) {
            aux += A[i][k] * B[k][j] + A[i][k + 1] * B[k + 1][j] +
                A[i][k + 2] * B[k + 2][j] + A[i][k + 3] * B[k + 3][j];
          }
          Result[i][j] = aux + A[i][PP] * B[PP][j] + A[i][PPP] * B[PPP][j] +
              A[i][PPPP] * B[PPPP][j];
        }
      }
    }
  }
}
