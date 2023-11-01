package matrix_multiplication;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Strassen {

  public void StrassenNaiv(double[][] A, double[][] B, double[][] Result, int N,
      int P, int M) {
    int MaxSize, k, m, NewSize, i, j;
    MaxSize = Math.max(N, P);
    MaxSize = Math.max(MaxSize, M);
    if (MaxSize < 16) {
      MaxSize = 16; // de lo contrario no es posible calcular k
    }
    k = (int) (Math.log(MaxSize) / Math.log(2)) - 4;
    m = (int) (MaxSize * Math.pow(2, -k)) + 1;
    NewSize = m * (int) Math.pow(2, k);

    // Agregar filas y columnas de ceros para utilizar el algoritmo de Strassen
    double[][] NewA = new double[NewSize][NewSize];
    double[][] NewB = new double[NewSize][NewSize];
    double[][] AuxResult = new double[NewSize][NewSize];

    for (i = 0; i < NewSize; i++) {
      for (j = 0; j < NewSize; j++) {
        NewA[i][j] = 0.0;
        NewB[i][j] = 0.0;
      }
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < P; j++) {
        NewA[i][j] = A[i][j];
      }
    }

    for (i = 0; i < P; i++) {
      for (j = 0; j < M; j++) {
        NewB[i][j] = B[i][j];
      }
    }

    StrassenNaivStep(NewA, NewB, AuxResult, NewSize, m);

    // Extraer el resultado
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        Result[i][j] = AuxResult[i][j];
      }
    }
  }

  public void StrassenNaivStep(double[][] A, double[][] B, double[][] Result,
      int N, int m) {
    int i, j, NewSize;
    if ((N % 2 == 0) && (N > m)) { // Uso recursivo de StrassenNaivStep
      NewSize = N / 2;

      // Descomponer A y B
      double[][] A11 = new double[NewSize][NewSize];
      double[][] A12 = new double[NewSize][NewSize];
      double[][] A21 = new double[NewSize][NewSize];
      double[][] A22 = new double[NewSize][NewSize];
      double[][] B11 = new double[NewSize][NewSize];
      double[][] B12 = new double[NewSize][NewSize];
      double[][] B21 = new double[NewSize][NewSize];
      double[][] B22 = new double[NewSize][NewSize];
      double[][] ResultPart11 = new double[NewSize][NewSize];
      double[][] ResultPart12 = new double[NewSize][NewSize];
      double[][] ResultPart21 = new double[NewSize][NewSize];
      double[][] ResultPart22 = new double[NewSize][NewSize];
      double[][] Helper1 = new double[NewSize][NewSize];
      double[][] Helper2 = new double[NewSize][NewSize];
      double[][] Aux1 = new double[NewSize][NewSize];
      double[][] Aux2 = new double[NewSize][NewSize];
      double[][] Aux3 = new double[NewSize][NewSize];
      double[][] Aux4 = new double[NewSize][NewSize];
      double[][] Aux5 = new double[NewSize][NewSize];
      double[][] Aux6 = new double[NewSize][NewSize];
      double[][] Aux7 = new double[NewSize][NewSize];

      // Llenar las nuevas matrices
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          A11[i][j] = A[i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          A12[i][j] = A[i][NewSize + j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          A21[i][j] = A[NewSize + i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          A22[i][j] = A[NewSize + i][NewSize + j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          B11[i][j] = B[i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          B12[i][j] = B[i][NewSize + j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          B21[i][j] = B[NewSize + i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          B22[i][j] = B[NewSize + i][NewSize + j];
        }
      }

      // Calcular las siete variables auxiliares
      plus(A11, A22, Helper1, NewSize, NewSize);
      plus(B11, B22, Helper2, NewSize, NewSize);
      StrassenNaivStep(Helper1, Helper2, Aux1, NewSize, m);
      plus(A21, A22, Helper1, NewSize, NewSize);
      StrassenNaivStep(Helper1, B11, Aux2, NewSize, m);
      minus(B12, B22, Helper1, NewSize, NewSize);
      StrassenNaivStep(A11, Helper1, Aux3, NewSize, m);
      minus(B21, B11, Helper1, NewSize, NewSize);
      StrassenNaivStep(A22, Helper1, Aux4, NewSize, m);
      plus(A11, A12, Helper1, NewSize, NewSize);
      StrassenNaivStep(Helper1, B22, Aux5, NewSize, m);
      minus(A21, A11, Helper1, NewSize, NewSize);
      plus(B11, B12, Helper2, NewSize, NewSize);
      StrassenNaivStep(Helper1, Helper2, Aux6, NewSize, m);
      minus(A12, A22, Helper1, NewSize, NewSize);
      plus(B21, B22, Helper2, NewSize, NewSize);
      StrassenNaivStep(Helper1, Helper2, Aux7, NewSize, m);

      // Calcular las cuatro partes del resultado
      plus(Aux1, Aux4, ResultPart11, NewSize, NewSize);
      minus(ResultPart11, Aux5, ResultPart11, NewSize, NewSize);
      plus(ResultPart11, Aux7, ResultPart11, NewSize, NewSize);
      plus(Aux3, Aux5, ResultPart12, NewSize, NewSize);
      plus(Aux2, Aux4, ResultPart21, NewSize, NewSize);
      plus(Aux1, Aux3, ResultPart22, NewSize, NewSize);
      minus(ResultPart22, Aux2, ResultPart22, NewSize, NewSize);
      plus(ResultPart22, Aux6, ResultPart22, NewSize, NewSize);

      // Almacenar los resultados en la matriz "result"
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          Result[i][j] = ResultPart11[i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          Result[i][NewSize + j] = ResultPart12[i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          Result[NewSize + i][j] = ResultPart21[i][j];
        }
      }
      for (i = 0; i < NewSize; i++) {
        for (j = 0; j < NewSize; j++) {
          Result[NewSize + i][NewSize + j] = ResultPart22[i][j];
        }
      }
    } else {
      NaiveStandard(A, B, Result, N, N, N);
    }
  }

  public static void plus(double[][] A, double[][] B, double[][] Result, int N,
      int M) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        Result[i][j] = A[i][j] + B[i][j];
      }
    }
  }

  public void minus(double[][] A, double[][] B, double[][] Result, int N,
      int M) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        Result[i][j] = A[i][j] - B[i][j];
      }
    }
  }

  public void NaiveStandard(double[][] A, double[][] B, double[][] Result,
      int N, int P, int M) {
    double aux;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        aux = 0.0;
        for (int k = 0; k < P; k++) {
          aux += A[i][k] * B[k][j];
        }
        Result[i][j] = aux;
      }
    }
  }

  public static double[][] loadMatrix(String filePath, int rows, int cols) {
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
}
