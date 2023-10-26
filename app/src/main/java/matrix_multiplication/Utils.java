package matrix_multiplication;

public class Utils {

  public void pretty_print(int[][] matrix) {
    for (int i = 0; i < matrix.length; i++) {
      System.out.print("[");
      for (int j = 0; j < matrix[0].length; j++) {
        System.out.print(matrix[i][j]);
        if (j != matrix[0].length - 1) {
          System.out.print(", ");
        }
      }
      System.out.println("]");
    }
  }

  public void pretty_print(double[][] matrix) {
    for (int i = 0; i < matrix.length; i++) {
      System.out.print("[");
      for (int j = 0; j < matrix[0].length; j++) {
        System.out.print(matrix[i][j]);
        if (j != matrix[0].length - 1) {
          System.out.print(", ");
        }
      }
      System.out.println("]");
    }
  }

  public int[] read_file(String path) {
    int[] arr = null;

    try {
      File file = new File(path);
      Scanner scanner = new Scanner(file);
      scanner.useDelimiter(",\\s*");

      int cantidadElementos = 0;
      while (scanner.hasNextInt()) {
        cantidadElementos++;
        scanner.nextInt();
      }

      scanner.close(); // NOTE: maybe not necessary
      scanner = new Scanner(file);
      scanner.useDelimiter(",\\s*");

      arr = new int[cantidadElementos];

      for (int i = 0; i < cantidadElementos; i++) {
        arr[i] = scanner.nextInt();
      }

      scanner.close();
    } catch (FileNotFoundException e) {
      System.out.println("File not found: " + path);
      e.printStackTrace();
    }

    return arr;
  }
}
