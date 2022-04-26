public class Test1 {
    static int r = 128;
    static int Nx = r;
    static int Ny = r;
    static int type = 2;

    static double aX = 1.0;
    static double bX = 2.0;
    static double cY = 1.0;
    static double dY = 2.0;

    public static double a[][][] = new double[Nx + 1][Ny + 1][Ny + 1];
    public static double b[][][] = new double[Nx + 1][Ny + 1][Ny + 1];
    public static double c[][][] = new double[Nx + 1][Ny + 1][Ny + 1];
    public static double v[][] = new double[Nx + 1][Ny + 1];
    public static double f[][] = new double[Nx + 1][Ny + 1];

    public static double x[] = new double[Nx + 1];
    public static double y[] = new double[Ny + 1];
    public static double helpX[] = new double[Nx];
    public static double helpY[] = new double[Ny];
    public static double hx;
    public static double hy;


    static double k1(double x, double y, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1 -> x + y + 1;
            case 2 -> x * x + y + 1;
            default -> 0;
        };
    }

    static double k2(double x, double y, int type) {

        return switch (type) {
            case 0 -> 1.0;
            case 1 -> x + y + 1;
            case 2 -> x + y * y + 1;
            default -> 0;
        };
    }

    static double g1(double y, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1, 2 -> y + aX;
            default -> 0;
        };
    }

    static double g2(double y, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1 -> 4 * y + 9;
            case 2 -> 3 * bX + 4 * y + bX * bX + 1;
            default -> 0;
        };
    }

    static double g3(double x, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1 -> 2 * x + 1;
            case 2 -> 2 * x + 3 * cY - cY * cY - 1;
            default -> 0;
        };
    }

    static double g4(double x, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1, 2 -> dY + x;
            default -> 0;
        };
    }

    static double F(double x, double y, int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1 -> -2.0;
            case 2 -> -2 * x - 2 * y;
            default -> 0;
        };
    }

    static double kappa3(int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1, 2 -> 3.0;
            default -> 0;
        };
    }

    static double kappa2(int type) {
        return switch (type) {
            case 0 -> 1.0;
            case 1, 2 -> 3.0;
            default -> 0;
        };
    }

    static double U(double x, double y, int type) {

        return switch (type) {
            case 0 -> 1.0;
            case 1, 2 -> x + y;
            default -> 0;
        };
    }

    static void printAllMatrix(double[][][] tmpM) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                for (int k = 0; k < Nx; k++) {
                    System.out.print(tmpM[j][i][k] + " ");
                }
                System.out.println();
            }
            System.out.println("******************");
        }
    }

    public static double[][] inversion(double[][] array, int N) {
        double temp;
        double[][] E = new double[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                E[i][j] = 0f;

                if (i == j) {
                    E[i][j] = 1f;
                }
            }
        for (int k = 0; k < N; k++) {
            temp = array[k][k];
            for (int j = 0; j < N; j++) {
                array[k][j] /= temp;
                E[k][j] /= temp;
            }
            for (int i = k + 1; i < N; i++) {
                temp = array[i][k];
                for (int j = 0; j < N; j++) {
                    array[i][j] -= array[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int k = N - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                temp = array[i][k];
                for (int j = 0; j < N; j++) {
                    array[i][j] -= array[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int i = 0; i < N; i++) {
            System.arraycopy(E[i], 0, array[i], 0, N);
        }
        return array;
    }

    public static double[][] multiplicationOfMatrix(double[][] firstM, double[][] secondM) {
        double[][] res = new double[firstM.length][firstM.length];
        for (int i = 0; i < firstM.length; i++) {
            for (int j = 0; j < firstM.length; j++) {
                for (int k = 0; k < firstM.length; k++) {
                    res[i][j] += firstM[i][k] * secondM[k][j];
                }
            }
        }
        return res;
    }

    public static double[] multiplicationMatrixAndVector(double[][] matrix, double[] vector) {
        double[] res = new double[matrix.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < vector.length; ++j)
                res[i] += matrix[i][j] * vector[j];
        }
        return res;
    }

    public static double[][] minusOne(double[][] matrix) {
        double[][] res = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[0].length; ++j)
                res[i][j] = -1 * matrix[i][j];
        }
        return res;
    }

    public static double[][] plusMatrix(double[][] m1, double[][] m2) {
        double[][] res = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; ++i) {
            for (int j = 0; j < m1[0].length; ++j)
                res[i][j] = m1[i][j] + m2[i][j];
        }
        return res;
    }

    public static double[] minusVector(double[] m1, double[] m2) {
        double[] res = new double[m1.length];
        for (int i = 0; i < m1.length; ++i) {
            res[i] = m1[i] - m2[i];
        }
        return res;
    }

    public static double[] plusVector(double[] m1, double[] m2) {
        double[] res = new double[m1.length];
        for (int i = 0; i < m1.length; ++i) {
            res[i] = m1[i] + m2[i];
        }
        return res;
    }

    public static void printMatrix(double[][] m1) {
        for (double[] doubles : m1) {
            for (int j = 0; j < m1[0].length; ++j) {
                System.out.print(doubles[j] + " ");
            }
            System.out.println();
        }
    }

    public static double[][] solveSystem(double[][][] massA, double[][][] massB, double[][][] massC, double[][] func, int x, int y) {
        double[][][] alpha = new double[y][x][y];
        double[][] beta = new double[x][y];
        double[][] res = new double[x][y];
        double[][] tmpBeta;
        double[][] tmpRes;

        tmpRes = inversion(massC[0], x);
        double[][] tmpRes2 = minusOne(tmpRes);
        alpha[0] = multiplicationOfMatrix(tmpRes2, massB[0]);

        tmpBeta = inversion(massC[0], x);
        beta[0] = multiplicationMatrixAndVector(tmpBeta, func[0]);
        for (int i = 1; i < y; i++) {
            double[][] ttt = multiplicationOfMatrix(massA[i - 1], alpha[i - 1]);
            double[][] tmpResCom = plusMatrix(ttt, massC[i - 1]);
            double[][] tmpResCom1 = inversion(tmpResCom, x);
            double[][] tmpRes1 = minusOne(tmpResCom1);
            alpha[i] = multiplicationOfMatrix(tmpRes1, massB[i - 1]);

            double[][] tmpResCom3 = plusMatrix(multiplicationOfMatrix(massA[i - 1], alpha[i - 1]), massC[i - 1]);
            double[][] tmpResCom4 = inversion(tmpResCom3, x);
            double[] tmpBeta1 = minusVector(func[i - 1], multiplicationMatrixAndVector(massA[i - 1], beta[i - 1]));
            beta[i] = multiplicationMatrixAndVector(tmpResCom4, tmpBeta1);
        }

        double[][] tmpResCom = plusMatrix(multiplicationOfMatrix(massA[y - 1], alpha[y - 1]), massC[y - 1]);
        double[][] tmpResCom1 = inversion(tmpResCom, x);
        double[] tmpBeta1 = minusVector(func[y - 1], multiplicationMatrixAndVector(massA[y - 1], beta[y - 1]));
        res[y - 1] = multiplicationMatrixAndVector(tmpResCom1, tmpBeta1);

        for (int i = y - 2; i >= 0; i--) {
            res[i] = plusVector(multiplicationMatrixAndVector(alpha[i + 1], res[i + 1]), beta[i + 1]);
        }

        return res;
    }

    public static void main(String[] args) {
        //setka
        double stepX = (bX - aX) / Nx;
        double stepY = (dY - cY) / Ny;
        for (int i = 0; i < Nx + 1; i++) {
            x[i] = aX + stepX * i;
        }
        helpX[0] = aX + stepX / 2;
        for (int i = 1; i < Nx; i++) {
            helpX[i] = x[i] + stepX / 2;
        }
//        helpX[Nx] = x[Nx - 1] - stepX / 2;
        for (int i = 0; i < Ny + 1; i++) {
            y[i] = cY + stepY * i;
        }
        helpY[0] = cY + stepY / 2;
        for (int i = 1; i < Ny; i++) {
            helpY[i] = y[i] + stepY / 2;
        }
//        helpY[Nx] = y[Ny - 1] - stepY / 2;
        hx = stepX;
        hy = stepY;

        //obnulit
        for (int j = 0; j < Ny + 1; j++) {
            for (int i = 0; i < Nx + 1; i++) {
                for (int k = 0; k < Nx + 1; k++) {
                    a[j][i][k] = 0;
                    b[j][i][k] = 0;
                    c[j][i][k] = 0;
                }
            }
        }

        // make matrices
        // border down
        int j = 0;
        a[j][0][0] = 0;
        b[j][0][0] = 0;
        c[j][0][0] = 1;
        f[j][0] = g1(y[j], type);
        for (int i = 1; i < Nx; i++) {
            a[j][i][i] = 0;     // coefff i j - 1
            c[j][i][i - 1] = -(hy / 2 * k1(helpX[i - 1], y[j], type)) / hx; // i-1, j
            c[j][i][i] = hy / 2 * (k1(helpX[i], y[j], type) / hx + k1(helpX[i - 1], y[j], type) / hx) +
                    hx * (k2(x[i], helpY[j], type) / hy + kappa3(type)); // i j
            c[j][i][i + 1] = -(hy / 2 * k1(helpX[i], y[j], type)) / hx; // i+1, j
            b[j][i][i] = -(hx * k2(x[i], helpY[j], type)) / hy; // i, j+1

            f[j][i] = hx * hy / 2 * F(x[i], y[j], type) + hx * g3(x[i], type);
        }
        // angle dot
        j = 0;
        int I = Nx;
        a[j][I][I] = 0;
        b[j][I][I] = -(hx / 2 * (k2(x[I], helpY[j], type) / hy));
        c[j][I][I - 1] = -(hy / 2 * k1(helpX[I - 1], y[j], type) / hx);
        c[j][I][I] = hy / 2 * (kappa2(type) + k1(helpX[I - 1], y[j], type) / hx) +
                hx / 2 * (k2(x[I], helpY[j], type) / hy + kappa3(type));
        f[j][I] = hx * hy / 4 * F(x[I], y[j], type) + hx / 2 * g3(x[I], type) + hy / 2 * g2(y[j], type);

        //plain part
        for (j = 1; j < Ny; j++) {
            a[j][0][0] = 0;
            b[j][0][0] = 0;
            c[j][0][0] = 1;
            f[j][0] = g1(y[j], type);
            for (int i = 1; i < Nx; i++) {
                a[j][i][i] = -(hx * k2(x[i], helpY[j - 1], type)) / hy;     // coefff i j - 1
                c[j][i][i - 1] = -(hy * k1(helpX[i - 1], y[j], type)) / hx; // i-1, j
                c[j][i][i] = hy * (k1(helpX[i], y[j], type) / hx + k1(helpX[i - 1], y[j], type) / hx) +
                        hx * (k2(x[i], helpY[j], type) / hy + k2(x[i], helpY[j - 1], type) / hy); // i j
                c[j][i][i + 1] = -(hy * k1(helpX[i], y[j], type)) / hx; // i+1, j
                b[j][i][i] = -(hx * k2(x[i], helpY[j], type)) / hy; // i, j+1

                f[j][i] = hx * hy * F(x[i], y[j], type);
            }

            // right border
            I = Nx;
            a[j][I][I] = -(hx / 2 * k2(x[I], helpY[j - 1], type)) / hy;     // coefff i j - 1
            c[j][I][I - 1] = -(hy * k1(helpX[I - 1], y[j], type)) / hx; // i-1, j
            c[j][I][I] = hy * (kappa2(type) + k1(helpX[I - 1], y[j], type) / hx) +
                    hx / 2 * (k2(x[I], helpY[j], type) / hy + k2(x[I], helpY[j - 1], type) / hy);
            b[j][I][I] = -(hx / 2 * k2(x[I], helpY[j], type)) / hy; // i, j+1

            f[j][I] = hx * hy / 2 * F(x[I], y[j], type) + hy * g2(y[j], type);

        }

        // upper border
        j = Ny;
        for (int i = 0; i < Nx + 1; i++) {
            a[j][i][i] = 0;
            b[j][i][i] = 0;
            c[j][i][i] = 1;
            f[j][i] = g4(x[i], type);
        }
        printAllMatrix(c);
        double[][] answer = solveSystem(a, b, c, f, Nx + 1, Ny + 1);
        for (int i = 0; i < Nx + 1; i++) {
            for (int k = 0; k < Ny + 1; k++) {
                System.out.print(answer[i][k] + " ");
            }
            System.out.println();
        }
        System.out.println("*******************************");
        double max = Math.abs(U(x[0], y[0], type) - answer[0][0]);
        for (int i = 0; i < Nx + 1; i++) {
            for (int k = 0; k < Ny + 1; k++) {
                double tmp = Math.abs(U(x[i], y[k], type) - answer[i][k]);
                System.out.print(tmp + " ");
                if (tmp > max) {
                    max = tmp;
                }
            }
            System.out.println();
        }
        System.out.println("*******************************");
        System.out.println("MAX: " + max);
    }

}
