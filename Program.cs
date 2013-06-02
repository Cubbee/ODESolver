using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ODESolver
{
    class Program
    {
        static double b(double lambda, double t)
        {
            //return lambda*lambda*lambda + t;
            return 0.25;
        }

        static double dBLambda(double lambda, double t)
        {
            //return 3*lambda*lambda;
            return 0;
        }

        static double d2BLambda(double lambda, double t)
        {
            //return 6*lambda;
            return 0;
        }

        static double a(double lambda, double t)
        {
            //return lambda*lambda + t;
            return 2 * lambda;
        }

        static double dALambda(double lambda, double t)
        {
            //return 2*lambda;
            return 2;
        }

        static double alpha(double labmda, double t, double h, double l)
        {
            return l / (2.0 * h * h) * b(labmda, t) + l / (2.0 * h) * (dBLambda(labmda, t) - a(labmda, t));
        }

        static double gamma(double labmda, double t, double h, double l)
        {
            return l / (2.0 * h * h) * b(labmda, t) - l / (2.0 * h) * (dBLambda(labmda, t) - a(labmda, t));
        }

        static double delta(double labmda, double t, double h, double l)
        {
            return 1 + l / (h * h) * b(labmda, t) - l * ((1 / 2.0) * d2BLambda(labmda, t) - dALambda(labmda, t));
        }

        static double beta(double labmda, double t, double h, double l)
        {
            return 1 - l / (h * h) * b(labmda, t) + l * ((1 / 2.0) * d2BLambda(labmda, t) - dALambda(labmda, t));
        }

        static bool checkTridiagonalMatrixAlgorithmStability(double A, double B, double C, double E)
        {
            return (A > 0) && (C > 0) && (B >= A + C) && (E >= 0 && E < 1);
        }

        static bool checkExplicitAlgorithmStability(double b, double l, double h)
        {
            return (b > 0) && ((l / (h * h)) < (1 / b));
        }

        static void Main(string[] args)
        {
            double c = 0;
            double d = 2.2;

            int N = 220;
            int M = 10000;
            List<List<double>> grid = new List<List<double>>();

            for (int i = 0; i < M; i++)
            {
                List<double> temp = new List<double>(N);

                for (int j = 0; j < N; j++)
                {
                    temp.Add(0);
                }

                grid.Add(temp);
            }
            
            double h = (d - c) / N;
            double l = 0.0001;
            
            ////initial distribution P(x) normal
            //double mx = (c + d) / 2;
            //double sigma = (mx - c) / 4;

            //for (int i = 1; i < grid[0].Count - 1; i++)
            //{
            //    grid[0][i] = 1 / (sigma * Math.Sqrt(2 * Math.PI)) *
            //                 Math.Exp(-Math.Pow(c + i * h - mx, 2) / (2 * Math.Pow(sigma, 2)));
            //}

            //initial delta distibution
            double initialX = 1;
            int n = 140;
            for (int i = 1; i < grid[0].Count - 1; i++)
            {
                grid[0][i] = n / Math.Sqrt(Math.PI) *
                             Math.Exp(-(n * (c + i * h - initialX)) * (n * (c + i * h - initialX)));
            }

            //normalization
            //double sum = 0;
            //for (int i = 0; i < grid[0].Count; i++)
            //{
            //    sum += grid[0][i];
            //}

            //for (int i = 0; i < grid[0].Count; i++)
            //{
            //    grid[0][i] /= sum;
            //}

            ExplicitMethod(N, M, c, h, l, grid);
            
            //ImplicitMethod(N, M, c, h, l, grid);

            using (var writer = new StreamWriter("output.csv"))
            {
                foreach (var line in grid)
                {
                    foreach (var point in line)
                    {
                        writer.Write(point + ";");
                    }
                    writer.WriteLine();
                }
            }

            Console.WriteLine("Finished");
            Console.ReadLine();
        }

        private static void ExplicitMethod(int N, int M, double c, double h, double l, List<List<double>> grid)
        {
            for (int j = 0; j < M - 1; j++)
            {
                //слой сетки
                double t = j * l;

                for (int i = 1; i < N - 1; i++)
                {
                    double lambda = c + i * h;

                    if (!checkExplicitAlgorithmStability(b(lambda, t), l, h))
                    {
                        Console.WriteLine("Boroda na shage: j: " + j.ToString() + " i: " + i.ToString());
                    }
                    grid[j + 1][i] = alpha(lambda, t, h, l) * grid[j][i + 1] + beta(lambda, t, h, l) * grid[j][i] +
                                     gamma(lambda, t, h, l) * grid[j][i - 1];
                }
            }
        }

        private static void ImplicitMethod(int N, int M, double c, double h, double l, List<List<double>> grid)
        {
            List<double> F = new List<double>(N);
            List<double> E = new List<double>(N);

            for (int i = 0; i < N; i++)
            {
                E.Add(0);
                F.Add(0);
            }

            //boundary conditions
            double chi1 = 0;
            double nu1 = 0;
            double chi2 = 0;
            double nu2 = 0;

            for (int j = 1; j < M; j++)
            {
                //слой сетки
                double t = j * l;

                //boundary conditions
                E[0] = chi1;
                F[0] = nu1;

                for (int i = 0; i < N - 1; i++)
                {
                    double lambda = c + i * h;

                    double C = alpha(lambda, t, h, l);
                    double B = delta(lambda, t, h, l);
                    double A = gamma(lambda, t, h, l);
                    double D = grid[j - 1][i];

                    E[i + 1] = C / (B - A * E[i]);
                    F[i + 1] = (A * F[i] + D) / (B - A * E[i]);

                    if (!checkTridiagonalMatrixAlgorithmStability(A, B, C, E[i + 1]))
                    {
                        Console.WriteLine("Boroda na shage: j: " + j.ToString() + " i: " + i.ToString());
                    }
                }

                grid[j][N - 1] = (nu2 + chi2 * F[N - 1]) / (1 - chi2 * E[N - 1]); // boundary conditions
                for (int i = N - 2; i >= 1; i--)
                {
                    grid[j][i] = E[i + 1] * grid[j][i + 1] + F[i + 1];
                }

                //нормировка в слое
                //double summ = 0;
                //for (int i = 0; i < N; i++)
                //{
                //    summ += grid[j][i];
                //}
                //for (int i = 0; i < N; i++)
                //{
                //    grid[j][i] /= summ;
                //}

                for (int i = 0; i < N; i++)
                {
                    E[i] = 0;
                    F[i] = 0;
                }
            }
        }
    }
}
