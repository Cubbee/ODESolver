using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ODESolver
{
    class TrigiagonalAlgorithmStabilityIndicator
    {
        public bool AGreaterThanZero { get; set; }
        public bool CGreaterThanZero { get; set; }
        public bool BGreaterThanSumOfAAndC { get; set; }
        public bool EGreaterOrEqualThanZero { get; set; }
        public bool ELesserThanOne { get; set; }

        public TrigiagonalAlgorithmStabilityIndicator()
        {
            AGreaterThanZero = false;
            CGreaterThanZero = false;
            BGreaterThanSumOfAAndC = false;
            EGreaterOrEqualThanZero = false;
            ELesserThanOne = false;
        }

        public void Reset()
        {
            AGreaterThanZero = false;
            CGreaterThanZero = false;
            BGreaterThanSumOfAAndC = false;
            EGreaterOrEqualThanZero = false;
            ELesserThanOne = false;   
        }
    }

    class ExplicitAlgorithmStabilityIndicator
    {
        public bool BGreaterThanZero {get; set;}
        public bool LOverHSquaredLesserThan1OverB { get; set; }

        public ExplicitAlgorithmStabilityIndicator()
        {
            BGreaterThanZero = false;
            LOverHSquaredLesserThan1OverB = false;
        }

        public void Reset()
        {
            BGreaterThanZero = false;
            LOverHSquaredLesserThan1OverB = false;
        }
    }

    class Program
    {
        static double b(double lambda, double t)
        {
            return lambda*lambda*lambda + t;
            //return 1;
        }

        static double dBLambda(double lambda, double t)
        {
            return 3*lambda*lambda;
            //return 0;
        }

        static double d2BLambda(double lambda, double t)
        {
            return 6*lambda;
            //return 0;
        }

        static double a(double lambda, double t)
        {
            //return lambda*lambda + t;
            return -1.1 * lambda;
        }

        static double dALambda(double lambda, double t)
        {
            //return 2*lambda;
            return -1.1;
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

        static bool checkTridiagonalMatrixAlgorithmStability(double A, double B, double C, double E, ref TrigiagonalAlgorithmStabilityIndicator indicator)
        {
            
            if (A > 0)
            {
                indicator.AGreaterThanZero = true;
            }

            if (C > 0)
            {
                indicator.CGreaterThanZero = true;
            }

            if ((B >= A + C))
            {
                indicator.BGreaterThanSumOfAAndC = true;
            }

            if (E >= 0)
            {
                indicator.EGreaterOrEqualThanZero = true;
            }

            if (E < 1)
            {
                indicator.ELesserThanOne = true;
            }

            return indicator.AGreaterThanZero && indicator.CGreaterThanZero && indicator.BGreaterThanSumOfAAndC && (indicator.EGreaterOrEqualThanZero && indicator.ELesserThanOne);
        }

        static bool checkExplicitAlgorithmStability(double b, double l, double h, ref ExplicitAlgorithmStabilityIndicator indicator)
        {
            if (b > 0)
            {
                indicator.BGreaterThanZero = true;
            }

            if ((l / (h * h)) < (1 / b))
            {
                indicator.LOverHSquaredLesserThan1OverB = true;
            }
            return indicator.BGreaterThanZero && indicator.LOverHSquaredLesserThan1OverB;
        }

        static void Main(string[] args)
        {
            double c = 0.001;
            double d = 2.5;

            int N = 50;
            int M = 1000;
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

            //initial distribution P(x) normal
            //double mx = (c + d) / 2;
            //double sigma = (mx - c) / 4;

            double mx = 1;
            double sigma = 1;

            for (int i = 1; i < grid[0].Count - 1; i++)
            {
                grid[0][i] = 1 / (sigma * Math.Sqrt(2 * Math.PI)) *
                             Math.Exp(-Math.Pow(c + i * h - mx, 2) / (2 * Math.Pow(sigma, 2)));
            }

            //initial delta distibution
            //double initialX = 0;
            //int n = 140;
            //for (int i = 1; i < grid[0].Count - 1; i++)
            //{
            //    grid[0][i] = n / Math.Sqrt(Math.PI) *
            //                 Math.Exp(-(n * (c + i * h - initialX)) * (n * (c + i * h - initialX)));
            //}

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

            ExplicitMethod(N, M, c, d, h, l, grid);
            
            //ImplicitMethod(N, M, c, d, h, l, grid);

            using (var writer = new StreamWriter("output.csv"))
            {
                for (int i = 0; i < grid.Count; i++)
                {
                    if (i % 100 == 0)
                    {
                        var line = grid[i];
                        writer.Write(l * i + ";");
                        foreach (var point in line)
                        {
                            writer.Write(point + ";");
                        }
                        writer.WriteLine();    
                    }
                    
                }
            }

            Console.WriteLine("Finished");
            Console.ReadLine();
        }

        private static void ExplicitMethod(int N, int M, double c, double d, double h, double l, List<List<double>> grid)
        {
            ExplicitAlgorithmStabilityIndicator indicator = new ExplicitAlgorithmStabilityIndicator();
            for (int j = 0; j < M - 1; j++)
            {
                //слой сетки
                double t = j * l;

                for (int i = 1; i < N - 1; i++)
                {
                    double lambda = c + i * h;

                    if (!checkExplicitAlgorithmStability(b(lambda, t), l, h, ref indicator))
                    {
                        Console.Write("Instability on step: j: " + j.ToString() + " i: " + i.ToString() + " ");
                        if (!indicator.BGreaterThanZero)
                        {
                            Console.Write("b < 0, b = {0} ", b(lambda, t));
                        }

                        if (!indicator.LOverHSquaredLesserThan1OverB)
                        {
                            Console.Write("(l / (h * h)) >= (1 / b), l / (h * h) = {0}, 1 / b = {1}", l / (h * h), 1 / b(lambda, t));
                        }
                        Console.WriteLine();
                    }
                    indicator.Reset();
                    grid[j + 1][i] = alpha(lambda, t, h, l) * grid[j][i + 1] + beta(lambda, t, h, l) * grid[j][i] +
                                     gamma(lambda, t, h, l) * grid[j][i - 1];
                    
                }

                
                grid[j + 1][0] = grid[j + 1][1] / (1 + h / b(c, t) * (2 * a(c, t) - dBLambda(c, t)));
                grid[j + 1][N - 1] = grid[j + 1][N - 2] / (1 - h * 1 / b(d, t) * (2 * a(d, t) - dBLambda(d, t)));
            }
        }

        private static void ImplicitMethod(int N, int M, double c, double d, double h, double l, List<List<double>> grid)
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
            TrigiagonalAlgorithmStabilityIndicator indicator = new TrigiagonalAlgorithmStabilityIndicator();
            for (int j = 1; j < M; j++)
            {
                //слой сетки
                double t = j * l;

                //boundary conditions
                E[0] = 0;
                E[0] = 1 / (1 + h / b( c, t ) * ( 2 * a( c, t ) - dBLambda(c, t)));
                //F[0] = 0;

                for (int i = 0; i < N - 1; i++)
                {
                    double lambda = c + i * h;

                    double C = alpha(lambda, t, h, l);
                    double B = delta(lambda, t, h, l);
                    double A = gamma(lambda, t, h, l);
                    double D = grid[j - 1][i];

                    E[i + 1] = C / (B - A * E[i]);
                    F[i + 1] = (A * F[i] + D) / (B - A * E[i]);

                    if (!checkTridiagonalMatrixAlgorithmStability(A, B, C, E[i + 1], ref indicator))
                    {
                        Console.Write("Instability on step: j: " + j.ToString() + " i: " + i.ToString() + " ");
                        if (!indicator.AGreaterThanZero)
                        {
                            Console.Write("A < 0, A = {0} ", A);
                        }
                        if (!indicator.CGreaterThanZero)
                        {
                            Console.Write("C < 0, C = {0} ", C);
                        }
                        if (!indicator.BGreaterThanSumOfAAndC)
                        {
                            Console.Write("B < A + C, B = {0}, A + C = {1} ", B, A + C);
                        }
                        if (!indicator.EGreaterOrEqualThanZero)
                        {
                            Console.Write("E < 0, E = {0} ", E);
                        }
                        if (!indicator.ELesserThanOne)
                        {
                            Console.Write("E > 1, E = {0} ", E);
                        }
                        Console.WriteLine();
                    }

                    indicator.Reset();
                }

                //grid[j][N - 1] = (nu2 + chi2 * F[N - 1]) / (1 - chi2 * E[N - 1]); // boundary conditions
                grid[j][N - 1] = F[N - 1] / (1 - E[N - 1] - h / b(d, t) * (2 * a(d, t) - dBLambda(d, t))); // n - 1 или n - 2, c или d ?
                //grid[j][N - 1] = 0;
                for (int i = N - 2; i >= 0; i--)
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
