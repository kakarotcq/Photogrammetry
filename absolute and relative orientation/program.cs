using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using matrixcal;
using System.IO;

namespace relative_orientation
{
    class Program
    {
        public static double[,] IR(double fi, double w, double k)
        {
            double[,] R = new double[3, 3];
            R[0, 0] = Math.Cos(fi) * Math.Cos(k) - Math.Sin(fi) * Math.Sin(w) * Math.Sin(k); //a1
            R[0, 1] = -Math.Cos(fi) * Math.Sin(k) - Math.Sin(fi) * Math.Sin(w) * Math.Cos(k);//a2
            R[0, 2] = -Math.Sin(fi) * Math.Cos(w);//a3
            R[1, 0] = Math.Cos(w) * Math.Sin(k);//b1
            R[1, 1] = Math.Cos(w) * Math.Cos(k);//b2
            R[1, 2] = -Math.Sin(w);//b3
            R[2, 0] = Math.Sin(fi) * Math.Cos(k) + Math.Cos(fi) * Math.Sin(w) * Math.Sin(k);//c1
            R[2, 1] = -Math.Sin(fi) * Math.Sin(k) + Math.Cos(fi) * Math.Sin(w) * Math.Cos(k);//c2
            R[2, 2] = Math.Cos(fi) * Math.Cos(w);//c3

            return R;
        }
        public static double[,] forward(double phi1, double w1, double k1, double phi2, double w2, double k2, double[,] c1, double[,] c2,double Bx,double By,double Bz,int n,double x0,double y0)
        {
            double f = 40.9349;
            double[,] R1 = new double[3, 3];
            double[,] R2 = new double[3, 3];
            double[,] X = new double[3 * n,1];
            R1 = IR(phi1, w1, k1);
            R2 = IR(phi2, w2, k2);
            double Xs, Ys, Zs;
            Xs = Ys = Zs = 0;
            double Xsb = Bx;double Ysb = By;double Zsb = Bz;
            for (int i = 0; i < n; i++)
            {

                double[,] A = new double[4, 3];
                double[,] l = new double[4, 1];
                double[,] x = new double[3, 1];
                A[0, 0] = f * R1[0, 0] + (c1[i * 2, 0] - x0) * R1[0, 2];
                A[0, 1] = f * R1[1, 0] + (c1[i * 2, 0] - x0) * R1[1, 2];
                A[0, 2] = f * R1[2, 0] + (c1[i * 2, 0] - x0) * R1[2, 2];
                A[1, 0] = f * R1[0, 1] + (c1[i * 2 + 1, 0] - y0) * R1[0, 2];
                A[1, 1] = f * R1[1, 1] + (c1[i * 2 + 1, 0] - y0) * R1[1, 2];
                A[1, 2] = f * R1[2, 1] + (c1[i * 2 + 1, 0] - y0) * R1[2, 2];
                l[0, 0] = f * R1[0, 0] * Xs + f * R1[1, 0] * Ys + f * R1[2, 0] * Zs + (c1[i * 2, 0] - x0) * R1[0, 2] * Xs + (c1[i * 2, 0] - x0) * R1[1, 2] * Ys + (c1[i * 2, 0] - x0) * R1[2, 2] * Zs;
                l[1, 0] = f * R1[0, 1] * Xs + f * R1[1, 1] * Ys + f * R1[2, 1] * Zs + (c1[i * 2+1, 0] - y0) * R1[0, 2] * Xs + (c1[i * 2+1, 0] - y0) * R1[1, 2] * Ys + (c1[i * 2+1, 0] - y0) * R1[2, 2] * Zs;

                A[2, 0] = f * R2[0, 0] + (c2[i * 2, 0] - x0) * R2[0, 2];
                A[2, 1] = f * R2[1, 0] + (c2[i * 2, 0] - x0) * R2[1, 2];
                A[2, 2] = f * R2[2, 0] + (c2[i * 2, 0] - x0) * R2[2, 2];
                A[3, 0] = f * R2[0, 1] + (c2[i * 2 + 1, 0] - y0) * R2[0, 2];
                A[3, 1] = f * R2[1, 1] + (c2[i * 2 + 1, 0] - y0) * R2[1, 2];
                A[3, 2] = f * R2[2, 1] + (c2[i * 2 + 1, 0] - y0) * R2[2, 2];
                l[2, 0] = f * R2[0, 0] * Xsb + f * R2[1, 0] * Ysb + f * R2[2, 0] * Zsb + (c2[i * 2, 0] - x0) * R2[0, 2] * Xsb + (c2[i * 2, 0] - x0) * R2[1, 2] * Ysb + (c2[i * 2, 0] - x0) * R2[2, 2] * Zsb;
                l[3, 0] = f * R2[0, 1] * Xsb + f * R2[1, 1] * Ysb + f * R2[2, 1] * Zsb + (c2[i * 2 + 1, 0] - y0) * R2[0, 2] * Xsb + (c2[i * 2 + 1, 0] - y0) * R2[1, 2] * Ysb + (c2[i * 2 + 1, 0] - y0) * R2[2, 2] * Zsb;
                x= Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), l);
                X[i * 3,0] = x[0, 0];X[i * 3 + 1,0] = x[1, 0];X[i * 3 + 2,0] = x[2, 0];
            }

            /*for (int i = 0; i < n; i++)
            {
                double[,] tmp1 = new double[3, 1];
                double[,] tmp2 = new double[3, 1];
                tmp1[0, 0] = c1[i * 2,0];
                tmp1[1, 0] = c1[i * 2 + 1,0];
                tmp1[2, 0] = -f;
                tmp1 = Matrix.MultiplyMatrix(R1, tmp1);

                tmp2[0, 0] = c2[i * 2,0];
                tmp2[1, 0] = c2[i * 2 + 1,0];
                tmp2[2, 0] = -f;
                tmp2 = Matrix.MultiplyMatrix(R2, tmp2);
                double N = (Bx * tmp2[2, 0] - Bz * tmp2[0, 0]) / (tmp1[0, 0] * tmp2[2, 0] - tmp1[2, 0] * tmp2[0, 0]);
                double N1 = (Bx * tmp1[2, 0] - Bz * tmp1[0, 0]) / (tmp1[0, 0] * tmp2[2, 0] - tmp1[2, 0] * tmp2[0, 0]);
               // Console.Write("{0},{1},{2}\n", Xs + Bx + N1 * tmp2[0, 0], Ys + (N * tmp1[1, 0] + N1 * tmp2[1, 0] + By) / 2, Zs + Bz + N1 * tmp2[2, 0]);
                X[i * 3,0] = Xs + Bx + N1 * tmp2[0, 0];X[i * 3 + 1,0] = Ys + (N * tmp1[1, 0] + N1 * tmp2[1, 0] + By) / 2;X[i * 3 + 2,0] = Zs + Bz + N1 * tmp2[2, 0];
            }*/
            return X;
        }
       /* public static void Bundleadjust(double[,] Po, double[,] c1, double[,] c2,int n,int N,double x0,double y0,double f) //n控制点个数
        {
            double[,] A = new double[N * 4, 12];
            double[,] B = new double[N * 4, 3];
            double[,] l = new double[N * 4, 1];
            double Xs1, Ys1, Zs1; Xs1 = 756.0875; Ys1 = -131.6018; Zs1 = -15.0643;
            double Xs2,Ys2,Zs2; Xs2 = 3081.0581; Ys2 = -135.7630; Zs2 = 68.4453;
            double phi1, w1, k1;phi1 = 0.225967; w1 = 0.112503; k1 = -0.081294;
            double phi2, w2, k2;phi2 = -0.020207; w2 = 0.052306; k2 = -0.078815;
            double[,] R1 = IR(phi1,w1,k1);
            double[,] R2 = IR(phi2,w2,k2);
            for (int i = 0; i < N; i++) //n为控制点个数
            {
                if (i >= n)
                {
                    Po[i * 3, 0] = 0; Po[i * 3 + 1, 0] = 0; Po[i * 3 + 2, 0] = 0;
                }
                double Z = R1[2, 0] * (Po[i*3, 0] - Xs1) + R1[2, 1] * (Po[i*3+1, 0] - Ys1) + R1[2, 2] * (Po[i*3+2, 0] - Zs1);
                A[i * 4, 0] = (R1[0, 0] * f + R1[2, 0] * (c1[i * 2, 0] - x0)) / Z;//a11
                A[i * 4, 1] = (R1[0, 1] * f + R1[2, 1] * (c1[i * 2, 0] - x0)) / Z;//a12
                A[i * 4, 2] = (R1[0, 2] * f + R1[2, 2] * (c1[i * 2, 0] - x0)) / Z;//a13
                A[i * 4, 3] = (c1[i * 2 + 1, 0] - y0) * Math.Sin(w1) - ((c1[i * 2, 0] - x0) / f * ((c1[i * 2, 0] - x0) * Math.Cos(k1) - (c1[i * 2 + 1, 0] - y0) * Math.Sin(k1)) + f * Math.Cos(k1)) * Math.Cos(w1);//a14
                A[i * 4, 4] = -f * Math.Sin(k1) - (c1[i * 2, 0] - x0) / f * ((c1[i * 2, 0] - x0) * Math.Sin(k1) + (c1[i * 2 + 1, 0] - y0) * Math.Cos(k1));//a15
                A[i * 4, 5] = (c1[i * 2 + 1, 0] - y0);//a16
                A[i * 4, 6] = A[i * 4, 7] = A[i * 4, 8] = A[i * 4, 9] = A[i * 4, 10] = A[i * 4, 11] = 0;
                
                 

                A[i * 4 + 1, 0] = (R1[1, 0] * f + R1[2, 0] * (c1[i * 2 + 1, 0] - y0)) / Z;//a21
                A[i * 4 + 1, 1] = (R1[1, 1] * f + R1[2, 1] * (c1[i * 2 + 1, 0] - y0)) / Z;//a22
                A[i * 4 + 1, 2] = (R1[1, 2] * f + R1[2, 2] * (c1[i * 2 + 1, 0] - y0)) / Z;//a23
                A[i * 4 + 1, 3] = -(c1[i * 2, 0] - x0) * Math.Sin(w1) - ((c1[i * 2 + 1, 0] - y0) / f * ((c1[i * 2, 0] - x0) * Math.Cos(k1) - (c1[i * 2 + 1, 0] - y0) * Math.Sin(k1)) - f * Math.Sin(k1)) * Math.Cos(w1);//a24
                A[i * 4 + 1, 4] = -f * Math.Cos(k1) - (c1[i * 2 + 1, 0] - y0) / f * ((c1[i * 2, 0] - x0) * Math.Sin(k1) + (c1[i * 2 + 1, 0] - y0) * Math.Cos(k1));//a25
                A[i * 4 + 1, 5] = -(c1[i * 2, 0] - x0);//a26
                A[i * 4 + 1, 6] = A[i * 4 + 1, 7] = A[i * 4 + 1, 8] = A[i * 4 + 1, 9] = A[i * 4 + 1, 10] = A[i * 4 + 1, 11] = 0;
                
                double cx, cy;//像点坐标近似值
                cx = x0 - f * (R1[0, 0] * (Po[i*3,0] - Xs1) + R1[0, 1] * (Po[i*3+1, 0] - Ys1) + R1[0, 2] * (Po[i*3+2, 0] - Zs1)) / (R1[2, 0] * (Po[i, 0] - Xs1) + R1[2, 1] * (Po[i*3+1, 0] - Ys1) + R1[2, 2] * (Po[i*3+2, 0] - Zs1));
                cy = y0 - f * (R1[1, 0] * (Po[i*3, 0] - Xs1) + R1[1, 1] * (Po[i*3+1, 0] - Ys1) + R1[1, 2] * (Po[i*3+2, 0] - Zs1)) / (R1[2, 0] * (Po[i, 0] - Xs1) + R1[2, 1] * (Po[i*3+1, 0] - Ys1) + R1[2, 2] * (Po[i*3+2, 0] - Zs1));
                l[4 * i, 0] = c1[i * 2, 0] - cx;
                l[4 * i + 1, 0] = c1[i * 2 + 1, 0] - cy;

                Z = R2[2, 0] * (Po[i*3, 0] - Xs2) + R2[2, 1] * (Po[i*3+1, 0] - Ys2) + R2[2, 2] * (Po[i*3+2, 0] - Zs2);
                A[i * 4 + 2, 6] = (R2[0, 0] * f + R2[2, 0] * (c2[i * 2, 0] - x0)) / Z;//a11
                A[i * 4 + 2, 7] = (R2[0, 1] * f + R2[2, 1] * (c2[i * 2, 0] - x0)) / Z;//a12
                A[i * 4 + 2, 8] = (R2[0, 2] * f + R2[2, 2] * (c2[i * 2, 0] - x0)) / Z;//a13
                A[i * 4 + 2, 9] = (c2[i * 2 + 1, 0] - y0) * Math.Sin(w2) - ((c2[i * 2, 0] - x0) / f * ((c2[i * 2, 0] - x0) * Math.Cos(k2) - (c2[i * 2 + 1, 0] - y0) * Math.Sin(k2)) + f * Math.Cos(k2)) * Math.Cos(w2);//a14
                A[i * 4 + 2, 10] = -f * Math.Sin(k2) - (c2[i * 2, 0] - x0) / f * ((c2[i * 2, 0] - x0) * Math.Sin(k2) + (c2[i * 2 + 1, 0] - y0) * Math.Cos(k2));//a15
                A[i * 4 + 2, 11] = (c2[i * 2 + 1, 0] - y0);//a16
                A[i * 4 + 2, 0] = A[i * 4 + 2, 1] = A[i * 4 + 2, 2] = A[i * 4 + 2, 3] = A[i * 4 + 2, 4] = A[i * 4 + 2, 5] = 0;
                

                A[i * 4 + 3, 6] = (R2[1, 0] * f + R2[2, 0] * (c2[i * 2 + 1, 0] - y0)) / Z;//a21
                A[i * 4 + 3, 7] = (R2[1, 1] * f + R2[2, 1] * (c2[i * 2 + 1, 0] - y0)) / Z;//a22
                A[i * 4 + 3, 8] = (R2[1, 2] * f + R2[2, 2] * (c2[i * 2 + 1, 0] - y0)) / Z;//a23
                A[i * 4 + 3, 9] = -(c2[i * 2, 0] - x0) * Math.Sin(w2) - ((c2[i * 2 + 1, 0] - y0) / f * ((c2[i * 2, 0] - x0) * Math.Cos(k2) - (c2[i * 2 + 1, 0] - y0) * Math.Sin(k2)) - f * Math.Sin(k2)) * Math.Cos(w2);//a24
                A[i * 4 + 3, 10] = -f * Math.Cos(k2) - (c2[i * 2 + 1, 0] - y0) / f * ((c2[i * 2, 0] - x0) * Math.Sin(k2) + (c2[i * 2 + 1, 0] - y0) * Math.Cos(k2));//a25
                A[i * 4 + 3, 11] = -(c2[i * 2, 0] - x0);//a26
                A[i * 4 + 3, 0] = A[i * 4 + 3, 1] = A[i * 4 + 3, 2] = A[i * 4 + 3, 3] = A[i * 4 + 3, 4] = A[i * 4 + 3, 5] = 0;


                if (i >= n)
                {
                    B[i * 4, 0] = -A[i * 4, 0]; B[i * 4, 1] = -A[i * 4, 1]; B[i * 4, 2] = -A[i * 4, 2];
                    B[i * 4 + 1, 0] = -A[i * 4 + 1, 0]; B[i * 4 + 1, 1] = -A[i * 4 + 1, 2]; B[i * 4 + 1, 2] = -A[i * 4 + 1, 2];
                    B[i * 4 + 2, 0] = -A[i * 4 + 2, 6]; B[i * 4 + 2, 1] = -A[i * 4 + 2, 7]; B[i * 4 + 2, 2] = -A[i * 4 + 2, 8];
                    B[i * 4 + 3, 0] = -A[i * 4 + 3, 6]; B[i * 4 + 3, 1] = -A[i * 4 + 3, 7]; B[i * 4 + 3, 2] = -A[i * 4 + 3, 8];
                }
     
                cx = x0 - f * (R2[0, 0] * (Po[i*3, 0] - Xs2) + R2[0, 1] * (Po[i*3+1, 0] - Ys2) + R2[0, 2] * (Po[i*3+2, 0] - Zs2)) / (R2[2, 0] * (Po[i*3, 0] - Xs2) + R2[2, 1] * (Po[i*3+1, 0] - Ys2) + R2[2, 2] * (Po[i*3+2, 0] - Zs2));
                cy = y0 - f * (R2[1, 0] * (Po[i*3, 0] - Xs2) + R2[1, 1] * (Po[i*3+1, 0] - Ys2) + R2[1, 2] * (Po[i*3+2, 0] - Zs2)) / (R2[2, 0] * (Po[i*3, 0] - Xs2) + R2[2, 1] * (Po[i*3+1, 0] - Ys2) + R2[2, 2] * (Po[i*3+2, 0] - Zs2));
                l[4 * i + 2, 0] = c2[i * 2, 0] - cx;
                l[4 * i + 3, 0] = c2[i * 2 + 1, 0] - cy;
            }
            double[,] N11 = Matrix.MultiplyMatrix(Matrix.Transpose(A), A);

            double[,] N12 = new double[12, 3 * N]; //N列12*3矩阵
            double[,] AT = Matrix.Transpose(A);
            for (int j = n; j < N; j++)             //i行 j列 前n列0
            {
                double[,] tmp1 = new double[12, 4];
                double[,] tmp2 = new double[4, 3];
                for (int k = 0; k < 12; k++)
                    for (int m = 0; m < 4; m++)
                        tmp1[k, m] = AT[k, j * 4 + m];
                for (int k = 0; k < 4; k++)
                    for (int m = 0; m < 3; m++)
                        tmp2[k, m] = B[j * 4 + k, m];
                double[,] tmp = Matrix.MultiplyMatrix(tmp1, tmp2);
                for (int k = 0; k < 12; k++)
                    for (int m = 0; m < 3; m++)
                        N12[k, j * 3 + m] = tmp[k, m];
            }
            //double[,] N22 = Matrix.MultiplyMatrix(Matrix.Transpose(B), B);
            double[,] u1 = Matrix.MultiplyMatrix(Matrix.Transpose(A), l);
            double[,] u2 = new double[3 * N, 1];
            for (int i = n; i < N ; i++)//N行3*1矩阵
            {
                double[,] tmp1 = new double[3, 4];
                double[,] tmp2 = new double[4, 1];
                for (int k = 0; k < 3; k++)
                    for (int m = 0; m < 4; m++)
                        tmp1[k, m] = B[i * 4 + m, k];
                for (int k = 0; k < 4; k++)
                    tmp2[k, 0] = l[i * 4 + k, 0];
                double[,] tmp = Matrix.MultiplyMatrix(tmp1, tmp2);
                for (int j = 0; j < 3; j++)
                    u2[i * 3 + j, 0] = tmp[j, 0];
            }
            double[,] t = new double[12, 1];
            double[,] X = new double[3 * N, 1];
            double[,] t1 = new double[12, 3*N];//N12N22-1
            for (int i = n; i < N ; i++) //N个12*3的矩阵
            {
                double[,] tmp1 = new double[4, 3];
                double[,] tmp2 = new double[3, 4];
                double[,] tmp3 = new double[12, 3];
                for (int k = 0; k < 4; k++)
                    for (int m = 0; m < 3; m++)
                    {
                        tmp1[k, m] = B[4 * i + k,  m];
                        tmp2[m, k] = B[4 * i + k, m];
                    }
                for (int k = 0; k < 12; k++)
                    for (int m = 0; m < 3; m++)
                        tmp3[k, m] = N12[k, i * 3 + m];
                tmp3 = Matrix.MultiplyMatrix(tmp3, Matrix.Athwart(Matrix.MultiplyMatrix(tmp1, tmp2)));
                for (int k = 0; k < 12; k++)
                    for (int m = 0; m < 3; m++)
                        t1[k, i * 3 + m] = tmp3[k, m];
            }
            t = Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.SubMatrix(N11, Matrix.MultiplyMatrix(t1, Matrix.Transpose(N12)))), Matrix.SubMatrix(u1, Matrix.MultiplyMatrix(t1, u2)));
             double[,]t2= new double[3 * N, 1];
            //t2 = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Transpose(N12), Matrix.Athwart(N11)), N12);
            t2 = Matrix.SubMatrix(u2, Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Transpose(N12), Matrix.Athwart(N11)), u1));
            double judge = 0.1 / 60 * 2 * Math.PI / 360;
            while (Math.Abs(X[0, 0]) > 0.1 || Math.Abs(X[1, 0]) > 0.1 || Math.Abs(X[2, 0]) > 0.1 || Math.Abs(X[3, 0]) > judge || Math.Abs(X[4, 0]) > judge || Math.Abs(X[5, 0]) > judge || Math.Abs(X[6, 0]) > 0.1 || Math.Abs(X[7, 0]) > 0.1 || Math.Abs(X[8, 0]) > 0.1 || Math.Abs(X[9, 0]) > judge || Math.Abs(X[10, 0]) > judge || Math.Abs(X[11, 0]) > judge)
            {

            }

                Console.ReadLine();
        }*/

        static void Main(string[] args)
        {
            int n;//点数
            double x0, y0, f, m, K1, K2, p1, p2, a, B;//畸变系数
            StreamReader reader = new StreamReader("bundleadjustment_SCBA_Camera_Result.scbacmr");
            reader.ReadLine();//过掉第一行
            string line = reader.ReadLine();
            reader.Close();
            string[] tmp = line.Split('\t');
            //int j = 0;
            x0 = double.Parse(tmp[0]);
            y0 = double.Parse(tmp[1]);
            f = double.Parse(tmp[2]);
            //f *= 0.001;
            m = double.Parse(tmp[5]);
            K1 = double.Parse(tmp[6]);
            K2 = double.Parse(tmp[7]);
            p1 = double.Parse(tmp[8]);
            p2 = double.Parse(tmp[9]);
            a = double.Parse(tmp[10]);
            B = double.Parse(tmp[11]);
            reader = new StreamReader("bundleadjustment_SCBA_Point_Result.scbapts");
            line = reader.ReadLine();
            tmp = line.Split('\t');
            n = int.Parse(tmp[0]);
            double[,] C1 = new double[2 * n,1];
            double[,] C2 = new double[2 * n,1];
            double[,] X = new double[5, 1];//五个gai'zheng'shu
            double[,] A = new double[n, 5];
            double[,] P = new double[n * 3, 1];//控制点物方空间坐标
                                               // double[,] B = new double[n, 5];
            double[,] L = new double[n, 1];
            double[,] V = new double[n, 1];
            double[] Pname = new double[n];
            reader.ReadLine();
            for (int i = 0; i < n; i++)//将点位坐标放入数组中
            {
                line = reader.ReadLine();
                tmp = line.Split('\t');
                Pname[i] = double.Parse(tmp[0]);
                P[i * 3, 0] = double.Parse(tmp[1]);
                P[i * 3 + 1, 0] = double.Parse(tmp[2]);
                P[i * 3 + 2, 0] = double.Parse(tmp[3]);
                reader.ReadLine();
                line = reader.ReadLine();
                tmp = line.Split('\t');
                C1[i * 2,0] = double.Parse(tmp[2]);
                C1[i * 2 + 1,0] = double.Parse(tmp[3]);
                double r = Math.Sqrt(Math.Pow(C1[i * 2,0] - x0, 2) + Math.Pow(C1[i * 2 + 1,0] - y0, 2));
                double dx = (C1[i * 2,0] - x0) * (K1 * r * r + K2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(C1[i * 2,0] - x0, 2)) + 2 * p2 * (C1[i * 2,0] - x0) * (C1[i * 2 + 1,0] - y0) + a * (C1[i * 2,0] - x0) + B * (C1[i * 2 + 1,0] - y0);
                double dy = (C1[i * 2 + 1,0] - y0) * (K1 * r * r + K2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(C1[i * 2 + 1,0] - y0, 2) + 2 * p1 * (C1[i * 2,0] - x0) * (C1[i * 2 + 1,0] - y0));
                C1[i * 2,0] -= dx;
                C1[i * 2 + 1,0] -= dy;
                line = reader.ReadLine();
                tmp = line.Split('\t');

                C2[i * 2,0] = double.Parse(tmp[2]);
                C2[i * 2 + 1,0] = double.Parse(tmp[3]);
                r = Math.Sqrt(Math.Pow(C2[i * 2,0] - x0, 2) + Math.Pow(C2[i * 2 + 1,0] - y0, 2));
                dx = (C2[i * 2,0] - x0) * (K1 * r * r + K2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(C2[i * 2,0] - x0, 2)) + 2 * p2 * (C2[i * 2,0] - x0) * (C2[i * 2 + 1,0] - y0) + a * (C2[i * 2,0] - x0) + B * (C2[i * 2 + 1,0] - y0);
                dy = (C2[i * 2 + 1,0] - y0) * (K1 * r * r + K2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(C2[i * 2 + 1,0] - y0, 2) + 2 * p1 * (C2[i * 2,0] - x0) * (C2[i * 2 + 1,0] - y0));
                C2[i * 2,0] -= dx;
                C2[i * 2 + 1,0] -= dy;
            }
            reader.Close();
            //光束法
          // Bundleadjust(P, C1, C2, 4, n, x0, y0, f);

            double Bx = 500;
            const double phi1 = 0, w1 = 0, k1 = 0;//第一个相机的已知外方位元素0.235967	0.102503	-0.041294
            double By, Bz, phi2 = phi1, w2 = w1, k2 = k1, u, v;
            By = Bz = u = v = 0;
            StreamWriter Write = new StreamWriter("相对定向迭代过程");
            Write.WriteLine("第一张相片的三个角度值均为0，Bz=500");
            Write.WriteLine("{0}    {1} {2} {3} {4}", phi2, w2, k2, By, Bz);
            By = Bx*u; Bz =Bx*v;
            for (int i = 0; i < n; i++)
            {

                double[,] X1 = new double[3, 1];//第一个相片的像点的像空间辅助坐标
                double[,] X2 = new double[3, 1];//第二张相片的像点的像空间辅助坐标
               // double[,] mat = new double[3, 3];//算行列式用的辅助矩阵
                double[,] tmppp = new double[3, 1];
                double[,] R1 = new double[3, 3];
                double[,] R2 = new double[3, 3];
                double q, N, N1;
                tmppp[0, 0] = C1[i * 2,0];
                tmppp[1, 0] = C1[i * 2 + 1,0];
                tmppp[2, 0] = -40.9349;
                R1 = IR(phi1, w1, k1);
                X1 = Matrix.MultiplyMatrix(R1, tmppp);
                tmppp[0, 0] = C2[i * 2,0];
                tmppp[1, 0] = C2[i * 2 + 1,0];
                tmppp[2, 0] = -40.9349;
                R2 = IR(phi2, w2, k2);
                X2 = Matrix.MultiplyMatrix(R2, tmppp);            
                N = (Bx * X2[2, 0] - Bz * X2[0, 0]) / (X1[0, 0] * X2[2, 0] - X1[2, 0] * X2[0, 0]);
                N1 = (Bx * X1[2, 0] - Bz * X1[0, 0]) / (X1[0, 0] * X2[2, 0] - X1[2, 0] * X2[0, 0]);
                q = N * X1[1, 0] - N1 * X2[1, 0] - By; ;
                A[i, 0] = -X2[0, 0] * X2[1, 0] / X2[2, 0] * N1;
                A[i, 1] = -(X2[2, 0] + Math.Pow(X2[1, 0], 2) / X2[2, 0]) * N1;
                A[i, 2] = X2[0, 0] * N1;
                A[i, 3] = Bx; A[i, 4] = -X2[1, 0] / X2[2, 0]*Bx ;
                L[i, 0] = q;

            }
            //X = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Transpose(B), Matrix.Athwart(Matrix.MultiplyMatrix(A, Matrix.Transpose(A)))), B)), Matrix.Transpose(B)), Matrix.Athwart(Matrix.MultiplyMatrix(A, Matrix.Transpose(A)))), L);
            X = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), L);
            //By += X[0,0]; Bz += X[1, 0];phi2 += X[2, 0];w2 += X[3, 0];k2 += X[4,0];
            Write.WriteLine("改正值 {0}    {1} {2} {3} {4}", X[0, 0], X[1, 0], X[2, 0], X[3, 0], X[4, 0]);
            phi2 += X[0, 0]; w2 += X[1, 0]; k2 += X[2, 0]; u += X[3, 0]; v += X[4, 0];
            V = Matrix.SubMatrix(Matrix.MultiplyMatrix(A, X), L);
            // L = Matrix.AddMatrix(L, V);



            while (Math.Abs(X[0, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(X[1, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(X[2, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(X[3, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(X[4, 0]) > 0.3 * Math.Pow(10, -4))
            {
                By = Bx * Math.Tan(u); Bz = Math.Sqrt(Bx * Bx + By * By) * Math.Tan(v);
                Write.WriteLine("{0}    {1} {2} {3} {4}", phi2, w2, k2, By, Bz);
                for (int i = 0; i < n; i++)
                {

                    double[,] X1 = new double[3, 1];//第一个相片的像点的像空间辅助坐标
                    double[,] X2 = new double[3, 1];//第二张相片的像点的像空间辅助坐标
                    //double[,] mat = new double[3, 3];//算行列式用的辅助矩阵
                    double[,] tmppp = new double[3, 1];
                    double[,] R1 = new double[3, 3];
                    double[,] R2 = new double[3, 3];
                    double q, N, N1;
                    tmppp[0, 0] = C1[i * 2,0];
                    tmppp[1, 0] = C1[i * 2 + 1,0];
                    tmppp[2, 0] = -40.9349;
                    R1 = IR(phi1, w1, k1);
                    X1 = Matrix.MultiplyMatrix(R1, tmppp);
                    tmppp[0, 0] = C2[i * 2,0];
                    tmppp[1, 0] = C2[i * 2 + 1,0];
                    tmppp[2, 0] = -40.9349;
                    R2 = IR(phi2, w2, k2);
                    X2 = Matrix.MultiplyMatrix(R2, tmppp);
                 

                    //By = Bx * Math.Tan(u);Bz = Math.Sqrt(Bx * Bx + By * By ) * Math.Tan(v);
                    N = (Bx * X2[2, 0] - Bz * X2[0, 0]) / (X1[0, 0] * X2[2, 0] - X1[2, 0] * X2[0, 0]);
                    N1 = (Bx * X1[2, 0] - Bz * X1[0, 0]) / (X1[0, 0] * X2[2, 0] - X1[2, 0] * X2[0, 0]);
                    q = N * X1[1, 0] - N1 * X2[1, 0] - By;
                    A[i, 0] = -X2[0, 0] * X2[1, 0] / X2[2, 0] * N1;
                    A[i, 1] = -(X2[2, 0] + Math.Pow(X2[1, 0], 2) / X2[2, 0]) * N1;
                    A[i, 2] = X2[0, 0] * N1;
                    A[i, 3] = Bx; A[i, 4] = -X2[1, 0] / X2[2, 0]*Bx ;
                    L[i, 0] = q;
                    Console.Write("{0}\n", q);

                }
                //X = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Transpose(B), Matrix.Athwart(Matrix.MultiplyMatrix(A, Matrix.Transpose(A)))), B)), Matrix.Transpose(B)), Matrix.Athwart(Matrix.MultiplyMatrix(A, Matrix.Transpose(A)))), L);
                //By += X[0, 0]; Bz += X[1, 0]; phi2 += X[2, 0]; w2 += X[3, 0]; k2 += X[4, 0];
                X = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), L);
                //By += X[0,0]; Bz += X[1, 0];phi2 += X[2, 0];w2 += X[3, 0];k2 += X[4,0];
                phi2 += X[0, 0]; w2 += X[1, 0]; k2 += X[2, 0]; u += X[3, 0]; v += X[4, 0];
                Write.WriteLine("改正值 {0}    {1} {2} {3} {4}", X[0, 0], X[1, 0], X[2, 0], X[3, 0], X[4, 0]);
                V = Matrix.SubMatrix(Matrix.MultiplyMatrix(A, X), L);
                //L = Matrix.AddMatrix(L, V);
            }
            By = Bx * Math.Tan(u); Bz = Math.Sqrt(Bx * Bx + By * By) * Math.Tan(v);
            Write.WriteLine("{0}    {1} {2} {3} {4}", phi2, w2, k2, By, Bz);
            Write.Close();
            double[,] Q = Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A));
            double[,] t = Matrix.MultiplyMatrix(Matrix.Transpose(V), V);
            double m0 = Math.Sqrt(t[0, 0] / (2 * n - 6));
            //By = Bx * Math.Tan(u); Bz = Math.Sqrt(Bx * Bx + By * By) * Math.Tan(v);
            StreamWriter writer = new StreamWriter("相对定向结果及精度");
            writer.WriteLine("phi2={0}|{1}", phi2, Math.Sqrt(Q[0, 0]) * m0);
            writer.WriteLine("w2={0}|{1}", w2, Math.Sqrt(Q[1, 1]) * m0);
            writer.WriteLine("k2={0}|{1}", k2, Math.Sqrt(Q[2, 2]) * m0);
            writer.WriteLine("Bx={0}", Bx);
            writer.WriteLine("By={0}|{1}", By, Math.Sqrt(Q[3, 3]) * m0);
            writer.WriteLine("Bz={0}|{1}", Bz, Math.Sqrt(Q[4, 4]) * m0);
            writer.Close();

            //前方交会
            By = 9.58339;Bz = -100.899;phi2 = -0.246099;w2 = -0.0472591;k2 = -0.0737229;
            double[,] XX = new double[3 * n,1];//前方交会后的坐标
            XX = forward(phi1, w1, k1, phi2, w2, k2, C1, C2, Bx, By, Bz, n, x0, y0);
           
            
            //绝对定向

            double dX=0, dY=0, dZ=0, na=1, PHI=0, W=0, K=0;
            double[,] DX = new double[7, 1];
            double[,] V1 = new double[3*n, 1];
            double[,] L1 = new double[3*n, 1];
            double[,] A1 = new double[3 * n, 7];

            writer = new StreamWriter("绝对定向迭代过程");
            writer.WriteLine("dX dY dZ na phi w k ");
            writer.WriteLine("{0} {1} {2} {3} {4} {5} {6}", dX, dY, dZ, na, PHI, W, K);
            for (int i = 0; i < n; i++)
            {
                double[,] ttmp = new double[3, 1];
                ttmp[0, 0] = XX[i * 3,0];
                ttmp[1, 0] = XX[i*3+1,0];
                ttmp[2, 0] = XX[i*3+2,0];
                ttmp = Matrix.MultiplyMatrix(IR(PHI, W, K), ttmp);
                A1[i * 3, 0] = 1; A1[i * 3, 1] = 0;A1[i*3,2]=0 ; A1[i * 3, 3] = ttmp[0, 0]; A1[i * 3, 4] = -na * ttmp[2, 0];A1[i * 3, 5] = -na * ttmp[1, 0] * Math.Sin(PHI);A1[i * 3, 6] = -na * ttmp[1, 0] * Math.Cos(PHI) * Math.Cos(W) - na * ttmp[2, 0] * Math.Sin(W);
                A1[i * 3 + 1, 0] = 0; A1[i * 3 + 1, 1] = 1; A1[i*3+1,2]=0; A1[i * 3 + 1, 3] = ttmp[1, 0]; A1[i * 3 + 1, 4] = 0;A1[i * 3 + 1, 5] = na * ttmp[0, 0] * Math.Sin(PHI) - na * ttmp[2, 0] * Math.Cos(PHI);A1[i*3+1,6]=na*ttmp[0,0]*Math.Cos(PHI)*Math.Cos(W)+na*ttmp[2,0]*Math.Sin(PHI)*Math.Cos(W);
                A1[i * 3 + 2, 0] = 0; A1[i * 3 + 2, 1] = 0;A1[i*3+2,2]=1 ; A1[i * 3 + 2, 3] = ttmp[2, 0]; A1[i * 3 + 2, 4] = na * ttmp[0, 0];A1[i * 3 + 2,5] = na * ttmp[1, 0] * Math.Cos(PHI);A1[i * 3 + 2, 6] = na * ttmp[0, 0] * Math.Sin(W) - na * ttmp[1, 0] * Math.Sin(PHI) * Math.Cos(W);
                L1[i * 3, 0] = P[i * 3, 0] - dX - na * ttmp[0, 0];
                L1[i * 3 + 1, 0] = P[i * 3 + 1, 0] - dY - na * ttmp[1, 0];
                L1[i * 3 + 2, 0] = P[i * 3 + 2, 0] - dZ - na * ttmp[2, 0];
            }
            DX = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A1), A1)), Matrix.Transpose(A1)), L1);
            dX +=DX[0, 0];dY += DX[1, 0];dZ += DX[2, 0]; na = na+ DX[3, 0]; PHI += DX[4, 0];W += DX[5, 0];K += DX[6, 0];
            writer.WriteLine("改正值：{0} {1} {2} {3} {4} {5} {6}", DX[0, 0], DX[1, 0], DX[2, 0], DX[3, 0], DX[4, 0], DX[5, 0], DX[6, 0]);
            writer.WriteLine("{0} {1} {2} {3} {4} {5} {6}", dX, dY, dZ, na, PHI, W, K);
            V1 = Matrix.SubMatrix(Matrix.MultiplyMatrix(A1, DX), L1);
           // XX = Matrix.AddMatrix(XX, V1);
            while (Math.Abs(DX[0, 0]) > 0.1 || Math.Abs(DX[1, 0]) > 0.1 || Math.Abs(DX[2, 0]) > 0.1 || Math.Abs(DX[3, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(DX[4, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(DX[5, 0]) > 0.3 * Math.Pow(10, -4) || Math.Abs(DX[6, 0]) > 0.3 * Math.Pow(10, -4))
            {
                for (int i = 0; i < n; i++)
                {
                    double[,] ttmp = new double[3, 1];
                    ttmp[0, 0] = XX[i * 3,0];
                    ttmp[1, 0] = XX[i * 3 + 1,0];
                    ttmp[2, 0] = XX[i * 3 + 2,0];
                    ttmp = Matrix.MultiplyMatrix(IR(PHI, W, K), ttmp);
                    A1[i * 3, 0] = 1; A1[i * 3, 1] = 0; A1[i * 3, 2] = 0; A1[i * 3, 3] = ttmp[0, 0]; A1[i * 3, 4] = -na * ttmp[2, 0]; A1[i * 3, 5] = -na * ttmp[1, 0] * Math.Sin(PHI); A1[i * 3, 6] = -na * ttmp[1, 0] * Math.Cos(PHI) * Math.Cos(W) - na * ttmp[2, 0] * Math.Sin(W);
                    A1[i * 3 + 1, 0] = 0; A1[i * 3 + 1, 1] = 1; A1[i * 3 + 1, 2] = 0; A1[i * 3 + 1, 3] = ttmp[1, 0]; A1[i * 3 + 1, 4] = 0; A1[i * 3 + 1, 5] = na * ttmp[0, 0] * Math.Sin(PHI) - na * ttmp[2, 0] * Math.Cos(PHI); A1[i * 3 + 1, 6] = na * ttmp[0, 0] * Math.Cos(PHI) * Math.Cos(W) + na * ttmp[2, 0] * Math.Sin(PHI) * Math.Cos(W);
                    A1[i * 3 + 2, 0] = 0; A1[i * 3 + 2, 1] = 0; A1[i * 3 + 2, 2] = 1; A1[i * 3 + 2, 3] = ttmp[2, 0]; A1[i * 3 + 2, 4] = na * ttmp[0, 0]; A1[i * 3 + 2, 5] = na * ttmp[1, 0] * Math.Cos(PHI); A1[i * 3 + 2, 6] = na * ttmp[0, 0] * Math.Sin(W) - na * ttmp[1, 0] * Math.Sin(PHI) * Math.Cos(W);
                    L1[i * 3, 0] = P[i * 3, 0] - dX - na * ttmp[0, 0];
                    L1[i * 3 + 1, 0] = P[i * 3 + 1, 0] - dY - na * ttmp[1, 0];
                    L1[i * 3 + 2, 0] = P[i * 3 + 2, 0] - dZ - na * ttmp[2, 0];
                }
                DX = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A1), A1)), Matrix.Transpose(A1)), L1);
                dX += DX[0, 0]; dY += DX[1, 0]; dZ += DX[2, 0]; na = na  + DX[3, 0]; PHI += DX[4, 0]; W += DX[5, 0]; K += DX[6, 0];
                writer.WriteLine("改正值：{0} {1} {2} {3} {4} {5} {6}\n", DX[0, 0], DX[1, 0], DX[2, 0], DX[3, 0], DX[4, 0], DX[5, 0], DX[6, 0]);
                writer.WriteLine("{0} {1} {2} {3} {4} {5} {6}\n", dX, dY, dZ, na, PHI, W, K);
                V1 = Matrix.SubMatrix(Matrix.MultiplyMatrix(A1, DX), L1);
                //XX = Matrix.AddMatrix(XX, V1);
                //P = Matrix.AddMatrix(P, V1);
            }
            writer.Close();
            writer = new StreamWriter("绝对定向结果及精度 ");
            double[,] Q1 = Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A1), A1));
            double[,] t1 = Matrix.MultiplyMatrix(Matrix.Transpose(V1), V1);
            m0 = Math.Sqrt(t[0, 0] / (2 * n - 6));
            writer.WriteLine("dX={0}|{1}", dX, Math.Sqrt(Q1[0, 0]) * m0);
            writer.WriteLine("dY={0}|{1}", dY, Math.Sqrt(Q1[1, 1]) * m0);
            writer.WriteLine("dZ={0}|{1}", dZ, Math.Sqrt(Q1[2, 2]) * m0);
            writer.WriteLine("na={0}|{1}", na, Math.Sqrt(Q1[3, 3]) * m0);
            writer.WriteLine("phi={0}|{1}", PHI, Math.Sqrt(Q1[4, 4]) * m0);
            writer.WriteLine("w={0}|{1}", W, Math.Sqrt(Q1[5, 5]) * m0);
            writer.WriteLine("k={0}|{1}", K, Math.Sqrt(Q1[6, 6]) * m0);
            writer.Close();

            for (int i = 0; i < n; i++)
            {
                double[,] R = IR(PHI, W, K);
                double[,] ttmp = new double[3, 1];
                ttmp[0, 0] = XX[i * 3,0];
                ttmp[1, 0] = XX[i * 3 + 1,0];
                ttmp[2, 0] = XX[i * 3 + 2,0];
                ttmp =  Matrix.MultiplyMatrix(R, ttmp);
                ttmp[0, 0] = na * ttmp[0, 0] + dX;
                ttmp[1, 0] = na * ttmp[1, 0] + dY;
                ttmp[2, 0] = na * ttmp[2, 0] + dZ;
                Console.Write("点号：{0},{1},{2},{3}\n",Pname[i], ttmp[0, 0], ttmp[1, 0], ttmp[2, 0]);
            }
           Console.ReadLine();
        }
    }
}
