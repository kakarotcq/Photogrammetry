using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using matrixcal;

namespace photogram
{
    class Program
    {
        /* public void pcrec(double x,double y,double x0,double y0,out double dx,out double dy,double k1,double k2,double p1,double p2,double a,double B)
         {
             double r = Math.Sqrt(Math.Pow(x-x0,2)-Math.Pow(y-y0,2));
             dx = (x - x0)*(k1 * r * r + k2 * Math.Pow(r, 4))+p1*(r*r+2*Math.Pow(x-x0,2))+2*p2*(x-x0)*(y-y0)+a*(x-x0)+B*(y-y0);
             dy = (y - y0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(y - y0, 2) + 2 * p1 * (x - x0)*(y - y0));
         }*/
         //对于一个相机输入参数计算其外方位元素
        public static double[,] cal(double x0,double y0,double f,double k1,double k2,double p1,double p2,double a,double B,ref double Xs,ref double Ys,ref double Zs,ref double fi,ref double w,ref double k,int n,double[,]c,double[,]Po, StreamWriter writer)
        {
            double[,] A = new double[n / 2 * 2, 6];
            double[,] l = new double[n / 2 * 2, 1];
            double[,] V = new double[n / 2 * 2, 1];
            double[,] dX = new double[6, 1];//6个外方位元素改正值

            writer.WriteLine("{0}  {1}  {2}  {3}  {4}  {5}", Xs, Ys, Zs, fi, w, k);
            for (int i = 0; i < n / 2; i++)
            {
                double[,] R = IR(fi, w, k);
                double r = Math.Sqrt(Math.Pow(c[i * 2, 0] - x0, 2) + Math.Pow(c[i * 2 + 1, 0] - y0, 2));
                double dx = (c[i * 2, 0] - x0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(c[i * 2, 0] - x0, 2)) + 2 * p2 * (c[i * 2, 0] - x0) * (c[i * 2 + 1, 0] - y0) + a * (c[i * 2, 0] - x0) + B * (c[i * 2 + 1, 0] - y0);
                double dy = (c[i * 2 + 1, 0] - y0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(c[i * 2 + 1, 0] - y0, 2) + 2 * p1 * (c[i * 2, 0] - x0) * (c[i * 2 + 1, 0] - y0));
                //double tx0 = x0, ty0 = y0;//像点中心坐标
                //tx0 += dx; ty0 += dy;//更正像原点
                //CX[i * 2] = tx0; CX[i * 2 + 1] = ty0;
                c[i * 2, 0] -= dx; c[i * 2 + 1, 0] -= dy;
                double Z = R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs);
                A[i * 2, 0] = (R[0, 0] * f + R[2, 0] * (c[i * 2, 0] - x0)) / Z;//a11
                A[i * 2, 1] = (R[0, 1] * f + R[2, 1] * (c[i * 2, 0] - x0)) / Z;//a12
                A[i * 2, 2] = (R[0, 2] * f + R[2, 2] * (c[i * 2, 0] - x0)) / Z;//a13
                A[i * 2, 3] = (c[i * 2 + 1, 0] - y0) * Math.Sin(w) - ((c[i * 2, 0] - x0) / f * ((c[i * 2, 0] - x0) * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) * Math.Sin(k)) + f * Math.Cos(k)) * Math.Cos(w);//a14
                A[i * 2, 4] = -f * Math.Sin(k) - (c[i * 2, 0] - x0) / f * ((c[i * 2, 0] - x0) * Math.Sin(k) + (c[i * 2 + 1, 0] - y0) * Math.Cos(k));//a15
                A[i * 2, 5] = (c[i * 2 + 1, 0] - y0);//a16

                A[i * 2 + 1, 0] = (R[1, 0] * f + R[2, 0] * (c[i * 2 + 1, 0] - y0)) / Z;//a21
                A[i * 2 + 1, 1] = (R[1, 1] * f + R[2, 1] * (c[i * 2 + 1, 0] - y0)) / Z;//a22
                A[i * 2 + 1, 2] = (R[1, 2] * f + R[2, 2] * (c[i * 2 + 1, 0] - y0)) / Z;//a23
                A[i * 2 + 1, 3] = -(c[i * 2, 0] - x0) * Math.Sin(w) - ((c[i * 2 + 1, 0] - y0) / f * ((c[i * 2, 0] - x0) * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) * Math.Sin(k)) - f * Math.Sin(k)) * Math.Cos(w);//a24
                A[i * 2 + 1, 4] = -f * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) / f * ((c[i * 2, 0] - x0) * Math.Sin(k) + (c[i * 2 + 1, 0] - y0) * Math.Cos(k));//a25
                A[i * 2 + 1, 5] = -(c[i * 2, 0] - x0);//a26
                double cx, cy;//像点坐标近似值
                cx = x0 - f * (R[0, 0] * (Po[i, 0] - Xs) + R[0, 1] * (Po[i, 1] - Ys) + R[0, 2] * (Po[i, 2] - Zs)) / (R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs));
                cy = y0 - f * (R[1, 0] * (Po[i, 0] - Xs) + R[1, 1] * (Po[i, 1] - Ys) + R[1, 2] * (Po[i, 2] - Zs)) / (R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs));
                l[2 * i, 0] = c[i * 2, 0] - cx;
                l[2 * i + 1, 0] = c[i * 2 + 1, 0] - cy;
            }
            dX = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), l); //第一组结果
            Xs += dX[0, 0]; Ys += dX[1, 0]; Zs += dX[2, 0]; fi += dX[3, 0]; w += dX[4, 0]; k += dX[5, 0];
            writer.WriteLine("改正值 {0}  {1}  {2}  {3}  {4}  {5}", dX[0, 0], dX[1, 0], dX[2, 0], dX[3, 0], dX[4, 0], dX[5, 0]);
            // V = Matrix.SubMatrix(Matrix.MultiplyMatrix(A, dX), l);
            // c1 = Matrix.AddMatrix(c1, V);
            double judge = 0.1 / 60 * 2 * Math.PI / 360;
            while (Math.Abs(dX[0, 0]) > 0.1 || Math.Abs(dX[1, 0]) > 0.1 || Math.Abs(dX[2, 0]) > 0.1 || Math.Abs(dX[3, 0]) > judge || Math.Abs(dX[4, 0]) > judge || Math.Abs(dX[5, 0]) > judge)
            {
                V = Matrix.SubMatrix(Matrix.MultiplyMatrix(A, dX), l);
                c = Matrix.AddMatrix(c, V);//更正像点
                writer.WriteLine("{0}  {1}  {2}  {3}  {4}  {5}", Xs, Ys, Zs, fi, w, k);
                for (int i = 0; i < n / 2; i++)
                {

                    double[,] R = IR(fi, w, k);
                    //r = Math.Sqrt(Math.Pow(c1[i*2,0] - x0, 2) + Math.Pow(c1[i*2+1,0] - y0, 2));
                    //double dx = (c1[i*2,0] - x0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(c1[i*2,0] - x0, 2)) + 2 * p2 * (c1[i*2,0] - x0) * (c1[i*2+1,0] - y0) + a * (c1[i*2,0] - x0) + B * (c1[i*2+1,0] - y0);
                    //double dy = (c1[i*2+1,0] - y0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(c1[i*2+1,0] - y0, 2) + 2 * p1 * (c1[i*2,0] - x0) * (c1[i*2+1,0] - y0));
                    //double tx0 =CX[i*2], ty0 = CX[i*2+1];//像主点坐标 
                    //c1[i * 2, 0] -= dx;c1[i * 2 + 1, 0] -= dy;             
                    double Z = R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs);
                    A[i * 2, 0] = (R[0, 0] * f + R[2, 0] * (c[i * 2, 0] - x0)) / Z;//a11
                    A[i * 2, 1] = (R[0, 1] * f + R[2, 1] * (c[i * 2, 0] - x0)) / Z;//a12
                    A[i * 2, 2] = (R[0, 2] * f + R[2, 2] * (c[i * 2, 0] - x0)) / Z;//a13
                    A[i * 2, 3] = (c[i * 2 + 1, 0] - y0) * Math.Sin(w) - ((c[i * 2, 0] - x0) / f * ((c[i * 2, 0] - x0) * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) * Math.Sin(k)) + f * Math.Cos(k)) * Math.Cos(w);//a14
                    A[i * 2, 4] = -f * Math.Sin(k) - (c[i * 2, 0] - x0) / f * ((c[i * 2, 0] - x0) * Math.Sin(k) + (c[i * 2 + 1, 0] - y0) * Math.Cos(k));//a15
                    A[i * 2, 5] = (c[i * 2 + 1, 0] - y0);//a16

                    A[i * 2 + 1, 0] = (R[1, 0] * f + R[2, 0] * (c[i * 2 + 1, 0] - y0)) / Z;//a21
                    A[i * 2 + 1, 1] = (R[1, 1] * f + R[2, 1] * (c[i * 2 + 1, 0] - y0)) / Z;//a22
                    A[i * 2 + 1, 2] = (R[1, 2] * f + R[2, 2] * (c[i * 2 + 1, 0] - y0)) / Z;//a23
                    A[i * 2 + 1, 3] = -(c[i * 2, 0] - x0) * Math.Sin(w) - ((c[i * 2 + 1, 0] - y0) / f * ((c[i * 2, 0] - x0) * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) * Math.Sin(k)) - f * Math.Sin(k)) * Math.Cos(w);//a24
                    A[i * 2 + 1, 4] = -f * Math.Cos(k) - (c[i * 2 + 1, 0] - y0) / f * ((c[i * 2, 0] - x0) * Math.Sin(k) + (c[i * 2 + 1, 0] - y0) * Math.Cos(k));//a25
                    A[i * 2 + 1, 5] = -(c[i * 2, 0] - x0);//a26
                    double cx, cy;//像点坐标近似值
                    cx = x0 - f * (R[0, 0] * (Po[i, 0] - Xs) + R[0, 1] * (Po[i, 1] - Ys) + R[0, 2] * (Po[i, 2] - Zs)) / (R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs));
                    cy = y0 - f * (R[1, 0] * (Po[i, 0] - Xs) + R[1, 1] * (Po[i, 1] - Ys) + R[1, 2] * (Po[i, 2] - Zs)) / (R[2, 0] * (Po[i, 0] - Xs) + R[2, 1] * (Po[i, 1] - Ys) + R[2, 2] * (Po[i, 2] - Zs));
                    l[2 * i, 0] = c[i * 2, 0] - cx;
                    l[2 * i + 1, 0] = c[i * 2 + 1, 0] - cy;
                }
                dX = Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), l); //第一组结果
                Xs += dX[0, 0]; Ys += dX[1, 0]; Zs += dX[2, 0]; fi += dX[3, 0]; w += dX[4, 0]; k += dX[5, 0];//更新外方位元素
                writer.WriteLine("改正值 {0}  {1}  {2}  {3}  {4}  {5}", dX[0, 0], dX[1, 0], dX[2, 0], dX[3, 0], dX[4, 0], dX[5, 0]);
            }
            writer.WriteLine("{0}  {1}  {2}  {3}  {4}  {5}", Xs, Ys, Zs, fi, w, k);
            double[,] tmp = Matrix.MultiplyMatrix(Matrix.Transpose(V), V);
            double m0 = Math.Sqrt(tmp[0,0] / (2 * n / 2 - 6));
            double[,] Q = Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A));
           
            double[,] m = new double[6, 1];
            m[0, 0] = Math.Sqrt(Q[0, 0]) * m0;m[1, 0] = Math.Sqrt(Q[1, 1]) * m0; m[2, 0] = Math.Sqrt(Q[2, 2]) * m0;
            m[3, 0] = Math.Sqrt(Q[3, 3]) * m0;m[4, 0] = Math.Sqrt(Q[4, 4]) * m0;m[5, 0] = Math.Sqrt(Q[5, 5]) * m0;
            return m;
        }
        public static double[,] IR(double fi, double w, double k)
        {
            double[,] R = new double[3, 3];
            R[0, 0] = Math.Cos(fi) * Math.Cos(k) - Math.Sin(fi) * Math.Sin(w) * Math.Sin(k); //a1
            R[1, 0] = -Math.Cos(fi) * Math.Sin(k) - Math.Sin(fi) * Math.Sin(w) * Math.Cos(k);//a2
            R[2, 0] = -Math.Sin(fi) * Math.Cos(w);//a3
            R[0, 1] = Math.Cos(w) * Math.Sin(k);//b1
            R[1, 1] = Math.Cos(w) * Math.Cos(k);//b2
            R[2, 1] = -Math.Sin(w);//b3
            R[0, 2] = Math.Sin(fi) * Math.Cos(k) + Math.Cos(fi) * Math.Sin(w) * Math.Sin(k);//c1
            R[1, 2] = -Math.Sin(fi) * Math.Sin(k) + Math.Cos(fi) * Math.Sin(w) * Math.Cos(k);//c2
            R[2, 2] = Math.Cos(fi) * Math.Cos(w);//c3

            return R;
        }
       
        static void Main(string[] args)
        {
            string filename = "bundleadjustment_SCBA_Camera_Result.scbacmr";
            StreamReader reader = new StreamReader(filename);
             double x0,y0,f,m,k1, k2, p1, p2, a, B;
            reader.ReadLine();//过掉第一行
            string line = reader.ReadLine();
            reader.Close();
           /* line = line.Replace("\t"," ");
            string tmp = line.Replace("  "," ");
            while (tmp != line)
            {
                line = tmp;
                tmp = line.Replace("  "," ");
            }*/
            string[] tmp = line.Split('\t');
            //int j = 0;
            x0 = double.Parse(tmp[0]);
            y0 = double.Parse(tmp[1]);
            f = double.Parse(tmp[2]);
            //f *= 0.001;
            m = double.Parse(tmp[5]);
            k1 = double.Parse(tmp[6]);
            k2 = double.Parse(tmp[7]);
            p1 = double.Parse(tmp[8]);
            p2 = double.Parse(tmp[9]);
            a = double.Parse(tmp[10]);
            B = double.Parse(tmp[11]);

            Array.Clear(tmp, 0, tmp.Length);
            double fi, w, k;
            double Xs, Ys, Zs; //第一个相机
            double fib, wb, kb;
            double Xsb, Ysb, Zsb;//第二个相机
            
           // Zs = f/m;//初始Zs坐标
            //fi = w = k = 0; //三个外方位元素初始化

            int n;
            filename = "bundleadjustment_SCBA_Point_Result.scbapts"; //读取点文件
            reader = new StreamReader(filename);
            line = reader.ReadLine();
            tmp = line.Split('\t');
            n = int.Parse(tmp[0]);
           // double[,] A = new double[n / 2 * 2, 6];
            //double[,] l = new double[n / 2 * 2,1];
            //double[,] V = new double[n / 2 * 2,1];
            //double[] P1 = new double[6];//第一个相机的6个外方位元素
            //double[] P2 = new double[6];//第二个相机的6个外方位元素
            double[,] c1 = new double[n / 2*2,1];//第一个相机的像点坐标
            double[,] c2 = new double[n / 2*2,1];//第二个相机的像点坐标
            double[,] Po = new double[n / 2, 3];//每个点的物方坐标
           // double[] CX = new double[n / 2 * 2];//对于每个点更正畸变后的像原点坐标
          
            //double[,] R = new double[3, 3];
            //double[,] dX = new double[6,1];//6个外方位元素改正值
           // double SX=0, SY=0;//横纵坐标之和
           // double rx, ry;//像点坐标的近似值
            reader.ReadLine();//跳过一行                 756.0875	-131.6018	-15.0643	0.225967	0.112503	-0.081294
           // Xs = 756.0875; Ys = -131.6018; Zs = -15.0643;
            // fi = 0.225967; w = 0.112503; k = -0.081294;
            //fi = w = k = 0;
            for (int i = 0; i < n / 2; i++)
            {
                line = reader.ReadLine();
                tmp = line.Split('\t');
                double X, Y, Z;
                X = double.Parse(tmp[1]);
                Y = double.Parse(tmp[2]);
                Z = double.Parse(tmp[3]);
                Po[i, 0] = X;
                Po[i, 1] = Y;
                Po[i, 2] = Z;
                //SX += X;SY += Y;
                reader.ReadLine();//跳过一行
                Array.Clear(tmp,0,tmp.Length);
                line = reader.ReadLine();
                tmp = line.Split('\t');
                c1[i*2,0] = double.Parse(tmp[2]);
                c1[i*2+1,0] = double.Parse(tmp[3]);
                line = reader.ReadLine();
                tmp = line.Split('\t');
                c2[i*2,0] = double.Parse(tmp[2]);
                c2[i*2+1,0] = double.Parse(tmp[3]);
            }
            reader.Close();
            //读入外方位元素初值
            filename = "bundleadjustment_SCBA_Photo_Result.scbapht";
            reader = new StreamReader(filename);
            reader.ReadLine();
            line = reader.ReadLine();
            tmp = line.Split('\t');
            Xs = double.Parse(tmp[0]);Ys = double.Parse(tmp[1]);Zs = double.Parse(tmp[2]);
            fi = double.Parse(tmp[3]);w = double.Parse(tmp[4]);k = double.Parse(tmp[5]);
            StreamWriter Writer = new StreamWriter("result and precision");
            StreamWriter writer1 = new StreamWriter("第一张相片");
            double[,]M=cal(x0,y0,f,k1,k2,p1,p2,a,B,ref Xs,ref Ys,ref Zs,ref fi,ref w,ref k,n,c1,Po,writer1);
            writer1.Close();
            Writer.WriteLine("Xs={0},Ys={1},Zs={2},fi={3},w={4},k={5}", Xs, Ys, Zs, fi, w, k);
            Writer.WriteLine("mxs={0},mys={1},mzs={2},mfi={3},mw={4},mk={5}", M[0,0],M[1,0], M[2,0], M[3,0], M[4,0], M[5,0]);
            //Array.Clear(tmp, 0, tmp.Length);
            line = reader.ReadLine();
            reader.Close();
            tmp = line.Split('\t');
            Xsb = double.Parse(tmp[0]); Ysb = double.Parse(tmp[1]); Zsb = double.Parse(tmp[2]);
            fib = double.Parse(tmp[4]); wb = double.Parse(tmp[5]); kb = double.Parse(tmp[6]);
            writer1 = new StreamWriter("第二张相片");
            M =cal(x0, y0, f, k1, k2, p1, p2, a, B,ref Xsb,ref Ysb,ref Zsb,ref fib,ref wb,ref kb, n, c2, Po,writer1);
            writer1.Close();            
            Writer.WriteLine("Xsb={0},Ysb={1},Zsb={2},fib={3},wb={4},kb={5}", Xsb, Ysb, Zsb, fib, wb, kb);
            Writer.WriteLine("mxsb={0},mysb={1},mzsb={2},mfib={3},mwb={4},mkb={5}", M[0, 0], M[1, 0], M[2, 0], M[3, 0], M[4, 0], M[5, 0]);
            Writer.Close();
            //前方交会
            double[,] R1 = new double[3, 3];//第一个相机的旋转矩阵
            double[,] R2 = new double[3, 3];//第二个相机的旋转矩阵
            double[,] A = new double[4, 3];
            double[,] l = new double[4, 1];
            double[,] x = new double[3,1];
            R1 = IR(fi, w, k);
            R2 = IR(fib, wb, kb);
            StreamWriter writer = new StreamWriter("data.out");

            filename = "bundleadjustment_SCBA_Point_Result.scbapts";
            reader = new StreamReader(filename);
            reader.ReadLine();
            reader.ReadLine();
            for (int i = 0; i < n / 2; i++)
            {
                reader.ReadLine(); 
                reader.ReadLine();//跳过一行
                reader.ReadLine();
                reader.ReadLine();
                
            }
            for (int i =0; i < n-n/2-1; i++)
            {
                line=reader.ReadLine();//点号
                tmp = line.Split('\t');
                string Pname = tmp[0];
                reader.ReadLine();//跳过一行
                line = reader.ReadLine();
                tmp = line.Split('\t');
                c1[i * 2, 0] = double.Parse(tmp[2]);
                c1[i * 2 + 1, 0] = double.Parse(tmp[3]);
                double r = Math.Sqrt(Math.Pow(c1[i * 2, 0] - x0, 2) + Math.Pow(c1[i * 2 + 1, 0] - y0, 2));
                double dx = (c1[i * 2, 0] - x0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(c1[i * 2, 0] - x0, 2)) + 2 * p2 * (c1[i * 2, 0] - x0) * (c1[i * 2 + 1, 0] - y0) + a * (c1[i * 2, 0] - x0) + B * (c1[i * 2 + 1, 0] - y0);
                double dy = (c1[i * 2 + 1, 0] - y0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(c1[i * 2 + 1, 0] - y0, 2) + 2 * p1 * (c1[i * 2, 0] - x0) * (c1[i * 2 + 1, 0] - y0));
                c1[i * 2, 0] -= dx; c1[i * 2 + 1, 0] -= dy;
                line = reader.ReadLine();
                tmp = line.Split('\t');
                c2[i * 2, 0] = double.Parse(tmp[2]);
                c2[i * 2 + 1, 0] = double.Parse(tmp[3]); 
                r = Math.Sqrt(Math.Pow(c2[i * 2, 0] - x0, 2) + Math.Pow(c2[i * 2 + 1, 0] - y0, 2));
                dx = (c2[i * 2, 0] - x0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p1 * (r * r + 2 * Math.Pow(c2[i * 2, 0] - x0, 2)) + 2 * p2 * (c2[i * 2, 0] - x0) * (c2[i * 2 + 1, 0] - y0) + a * (c2[i * 2, 0] - x0) + B * (c2[i * 2 + 1, 0] - y0);
                dy = (c2[i * 2 + 1, 0] - y0) * (k1 * r * r + k2 * Math.Pow(r, 4)) + p2 * (r * r + 2 * Math.Pow(c2[i * 2 + 1, 0] - y0, 2) + 2 * p1 * (c2[i * 2, 0] - x0) * (c2[i * 2 + 1, 0] - y0));
                c2[i * 2, 0] -= dx; c2[i * 2 + 1, 0] -= dy;
                A[0, 0] = f * R1[0, 0] + (c1[i * 2, 0] - x0) * R1[2, 0];
                A[0, 1] = f * R1[0, 1] + (c1[i * 2, 0] - x0) * R1[2, 1];
                A[0, 2] = f * R1[0, 2] + (c1[i * 2, 0] - x0) * R1[2, 2];
                A[1, 0] = f * R1[1, 0] + (c1[i * 2 + 1, 0] - y0) * R1[2, 0];
                A[1, 1] = f * R1[1, 1] + (c1[i * 2 + 1, 0] - y0) * R1[2, 1];
                A[1, 2] = f * R1[1, 2] + (c1[i * 2 + 1, 0] - y0) * R1[2, 2];
                l[0, 0] = f * R1[0, 0] * Xs + f * R1[0, 1] * Ys + f * R1[0, 2] * Zs + (c1[i * 2, 0] - x0) * R1[2, 0] * Xs + (c1[i * 2, 0] - x0) * R1[2, 1] * Ys + (c1[i * 2, 0] - x0) * R1[2, 2] * Zs;
                l[1, 0] = f * R1[1, 0] * Xs + f * R1[1, 1] * Ys + f * R1[1, 2] * Zs + (c1[i * 2+1, 0] - y0) * R1[2, 0] * Xs + (c1[i * 2+1, 0] - y0) * R1[2, 1] * Ys + (c1[i * 2+1, 0] - y0) * R1[2, 2] * Zs;

                A[2, 0] = f * R2[0, 0] + (c2[i * 2, 0] - x0) * R2[2, 0];
                A[2, 1] = f * R2[0, 1] + (c2[i * 2, 0] - x0) * R2[2, 1];
                A[2, 2] = f * R2[0, 2] + (c2[i * 2, 0] - x0) * R2[2, 2];
                A[3, 0] = f * R2[1, 0] + (c2[i * 2 + 1, 0] - y0) * R2[2, 0];
                A[3, 1] = f * R2[1, 1] + (c2[i * 2 + 1, 0] - y0) * R2[2, 1];
                A[3, 2] = f * R2[1, 2] + (c2[i * 2 + 1, 0] - y0) * R2[2, 2];
                l[2, 0] = f * R2[0, 0] * Xsb + f * R2[0, 1] * Ysb + f * R2[0, 2] * Zsb + (c2[i * 2, 0] - x0) * R2[2, 0] * Xsb + (c2[i * 2, 0] - x0) * R2[2, 1] * Ysb + (c2[i * 2, 0] - x0) * R2[2, 2] * Zsb;
                l[3, 0] = f * R2[1, 0] * Xsb + f * R2[1, 1] * Ysb + f * R2[1, 2] * Zsb + (c2[i * 2 + 1, 0] - y0) * R2[2, 0] * Xsb + (c2[i * 2 + 1, 0] - y0) * R2[2, 1] * Ysb + (c2[i * 2 + 1, 0] - y0) * R2[2, 2] * Zsb;
                x= Matrix.MultiplyMatrix(Matrix.MultiplyMatrix(Matrix.Athwart(Matrix.MultiplyMatrix(Matrix.Transpose(A), A)), Matrix.Transpose(A)), l);
                writer.Write("点号：{0} 坐标{1} {2} {3}",Pname,x[0,0],x[1,0],x[2,0]);
                writer.WriteLine();
            }//读入另外一半点的像点坐标
            reader.Close();
            writer.Close();

            Console.ReadLine();

        }
    }
}
