using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace matrixcal
{
    public class Matrix
    {

        #region
        ///   <summary> 
        ///   求矩阵的秩，改变了原矩阵 
        ///   </summary> 
        ///   <param   name= "iMatrix "> </param> 
        ///   <returns> </returns> 
        //public   static   int   OrderCount(double[,]   iMatrix) 
        //{ 
        //    //int   i   =   0; 
        //    //for(i   =   0;i <3;i++) 
        //    //{ 
        //    //if(iMatrix[i,i]   !=   0) 
        //    //{ 
        //    //double   intTemp   =   iMatrix[i,i]; 
        //    //for(int   j   =   0;j <3;j++) 
        //    //{ 
        //    //iMatrix[i,j]   =   iMatrix[i,j]/intTemp; 
        //    //} 
        //    //} 
        //    //for(int   j   =   0;j <3;j++) 
        //    //{ 
        //    //if(j   ==   i) 
        //    //continue; 
        //    //double   intTemp   =   iMatrix[j,i]; 
        //    //for(int   k   =   0;k <3;k++) 
        //    //{ 
        //    //iMatrix[j,k]   =   iMatrix[j,k]   -   iMatrix[i,k]   *   intTemp; 
        //    //} 
        //    //} 
        //    //} 
        //    //for(i   =   0;i <3;i++) 
        //    //if(iMatrix[i,i]   ==   0) 
        //    //break; 
        //    //return   i; 
        //} 



        /////   <summary> 
        /////   求矩阵的秩，不改变原矩阵 
        /////   </summary> 
        /////   <param   name= "iMatrix "> </param> 
        /////   <returns> </returns> 
        //public   static   int   OrderCountNotChange(double[,]   iMatrix) 
        //{ 
        //    //double[,]   TempMatrix   =   new   double[3,3]; 
        //    //int   i   =   0; 
        //    //for(i   =   0;i <3;i++) 
        //    //for(int   j   =   0;j <3;j++) 
        //    //TempMatrix[i,j]   =   iMatrix[i,j]; 
        //    //for(i   =   0;i <3;i++) 
        //    //{ 
        //    //if(TempMatrix[i,i]   !=   0) 
        //    //{ 
        //    //double   intTemp   =   TempMatrix[i,i]; 
        //    //for(int   j   =   0;j <3;j++) 
        //    //{ 
        //    //TempMatrix[i,j]   =   TempMatrix[i,j]/intTemp; 
        //    //} 
        //    //} 
        //    //for(int   j   =   0;j <3;j++) 
        //    //{ 
        //    //if(j   ==   i) 
        //    //continue; 
        //    //double   intTemp   =   TempMatrix[j,i]; 
        //    //for(int   k   =   0;k <3;k++) 
        //    //{ 
        //    //TempMatrix[j,k]   =   TempMatrix[j,k]   -   TempMatrix[i,k]   *   intTemp; 
        //    //} 
        //    //} 
        //    //} 
        //    //for(i   =   0;i <3;i++) 
        //    //if(TempMatrix[i,i]   ==   0) 
        //    //break; 
        //    //return   i; 
        //} 
        #endregion

        ///   <summary> 
        ///   矩阵的转置 
        ///   </summary> 
        ///   <param   name= "iMatrix "> </param> 
        public static double[,] Transpose(double[,] iMatrix)
        {
            int row = iMatrix.GetLength(0);
            int column = iMatrix.GetLength(1);
            //double[,] iMatrix = new double[column, row];
            double[,] TempMatrix = new double[row, column];
            double[,] iMatrixT = new double[column, row];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < column; j++)
                {
                    TempMatrix[i, j] = iMatrix[i, j];
                }
            }
            for (int i = 0; i < column; i++)
            {
                for (int j = 0; j < row; j++)
                {
                    iMatrixT[i, j] = TempMatrix[j, i];
                }
            }
            return iMatrixT;

        }

        ///   <summary> 
        ///   矩阵的逆矩阵 
        ///   </summary> 
        ///   <param   name= "iMatrix "> </param> 
        public static double[,] Athwart(double[,] iMatrix)
        {
            int i = 0;
            int row = iMatrix.GetLength(0);
            double[,] MatrixZwei = new double[row, row * 2];
            double[,] iMatrixInv = new double[row, row];
            for (i = 0; i < row; i++)
            {
                for (int j = 0; j < row; j++)
                {
                    MatrixZwei[i, j] = iMatrix[i, j];
                }
            }
            for (i = 0; i < row; i++)
            {
                for (int j = row; j < row * 2; j++)
                {
                    MatrixZwei[i, j] = 0;
                    if (i + row == j)
                        MatrixZwei[i, j] = 1;
                }
            }

            for (i = 0; i < row; i++)
            {
                if (MatrixZwei[i, i] != 0)
                {
                    double intTemp = MatrixZwei[i, i];
                    for (int j = 0; j < row * 2; j++)
                    {
                        MatrixZwei[i, j] = MatrixZwei[i, j] / intTemp;
                    }
                }
                for (int j = 0; j < row; j++)
                {
                    if (j == i)
                        continue;
                    double intTemp = MatrixZwei[j, i];
                    for (int k = 0; k < row * 2; k++)
                    {
                        MatrixZwei[j, k] = MatrixZwei[j, k] - MatrixZwei[i, k] * intTemp;
                    }
                }
            }

            for (i = 0; i < row; i++)
            {
                for (int j = 0; j < row; j++)
                {
                    iMatrixInv[i, j] = MatrixZwei[i, j + row];
                }
            }
            return iMatrixInv;
        }

        ///   <summary> 
        ///   矩阵加法 
        ///   </summary> 
        ///   <param   name= "MatrixEin "> </param> 
        ///   <param   name= "MatrixZwei "> </param> 
        public static double[,] AddMatrix(double[,] MatrixEin, double[,] MatrixZwei)
        {
            double[,] MatrixResult = new double[MatrixEin.GetLength(0), MatrixZwei.GetLength(1)];
            for (int i = 0; i < MatrixEin.GetLength(0); i++)
                for (int j = 0; j < MatrixZwei.GetLength(1); j++)
                    MatrixResult[i, j] = MatrixEin[i, j] + MatrixZwei[i, j];
            return MatrixResult;
        }

        ///   <summary> 
        ///   矩阵减法 
        ///   </summary> 
        ///   <param   name= "MatrixEin "> </param> 
        ///   <param   name= "MatrixZwei "> </param> 
        public static double[,] SubMatrix(double[,] MatrixEin, double[,] MatrixZwei)
        {
            double[,] MatrixResult = new double[MatrixEin.GetLength(0), MatrixZwei.GetLength(1)];
            for (int i = 0; i < MatrixEin.GetLength(0); i++)
                for (int j = 0; j < MatrixZwei.GetLength(1); j++)
                    MatrixResult[i, j] = MatrixEin[i, j] - MatrixZwei[i, j];
            return MatrixResult;
        }

        ///   <summary> 
        ///   矩阵乘法 
        ///   </summary> 
        ///   <param   name= "MatrixEin "> </param> 
        ///   <param   name= "MatrixZwei "> </param> 
        public static double[,] MultiplyMatrix(double[,] MatrixEin, double[,] MatrixZwei)
        {
            double[,] MatrixResult = new double[MatrixEin.GetLength(0), MatrixZwei.GetLength(1)];
            for (int i = 0; i < MatrixEin.GetLength(0); i++)
            {
                for (int j = 0; j < MatrixZwei.GetLength(1); j++)
                {
                    for (int k = 0; k < MatrixEin.GetLength(1); k++)
                    {
                        MatrixResult[i, j] += MatrixEin[i, k] * MatrixZwei[k, j];
                    }
                }
            }
            return MatrixResult;
        }

        ///   <summary> 
        ///   矩阵对应行列式的值 
        ///   </summary> 
        ///   <param   name= "MatrixEin "> </param> 
        ///   <returns> </returns> 
        public static double ResultDeterminant(double[,] MatrixEin)
        {
            return MatrixEin[0, 0] * MatrixEin[1, 1] * MatrixEin[2, 2] + MatrixEin[0, 1] * MatrixEin[1, 2] * MatrixEin[2, 0] + MatrixEin[0, 2] * MatrixEin[1, 0] * MatrixEin[2, 1]
            - MatrixEin[0, 2] * MatrixEin[1, 1] * MatrixEin[2, 0] - MatrixEin[0, 1] * MatrixEin[1, 0] * MatrixEin[2, 2] - MatrixEin[0, 0] * MatrixEin[1, 2] * MatrixEin[2, 1];

        }



    }
}
