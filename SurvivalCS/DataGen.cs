using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SurvivalCS.Distributions;
using MathNet.Numerics.LinearAlgebra;

namespace SurvivalCS
{
    public class DataGen
    {
        public List<double> organicRecoveryDurations;
        public List<double> inorganicRecoverydurations;

        public Matrix<double> fSamples;
        public Matrix<double> fCensored;

        public void GenLogLogisticWFeatures()
        {
            LogLogistic ll = new LogLogistic(1.2, 300.0);
            List<double> llSamples = ll.GenerateSample(2000);

            var censorLvl = llSamples.Sum()/llSamples.Count();
            var t = CensorList(llSamples, censorLvl);

            int n_x = llSamples.Count(ti => ti > censorLvl);
            
            var x = NPOnes(n_x, censorLvl);

            List<double[]> fSamplesData = new List<double[]>();
            List<double[]> fCensoredData = new List<double[]>();

            for (int i = 0; i < t.Count(); i++)
            {
                fSamplesData.Add(new double[] { 1.0, 2.0, 3.0 });
            }

            for (int i = 0; i < x.Count(); i++)
            {
                fCensoredData.Add(new double[] { 1.0, 2.0, 3.0 });
            }

            ll = new LogLogistic(0.7, 80.0);
            llSamples = ll.GenerateSample(2000);

            censorLvl = llSamples.Sum() / llSamples.Count();
            var t1 = CensorList(llSamples, censorLvl);
            AppendToLst(t, t1);

            n_x = llSamples.Count(ti => ti > censorLvl);
            var x1 = NPOnes(n_x, censorLvl);
            AppendToLst(x, x1);

            for (int i = 0; i < t1.Count(); i++)
            {
                fSamplesData.Add(new double[] { 1.0, 4.0, 2.0 });
            }

            for (int i = 0; i < x1.Count(); i++)
            {
                fCensoredData.Add(new double[] { 1.0, 4.0, 2.0 });
            }

            this.organicRecoveryDurations = t;
            this.inorganicRecoverydurations = x;

            Matrix<double> fSamples = Matrix<double>.Build.DenseOfArray(CreateRectangularArray(fSamplesData));
            Matrix<double> fCensored = Matrix<double>.Build.DenseOfArray(CreateRectangularArray(fCensoredData));
            this.fSamples = fSamples;
            this.fCensored = fCensored;
        }

        public void GenTrivFeaturesData()
        {
            var ti = new List<double> { 1.0, 1.0};
            var xi = new List<double> { 1.0, 1.0};

            List<double[]> fSamplesData = new List<double[]>();
            List<double[]> fCensoredData = new List<double[]>();
            fSamplesData.Add(new double[] { 1.0, 1.0 });
            fSamplesData.Add(new double[] { 1.0, 1.0 });
            fCensoredData.Add(new double[] { 1.0, 1.0 });
            fCensoredData.Add(new double[] { 1.0, 1.0 });
            Matrix<double> fSamples = Matrix<double>.Build.DenseOfArray(CreateRectangularArray(fSamplesData));
            Matrix<double> fCensored = Matrix<double>.Build.DenseOfArray(CreateRectangularArray(fCensoredData));
            this.fSamples = fSamples;
            this.fCensored = fCensored;
            this.organicRecoveryDurations = ti;
            this.inorganicRecoverydurations = xi;
        }

        public static List<double> NPOnes(int n, double val)
        {
            List<double> res = new List<double>();

            for (int i = 0; i < n; i++)
            {
                res.Add(val);
            }

            return res;
        }

        public static void AppendToLst(List<double> left, List<double> right)
        {
            foreach (var val in right)
            {
                left.Add(val);
            }
        }

        public static List<double> CensorList(List<double> lst, double censor)
        {
            List<double> res = new List<double>();

            foreach (var item in lst)
            {
                if (item < censor)
                {
                    res.Add(item);
                }
            }

            return res;
        }

        /// <summary>
        /// Converts list of arrays to 2d array.
        /// Source - http://stackoverflow.com/questions/tagged/c%23
        /// </summary>
        /// <typeparam name="T">Usually a double</typeparam>
        /// <param name="arrays">The list of arrays.</param>
        /// <returns>A two dimensional array with values from the staggered array.</returns>
        public static T[,] CreateRectangularArray<T>(IList<T[]> arrays)
        {
            // TODO: Validation and special-casing for arrays.Count == 0
            int minorLength = arrays[0].Length;
            T[,] ret = new T[arrays.Count, minorLength];
            for (int i = 0; i < arrays.Count; i++)
            {
                var array = arrays[i];
                if (array.Length != minorLength)
                {
                    throw new ArgumentException
                        ("All arrays must be the same length");
                }

                for (int j = 0; j < minorLength; j++)
                {
                    ret[i, j] = array[j];
                }
            }

            return ret;
        }
    }
}
