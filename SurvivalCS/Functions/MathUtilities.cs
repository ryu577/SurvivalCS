using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SurvivalCS.Functions
{
    class MathUtilities
    {
        /// <summary>
        /// Numerically integrates a function and returns the result over
        /// an entire range of values. To be used in place of MathNet functions when 
        /// a whole range of values is needed.
        /// </summary>
        /// <param name="func">The function to be integrated.</param>
        /// <param name="upperBound">The upper bound over which the function is to be integrated.</param>
        /// <param name="delta">The delta. Should be small, preferably don't change from default.</param>
        /// <returns>A double array with the values of the numerically integrated function
        /// between zero and upperBound at intervals of delta</returns>
        public static double[] CumulativeIntegral(Func<double, double> func, double upperBound, double delta = 1e-4)
        {
            int numTerms = (int)(upperBound / delta);
            double[] result = new double[numTerms];

            result[0] = delta * func(delta / 2);
            for (int i = 1; i < numTerms; i++)
            {
                result[i] = result[i - 1] + (func(i + (delta / 2)) * delta);
            }

            return result;
        }

        /// <summary>
        /// Compute Log(1 + x) without losing precision for small values of x.
        /// </summary>
        /// <param name="x">Variable for the function.</param>
        /// <returns>Log(1 + x)</returns>
        public static double Log1P(double x)
        {
            if (Math.Abs(x) > 1e-4)
            {
                return Math.Log(1.0 + x);
            }

            return (1.0 + (-0.5 * x)) * x;
        }

        /// <summary>
        /// Computes the Pearson coefficient of correlation between two vectors.
        /// Taken from http://stackoverflow.com/questions/35702231/how-to-compute-pearson-correlation-between-2-given-vectors
        /// </summary>
        /// <param name="xs">The first list of doubles</param>
        /// <param name="ys">The second list of doubles</param>
        /// <returns>The Person coefficient of correlation as a double</returns>
        public static double Correlation(IEnumerable<double> xs, IEnumerable<double> ys)
        {
            // sums of x, y, x squared etc.
            double sx = 0.0;
            double sy = 0.0;
            double sxx = 0.0;
            double syy = 0.0;
            double sxy = 0.0;

            int n = 0;

            using (var enX = xs.GetEnumerator())
            {
                using (var enY = ys.GetEnumerator())
                {
                    while (enX.MoveNext() && enY.MoveNext())
                    {
                        double x = enX.Current;
                        double y = enY.Current;

                        n += 1;
                        sx += x;
                        sy += y;
                        sxx += x * x;
                        syy += y * y;
                        sxy += x * y;
                    }
                }
            }

            // covariation
            double cov = (sxy / n) - (((sx * sy) / n) / n);

            // standard error of x
            double sigmaX = Math.Sqrt((sxx / n) - (((sx * sx) / n) / n));

            // standard error of y
            double sigmaY = Math.Sqrt((syy / n) - (((sy * sy) / n) / n));

            // correlation is just a normalized covariation
            return (cov / sigmaX) / sigmaY;
        }
    }
}
