using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SurvivalCS.Functions
{
    /// <summary>
    /// Sigmoid function - maps real numbers to a bounded interval in the form of an S-curve.
    /// </summary>
    public class Sigmoid
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Sigmoid"/> class.
        /// </summary>
        /// <param name="upperBound">The maximum value the transformed result can take.</param>
        public Sigmoid(double upperBound)
        {
            this.UpperBound = upperBound;
        }

        /// <summary>
        /// Gets or sets the upper bound, the maximum value this Sigmoid can take.
        /// </summary>
        public double UpperBound { get; set; }

        /// <summary>
        /// Calculates the inverse of the sigmoid function
        /// </summary>
        /// <param name="x">The result of the sigmoid. Must be between 0 and 1.</param>
        /// <param name="u">The scale of the sigmoid.</param>
        /// <returns>The inverse of the sigmoid.</returns>
        public static double InverseSigmoid(double x, double u = 1)
        {
            return -Math.Log(u / x - 1);
        }

        /// <summary>
        /// Transforms a number from real space to sigmoid space (constrained between 0 and UpperBound).
        /// </summary>
        /// <param name="x">Value at which to evaluate.</param>
        /// <returns>The value of transformed function.</returns>
        public double Transform(double x)
        {
            return this.UpperBound / (1 + Math.Exp(-x));
        }


        public static double Transform(double x, double u)
        {
            return u / (1 + Math.Exp(-x));
        }

        /// <summary>
        /// The gradient of the sigmoid function
        /// </summary>
        /// <param name="x">Input value at which to evaluate</param>
        /// <returns>A double, the derivative</returns>
        public double Grad(double x)
        {
            double p = this.Transform(x) / this.UpperBound;
            return this.UpperBound * p * (1 - p);
        }
    }
}

