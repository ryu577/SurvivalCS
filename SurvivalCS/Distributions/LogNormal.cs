using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics.CodeAnalysis;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;

namespace SurvivalCS.Distributions
{
    public class LogNormal : BaseModel
    {
        /// <summary>
        /// The mean of the normal backing the lognormal
        /// </summary>
        private double mu;

        /// <summary>
        /// The variance of the normal backing the lognormal
        /// </summary>
        private double sigma;

        /// <summary>
        /// Initializes a new instance of the <see cref="LogNormal"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="interventionCount">Number of interventions.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public LogNormal(
            List<double> organicRecoveryDurations,
            double interventionCount,
            double interventionThreshold) : base(organicRecoveryDurations, interventionCount, interventionThreshold)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogNormal"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">A list of durations of transitions from unhealthy to ready.</param>
        /// <param name="inorganicRecoveryDurations">A list of durations for unhealthy to powering on.</param>
        public LogNormal(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations) : 
                base(organicRecoveryDurations, inorganicRecoveryDurations)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogNormal"/> class.
        /// </summary>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        public LogNormal(double mu, double sigma) : base(mu, sigma)
        {
            this.mu = this.Kappa = mu;
            this.sigma = this.Lambda = sigma;
        }

        /// <summary>
        /// Calculates the PDF of the lognormal distribution.
        /// </summary>
        /// <param name="t">The value at which the PDF is to be calculated.</param>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>A double which is the PDF of the log normal.</returns>
        public double PDF(double t, double mu, double sigma)
        {
            return MathNet.Numerics.Distributions.LogNormal.PDF(t, mu, sigma);
        }

        /// <summary>
        /// Calculates the PDF of the lognormal distribution given the specified value. 
        /// A wrapper for the version where the parameters are also needed as input.
        /// </summary>
        /// <param name="t">The value at which the PDF is to be calculated.</param>
        /// <returns>A double which is the value of the PDF.</returns>
        public override double PDF(double t)
        {
            return this.PDF(t, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the log of the PDF of the lognormal.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>A double which is the log pdf.</returns>
        public override double LogPdf(double t, double mu, double sigma)
        {
            return (-Math.Log(2 * Math.PI) / 2) - Math.Log(sigma) - Math.Log(t)
                - ((Math.Log(t) - mu) * (Math.Log(t) - mu) / (2 * sigma * sigma));
        }

        /// <summary>
        /// Calculates the cumulative density function of the log normal distribution.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>A double which is the CDF.</returns>
        public double CDF(double t, double mu, double sigma)
        {
            var t1 = Math.Log(t);
            return MathNet.Numerics.Distributions.Normal.CDF(mu, sigma, t1);
        }

        /// <summary>
        /// A wrapper for the inverse CDF for the distribution that uses the shape and scale for the instance.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public override double InverseCDF(double p)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the survival function, which is the area under the PDF from some value to infinity.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>A double which is the survival function.</returns>
        public double Survival(double t, double mu, double sigma)
        {
            return 1 - this.CDF(t, mu, sigma);
        }

        /// <summary>
        /// A wrapper for the survival function that feeds in parameters automatically.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <returns>A double which is the survival function.</returns>
        public override double Survival(double t)
        {
            return this.Survival(t, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the log of the survival function. TODO: calculate this explicitly to make it more efficient.
        /// </summary>
        /// <param name="t">The value at which the log survival is to be calculated.</param>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>The log of the survival function at a particular value.</returns>
        public override double LogSurvival(double t, double mu, double sigma)
        {
            return Math.Log(this.Survival(t, mu, sigma));
        }

        /// <summary>
        /// Calculate the gradient of the log of PDF at a given point.
        /// </summary>
        /// <param name="t">The point to evaluate the gradient.</param>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The variance.</param>
        /// <returns>The gradient of log PDF at t.</returns>
        public override Vector<double> GradLPDF(double t, double mu, double sigma)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculate the gradient of log of the survival function at a given point.
        /// </summary>
        /// <param name="x">The point to evaluate the gradient.</param>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The variance.</param>
        /// <returns>The gradient of the log survival function at x.</returns>
        public override Vector<double> GradLSurvival(double x, double mu, double sigma)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the gradient of the Log-likelihood function for the Lognormal without features.
        /// </summary>
        /// <param name="mu">The mean of the normal backing the lognormal</param>
        /// <param name="sigma">The variance of the normal backing the lognormal</param>
        /// <returns>A double array - the gradients for mu and sigma</returns>
        public override Vector<double> GradLL(double mu, double sigma)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            var z = x.Select(xi => (Math.Log(xi) - mu) / sigma);
            double delmu = t.Sum(ti => (Math.Log(ti) - mu)) / (sigma * sigma) + z.Sum(zi => Normal.PDF(0, 1, zi) 
                / Normal.CDF(0, 1, -zi) / sigma);

            double delsigma1 = -(n / sigma);
            double delsigma2 = t.Sum(ti => (Math.Log(ti) - mu) * (Math.Log(ti) - mu)) / (sigma * sigma * sigma);
            double delsigma3 = z.Sum(zi => zi * Normal.PDF(0, 1, zi) / Normal.CDF(0, 1, -zi)) / sigma;

            double delsigma = delsigma1 + delsigma2 + delsigma3;

            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delmu, delsigma });
            return res;
        }
    }
}
