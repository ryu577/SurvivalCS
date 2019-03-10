using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics.CodeAnalysis;
using MathNet.Numerics.LinearAlgebra;


namespace SurvivalCS.Distributions
{
    public class LogLogistic : BaseModel
    {
        /// <summary>
        /// The shape parameter of the loglogistic
        /// </summary>
        private double beta;

        /// <summary>
        /// The scale parameter of the loglogistic
        /// </summary>
        private double alpha;

        public double Mean
        {
            get
            {
                return this.alpha / this.beta * Math.PI / Math.Sin(Math.PI/this.beta);
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogLogistic"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="interventionCount">Number of interventions.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public LogLogistic(
            List<double> organicRecoveryDurations,
            double interventionCount,
            double interventionThreshold) : base(organicRecoveryDurations, interventionCount, interventionThreshold)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogLogistic"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">A list of durations of transitions from unhealthy to ready.</param>
        /// <param name="inorganicRecoveryDurations">A list of durations for unhealthy to powering on.</param>
        public LogLogistic(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations) : 
            base(organicRecoveryDurations, inorganicRecoveryDurations)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogLogistic"/> class.
        /// </summary>
        /// <param name="beta">The shape parameter of the log logistic</param>
        /// <param name="alpha">The scale parameter of the log logistic</param>
        public LogLogistic(double beta, double alpha) : base(beta, alpha)
        {
            this.alpha = alpha;
            this.beta = beta;
            this.ShapeUpperBound = Math.Max(beta * 4, 10);
            this.ScaleUpperBound = Math.Max(alpha * 4, 1000);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LogLogistic"/> class with features.
        /// </summary>
        /// <param name="organicRecoveryDurations">Organic recover time data.</param>
        /// <param name="inorganicRecoveryDurations">Inorganic recover time data.</param>
        /// <param name="organicRecoveryFeatureData">Organic recover feature data.</param>
        /// <param name="inorganicRecoveryFeatureData">Inorganic recover feature data.</param>
        public LogLogistic(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations,
            Matrix<double> organicRecoveryFeatureData,
            Matrix<double> inorganicRecoveryFeatureData) : 
        base(organicRecoveryDurations, inorganicRecoveryDurations, organicRecoveryFeatureData, 
                inorganicRecoveryFeatureData)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Calculates the PDF of the <see cref="LogLogistic"/>.
        /// </summary>
        /// <param name="t">The value at which the PDF is to be calculated.</param>
        /// <param name="beta">The shape parameter of the log logistic.</param>
        /// <param name="alpha">The scale parameter of the log logistic.</param>
        /// <returns>A double which is the PDF of the log logistic at t.</returns>
        public double PDF(double t, double beta, double alpha)
        {
            return (beta / alpha) * Math.Pow(t / alpha, beta - 1) / Math.Pow(1 + Math.Pow(t / alpha, beta), 2);
        }

        /// <summary>
        /// Calculates the PDF of the log logostic distribution given the specified value. 
        /// A wrapper for the version where the parameters are also needed as input.
        /// </summary>
        /// <param name="t">The value at which the PDF is to be calculated.</param>
        /// <returns>A double which is the value of the PDF.</returns>
        public override double PDF(double t)
        {
            return this.PDF(t, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the log of the PDF of the <see cref="LogLogistic"/>.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="beta">The shape parameter of the log logistic</param>
        /// <param name="alpha">The scale parameter of the log logistic</param>
        /// <returns>A double which is the log pdf.</returns>
        public override double LogPdf(double t, double beta, double alpha)
        {
            return Math.Log(beta) - Math.Log(alpha) + ((beta - 1) * (Math.Log(t) - Math.Log(alpha)))
                - (2 * Math.Log(1 + Math.Pow(t / alpha, beta)));
        }

        /// <summary>
        /// A wrapper for the log pdf that uses the instancs shape and scale parameters to get the log pdf.
        /// </summary>
        /// <param name="t">The value at which the log of the pdf is to be calculated.</param>
        /// <returns>The logarithm of the pdf function.</returns>
        public double LogPdf(double t)
        {
            return this.LogPdf(t, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the cumulative density function of the <see cref="LogLogistic"/> distribution.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="beta">The shape parameter of the log logistic.</param>
        /// <param name="alpha">The scale parameter of the log logistic.</param>
        /// <returns>A double which is the CDF.</returns>
        public double CDF(double t, double beta, double alpha)
        {
            return 1 / (1 + Math.Pow(t / alpha, -beta));
        }

        /// <summary>
        /// Calculates the inverse CDF for the Log-logistic distribution.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <param name="beta">The shape parameter for the Log-logistic distribution.</param>
        /// <param name="alpha">The scale parameter for the Log-logistic distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public double InverseCDF(double p, double beta, double alpha)
        {
            return alpha * Math.Pow(1 / p - 1, -1 / beta);
        }

        /// <summary>
        /// A wrapper for the inverse CDF for the Log-logistic distribution that uses the shape and scale for the instance.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public override double InverseCDF(double p)
        {
            return this.InverseCDF(p, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the survival function, which is the area under the PDF from t to infinity.
        /// </summary>
        /// <param name="t">The value at which the survival is to be calculated.</param>
        /// <param name="beta">The shape parameter of the log logistic.</param>
        /// <param name="alpha">The scale parameter of the log logistic.</param>
        /// <returns>A double which is the survival function.</returns>
        public double Survival(double t, double beta, double alpha)
        {
            return 1 - this.CDF(t, beta, alpha);
        }

        /// <summary>
        /// A wrapper for the survival function that takes parameters from current instance.
        /// </summary>
        /// <param name="x">The value at which the survival is to be calculated.</param>
        /// <returns>A double which is the survival function.</returns>
        public override double Survival(double x)
        {
            return this.Survival(x, this.beta, this.alpha);
        }

        /// <summary>
        /// Calculates the log of the survival function. TODO: calculate this explicitly to make it more efficient.
        /// </summary>
        /// <param name="t">The value at which the log survival is to be calculated.</param>
        /// <param name="beta">The shape parameter of the log logistic.</param>
        /// <param name="alpha">The scale parameter of the log logistic.</param>
        /// <returns>The log of the survival function at a particular value.</returns>
        public override double LogSurvival(double t, double beta, double alpha)
        {
            // TODO: Calculate the logarithm of the function explicitly for better speed.
            return Math.Log(this.Survival(t, beta, alpha));
        }

        /// <summary>
        /// A wrapper for the logsurvival function that uses the alpha and beta parameters from the instance.
        /// </summary>
        /// <param name="t">The value at which the log survival is to be calculated.</param>
        /// <returns>The logarithm of the survival function.</returns>
        public double LogSurvival(double t)
        {
            return this.LogSurvival(t, this.beta, this.alpha);
        }

        
        /// <summary>
        /// Calculates the log of the likelihood function given the censored and un-censored data.
        /// </summary>
        /// <returns>A double which is the log of the likelihood.</returns>
        public double LogLikelihood()
        {
            // double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            return t.Sum(ti => this.LogPdf(ti)) + x.Sum(xi => this.LogSurvival(xi));
        }

        /// <summary>
        /// Calculates the gradient of the Log-likelihood function for the Loglogistic without features.
        /// </summary>
        /// <param name="beta">The shape parameter of the Loglogistic.</param>
        /// <param name="alpha">The scale parameter of the Loglogisti.c</param>
        /// <returns>A double array - the gradients for beta and alpha.</returns>
        public override Vector<double> GradLL(double beta, double alpha)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;

            double delalpha = -((n * beta) / alpha) + (((2 * beta) / Math.Pow(alpha, beta + 1)) *
                t.Sum(ti => Math.Pow(ti, beta) / (1 + Math.Pow(ti / alpha, beta)))) 
                + ((beta / Math.Pow(alpha, beta + 1))
                * x.Sum(xi => Math.Pow(xi, beta) / (1 + Math.Pow(xi / alpha, beta))));

            double delbeta = (n / beta) - (n * Math.Log(alpha)) + t.Sum(ti => Math.Log(ti))
                - (2 * t.Sum(ti => Math.Pow(ti / alpha, beta) / (1 + Math.Pow(ti / alpha, beta)) 
                    * Math.Log(ti / alpha)))
                - x.Sum(xi => (Math.Pow(xi / alpha, beta) / (1 + Math.Pow(xi / alpha, beta))) 
                    * Math.Log(xi / alpha));

            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delbeta, delalpha });
            return res;
        }

        /// <summary>
        /// Calculate the gradient of the log of PDF at a certain point.
        /// </summary>
        /// <param name="t">The point where we calculate the gradient.</param>
        /// <param name="beta">The shape parameter.</param>
        /// <param name="alpha">The scale parameter.</param>
        /// <returns>The gradient of the log-likelihood function at t.</returns>
        public override Vector<double> GradLPDF(double t, double beta, double alpha)
        {
            double temp = Math.Pow(t / alpha, beta);
            double delbeta = 1 / beta + Math.Log(t / alpha) * (1 - 2 * temp / (1 + temp));
            double delalpha = -beta / alpha + 2 * beta * temp / (1 + temp) / alpha;
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delbeta, delalpha });
            return res;
        }

        /// <summary>
        /// Calculate the gradient of log of the survival function at a single point.
        /// </summary>
        /// <param name="x">The point you want to evaluate the gradient.</param>       
        /// <param name="beta">The shape parameter.</param>
        /// <param name="alpha">The scale parameter.</param>
        /// <returns>The gradient of the survival function at x.</returns>
        public override Vector<double> GradLSurvival(double x, double beta, double alpha)
        {
            double temp = Math.Pow(x / alpha, beta);
            double delbeta = -temp / (1 + temp) * Math.Log(x / alpha);
            double delalpha = beta / alpha * temp / (1 + temp);
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delbeta, delalpha });
            return res;
        }
    }
}
