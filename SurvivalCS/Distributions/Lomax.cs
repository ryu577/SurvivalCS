using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics.CodeAnalysis;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.RootFinding;
using SurvivalCS.Functions;

namespace SurvivalCS.Distributions
{
    public class Lomax : BaseModel
    {
        /// <summary>
        /// Regularizers for fitting.
        /// </summary>
        private double alpha, kappa0, lambda0;

        /// <summary>
        /// Initializes a new instance of the <see cref="Lomax"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="interventionCount">Number of interventions.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public Lomax(
            List<double> organicRecoveryDurations,
            double interventionCount,
            double interventionThreshold) : base(organicRecoveryDurations, interventionCount, interventionThreshold)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Lomax"/> class.
        /// </summary>
        /// <param name="kappa">Kappa of Lomax distribution (shape parameter).</param>
        /// <param name="lambda">Lambda of Lomax distribution (scale parameter).</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public Lomax(double kappa, double lambda, double interventionThreshold = 600) : 
            base(kappa, lambda, interventionThreshold)
        {
            this.ShapeUpperBound = Math.Max(kappa * 4, 10);
            this.ScaleUpperBound = Math.Max(lambda * 4, 1);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Lomax"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">An array of durations waited before rebooting.</param>
        public Lomax(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations) : 
                    base(organicRecoveryDurations, inorganicRecoveryDurations)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Lomax"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">Organic recover time data.</param>
        /// <param name="inorganicRecoveryDurations">Inorganic recover time data.</param>
        /// <param name="organicRecoveryFeatureData">Organic recover feature data.</param>
        /// <param name="inorganicRecoveryFeatureData">Inorganic recover feature data.</param>
        public Lomax(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations,
            Matrix<double> organicRecoveryFeatureData,
            Matrix<double> inorganicRecoveryFeatureData) : 
            base(organicRecoveryDurations, inorganicRecoveryDurations, 
                organicRecoveryFeatureData, inorganicRecoveryFeatureData)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1;
        }

        public Lomax(
            List<double> organicRecoveryDurations) :
            base(organicRecoveryDurations)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1;
        }

        /// <summary>
        /// Calculates the optimum threshold.
        /// </summary>
        /// <param name="kappa">Kappa of Lomax distribution.</param>
        /// <param name="lambda">Lambda of Lomax distribution.</param>
        /// <param name="interventionCost">Cost of intervention.</param>
        /// <returns>The optimum intervention threshold.</returns>
        public static double GetOptimumThreshold(double kappa, double lambda, double interventionCost)
        {
            return (interventionCost * kappa) - (1 / lambda);
        }

        /// <summary>
        /// Calculates the expected downtime at a given intervention threshold.
        /// </summary>
        /// <param name="kappa">Kappa of Lomax distribution.</param>
        /// <param name="lambda">Lambda of Lomax distribution.</param>
        /// <param name="interventionCost">Cost of intervention.</param>
        /// <param name="interventionThreshold">Intervention threshold</param>
        /// <returns>The expected downtime at given intervention threshold.</returns>
        public static double GetExpectedDowntime(double kappa, double lambda, double interventionCost, 
            double interventionThreshold)
        {
            var cdf = 1 - Math.Pow(1 + (lambda * interventionThreshold), -kappa);

            return (cdf / (lambda * (kappa - 1))) + ((-interventionThreshold * kappa / (kappa - 1)) * (1 - cdf)) +
                   ((interventionThreshold + interventionCost) * (1 - cdf));
        }

        /// <summary>
        /// Calculates the expected downtime reduction at the optimum threshold.
        /// </summary>
        /// <param name="kappa">Kappa of Lomax distribution.</param>
        /// <param name="lambda">Lambda of Lomax distribution.</param>
        /// <param name="currentInterventionThreshold">The current intervention threshold.</param>
        /// <param name="optimumInterventionThreshold">The optimum intervention threshold.</param>
        /// <param name="interventionCost">The cost of intervention.</param>
        /// <returns>The downtime reduction at optimum threshold.</returns>
        public static double GetDowntimeReduction(
            double kappa,
            double lambda,
            double currentInterventionThreshold,
            double optimumInterventionThreshold,
            double interventionCost)
        {
            var currentDowntime = GetExpectedDowntime(kappa, lambda, interventionCost, currentInterventionThreshold);
            var optimumDowntime = GetExpectedDowntime(kappa, lambda, interventionCost, optimumInterventionThreshold);

            return currentDowntime - optimumDowntime;
        }

        /// <summary>
        /// Probability density function of the Lomax.
        /// </summary>
        /// <param name="x">The value at which PDF is to be evaluated.</param>
        /// <param name="k">Kappa (shape) of Lomax distribution.</param>
        /// <param name="lmb">Lambda (scale) of Lomax distribution.</param>
        /// <returns>The density function.</returns>
        public double PDF(double x, double k, double lmb)
        {
            return lmb * k / Math.Pow(1 + (lmb * x), k + 1);
        }

        /// <summary>
        /// Probability density function of Lomax.
        /// </summary>
        /// <param name="x">The value at which the density is to be calculated.</param>
        /// <returns>The value of the density function.</returns>
        public override double PDF(double x)
        {
            return this.PDF(x, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// Calculates the Log pdf of the Lomax distribution.
        /// </summary>
        /// <param name="t">The number at which the log-survival is to be calculated.</param>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>A double which is the log-pdf.</returns>
        public override double LogPdf(double t, double k, double lmb)
        {
            return (Math.Log(k) + Math.Log(lmb)) - ((k + 1) * MathUtilities.Log1P(lmb * t));
        }

        /// <summary>
        /// Calculates the inverse CDF for the Lomax distribution.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <param name="k">The shape parameter for the Lomax distribution.</param>
        /// <param name="lmb">The scale parameter for the Lomax distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public double InverseCDF(double p, double k, double lmb)
        {
            return (Math.Pow(1 / (1 - p), 1 / k) - 1) / lmb;
        }

        /// <summary>
        /// A wrapper for the inverse CDF for the Lomax distribution that uses the shape and scale for the instance.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public override double InverseCDF(double p)
        {
            return this.InverseCDF(p, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// A wrapper for the survival function.
        /// </summary>
        /// <param name="x">The value at which the survival is to be calculated.</param>
        /// <returns>A double which is the survival function at the given value.</returns>
        public override double Survival(double x)
        {
            return Math.Exp(this.LogSurvival(x, this.Kappa, this.Lambda));
        }

        /// <summary>
        /// Calculates the Log survival of the Lomax distribution.
        /// </summary>
        /// <param name="t">The number at which the log-survival is to be calculated.</param>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>A double which is the log-survival.</returns>
        public override double LogSurvival(double t, double k, double lmb)
        {
            return -k * MathUtilities.Log1P(lmb * t);
        }

        /// <summary>
        /// Calculate the gradient of the log of PDF at a given point.
        /// </summary>
        /// <param name="t">The point we want to evaluate the gradient.</param>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>The gradient of the log PDF at t.</returns>
        public override Vector<double> GradLPDF(double t, double k, double lmb)
        {
            double delK = (1 / k) - Math.Log(1 + (lmb * t));
            double delLmb = (1 / lmb) - ((k + 1) * t) / (1 + (lmb * t));
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculate the gradient of log of the survival function at a given point.
        /// </summary>
        /// <param name="x">The point to evaluate the gradient.</param>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>The gradient of survival function at x.</returns>
        public override Vector<double> GradLSurvival(double x, double k, double lmb)
        {
            double delK = -Math.Log(1 + lmb * x);
            double delLmb = -k * x / (1 + lmb * x);
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculates the gradient of the Log-likelihood function for the Lomax without features.
        /// </summary>
        /// <param name="k">The shape parameter of the Lomax distribution</param>
        /// <param name="lmb">The scale parameter of the Lomax distribution</param>
        /// <returns>A double array - the gradients for kappa and lambda</returns>
        public override Vector<double> GradLL(double k, double lmb)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;

            double delK = (n / k) - t.Sum(ti => MathUtilities.Log1P(lmb * ti));

            if (x.Count > 0)
            {
                delK += -x.Sum(xi => MathUtilities.Log1P(lmb * xi));
            }

            double delLmb = (n / lmb) - ((k + 1) * t.Sum(ti => ti / (1 + (lmb * ti))));

            if (x.Count > 0)
            {
                delLmb  += - (k * x.Sum(xi => xi / (1 + (lmb * xi))));
            }

            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Hessian for the log-likelihood of the Lomax distribution. Useful in Newton Rhapson method.
        /// To be used to verify the results from bisection method. Or when it fails.
        /// </summary>
        /// <param name="k">The shape parameter of the Lomax</param>
        /// <param name="lmb">The scale parameter of the Lomax</param>
        /// <returns>The matrix which is the hessian (matrix of second derivatives)</returns>
        public Matrix<double> HessianLL(double k, double lmb)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;

            double delKSq = -n / Math.Pow(k, 2);

            double delLmbSq = -(n / Math.Pow(lmb, 2)) + ((k + 1) * t.Sum(ti => Math.Pow(ti / (1 + (lmb * ti)), 2)))
                + (x == null ? 0 : (k * x.Sum(xi => Math.Pow(xi / (1 + (lmb * xi)), 2))));

            double delLmbDelK = -t.Sum(ti => ti / (1 + (lmb * ti))) - (x == null ? 0 : x.Sum(xi => xi 
                                / (1 + (lmb * xi))));

            Matrix<double> hessian = Matrix<double>.Build.Dense(2, 2);
            hessian[0, 0] = delKSq;
            hessian[1, 1] = delLmbSq;
            hessian[0, 1] = delLmbDelK;
            hessian[1, 0] = delLmbDelK;

            return hessian;
        }

        /// <summary>
        /// Derivative of expected downtime.
        /// </summary>
        /// <param name="interventionThreshold">Intervention threshold.</param>
        /// <param name="interventionCost">Intervention cost.</param>
        /// <returns>The derivative of expected downtime.</returns>
        public double DerivativeExpectedDowntime(double interventionThreshold, double interventionCost)
        {
            const double Eps = 1e-5;
            var rs = GetExpectedDowntime(this.Kappa, this.Lambda, interventionCost, interventionThreshold + Eps);
            var ls = GetExpectedDowntime(this.Kappa, this.Lambda, interventionCost, interventionThreshold - Eps);

            return (rs - ls) / (2 * Eps);
        }

        /// <summary>
        /// Calculate metrics for a given intervention cost.
        /// </summary>
        /// <param name="interventionCost">Intervention cost.</param>
        public void CalculateMetrics(double interventionCost)
        {
            this.OptimumThreshold = GetOptimumThreshold(this.Kappa, this.Lambda, interventionCost);
            this.DowntimeReduction = GetDowntimeReduction(
                this.Kappa,
                this.Lambda,
                this.InterventionThreshold,
                this.OptimumThreshold,
                interventionCost);
        }

        /// <summary>
        /// Optimization function.
        /// </summary>
        /// <param name="lambda">Lambda parameter.</param>
        /// <returns>Derivative of Kappa(Lambda).</returns>
        public double DelKappaFunc(double lambda)
        {
            var k1 = this.KappaFunc(lambda);

            double n = this.OrganicRecoveryCount;
            var term1 = this.OrganicRecoveryDurations.Sum(t => t/(1+lambda*t));
            var term2 = 0.0;

            if (this.InorganicRecoveryDurations != null)
            {
                term2 = this.InorganicRecoveryDurations.Sum(x => x / (1 + lambda * x));
            }

            var k2 = (n - k1 * term2) / lambda / term1-1;

            return (k1 - k2);
        }

        /// <summary>
        /// Kappa(Lambda) function.
        /// </summary>
        /// <param name="lambda">Lambda parameter.</param>
        /// <returns>Kappa parameter.</returns>
        public double KappaFunc(double lambda)
        {
            double n = this.OrganicRecoveryCount;

            var term1 = this.OrganicRecoveryDurations.Sum(t => MathUtilities.Log1P(lambda * t));
            var term2 = 0.0;

            if (this.InorganicRecoveryDurations != null)
            {
                term2 = this.InorganicRecoveryDurations.Sum(x => MathUtilities.Log1P(lambda * x));
            }

            return n / (term1 + term2);
        }
        
        /// <summary>
        /// Fits the model data with and without regularizers.
        /// </summary>
        /// <param name="alphaParam">Regularizer strength.</param>
        /// <param name="kappa0Param">Start point for Kappa.</param>
        /// <param name="lambda0Param">Start point for Lambda.</param>
        public void FitParameters(double alphaParam = 0, double kappa0Param = 0,
            double lambda0Param = 0)
        {
            // First assign 0 to kappa and lambda as a sign of this method being called.
            this.Kappa = 0.0;
            this.Lambda = 0.0;

            this.alpha = alphaParam;
            this.kappa0 = kappa0Param;
            this.lambda0 = lambda0Param;

            // Try to find the root, if cannot assign Lambda and Kappa to -1.
            try
            {
                this.Lambda = Bisection.FindRoot(this.DelKappaFunc, 0.0001, 1000.0, 1e-6, 10000);
                this.Kappa = this.KappaFunc(this.Lambda);
            }
            catch (Exception e)
            {
                string s = e.Message;
                this.Kappa = -1.0;
                this.Lambda = -1.0;
            }
        }
    }
}
