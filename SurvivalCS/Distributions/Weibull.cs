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
    public class Weibull : BaseModel
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Weibull"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="interventionCount">Number of interventions.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public Weibull(
            List<double> organicRecoveryDurations,
            double interventionCount,
            double interventionThreshold) : base(organicRecoveryDurations, interventionCount, interventionThreshold)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Weibull"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">An array of durations waited before rebooting.</param>
        public Weibull(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations) : base(organicRecoveryDurations, inorganicRecoveryDurations)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Weibull"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">An array of durations of inorganic transitions to 
        /// rebooting states.</param>
        /// <param name="organicRecoveryFeatureData">The features corresponding to samples of organic 
        /// recovery in a matrix.
        /// The number of rows should be the same as organicRecoveryDurations.</param>
        /// <param name="inorganicRecoveryFeatureData">The features corresponding to samples of inorganic 
        /// transitions to rebooting in a matrix.
        /// The number of rows should be the same as inorganicRecoveryDurations.</param>
        public Weibull(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations,
            Matrix<double> organicRecoveryFeatureData,
            Matrix<double> inorganicRecoveryFeatureData) : base(organicRecoveryDurations, 
                inorganicRecoveryDurations, organicRecoveryFeatureData, inorganicRecoveryFeatureData)
        {
            this.ShapeUpperBound = 10;
            this.ScaleUpperBound = 1000;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Weibull"/> class. Sets only the shape and scale
        /// without filling in any data.
        /// </summary>
        /// <param name="kappa">The shape parameter of the Weibull</param>
        /// <param name="lambda">The scale parameter of the Weibull</param>
        public Weibull(double kappa, double lambda) : base(kappa, lambda)
        {
            this.ShapeUpperBound = Math.Max(kappa * 4, 10);
            this.ScaleUpperBound = Math.Max(lambda * 4, 1000);
        }

        /// <summary>
        /// Calculates the optimal threshold based on Weibull distribution.
        /// </summary>
        /// <param name="kappa">The Kappa of the Weibull distribution</param>
        /// <param name="lambda">The Lambda of the Weibull distribution.</param>
        /// <param name="interventionCost">The reboot cost</param>
        /// <returns>The optimum Weibull threshold.</returns>
        public static double GetOptimumThreshold(double kappa, double lambda, double interventionCost)
        {
            return (Math.Pow(lambda, kappa / (kappa - 1)) / Math.Pow(interventionCost * kappa, 1 / (kappa - 1)));
        }

        /// <summary>
        /// The probability density function of Weibull.
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <returns>The Weibull density.</returns>
        public double PDF(double x, double k, double lmb)
        {
            return (k / lmb) * Math.Pow(x / lmb, k - 1) * Math.Exp(-Math.Pow(x / lmb, k));
        }

        /// <summary>
        /// Gets the PDf of the Weibull distribution.
        /// </summary>
        /// <param name="t">The value at which the PDF is to be calculated.</param>
        /// <returns>The Weibull PDF.</returns>
        public override double PDF(double t)
        {
            return this.PDF(t, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// The logarithm of the PDF of the Weibull distribution.
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <returns>The logrithm of the PDF at x.</returns>
        public override double LogPdf(double x, double k, double lmb)
        {
            return Math.Log(k) - (k * Math.Log(lmb)) + ((k - 1) * Math.Log(x)) - Math.Pow(x / lmb, k);
        }

        /// <summary>
        /// The gradient of the Pdf function.
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <param name="pdf">If pre-computed PDF is available, pass it here to save some computation.</param>
        /// <returns>Gradient of the pdf function (double array with 2 elements).</returns>
        public Vector<double> GradPDF(double x, double k, double lmb, double pdf = -1)
        {
            // If the pdf is not provided, calculate it. Otherwise, we can save the computation.
            if (pdf == -1)
            {
                pdf = this.PDF(x, k, lmb);
            }

            double delK = pdf * (((-Math.Pow(x / lmb, k) + 1) * Math.Log(x / lmb)) + (1 / k));
            double delLmb = (1 - Math.Pow(x / lmb, k)) * (-k / lmb) * pdf;
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculates the inverse CDF for the Weibull distribution.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <param name="k">The shape parameter for the Weibull distribution.</param>
        /// <param name="lmb">The scale parameter for the Weibull distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public double InverseCDF(double p, double k, double lmb)
        {
            return lmb * Math.Pow(-Math.Log(1 - p), 1 / k);
        }

        /// <summary>
        /// A wrapper for the inverse CDF for the Weibull distribution that uses the shape and scale for the instance.
        /// </summary>
        /// <param name="p">A number between 0 and 1, typically generated from a uniform distribution.</param>
        /// <returns>The value of the inverse CDF at p.</returns>
        public override double InverseCDF(double p)
        {
            return this.InverseCDF(p, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// The survival function for Weibull.
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <returns>The value of the survival function.</returns>
        public double Survival(double x, double k, double lmb)
        {
            return Math.Exp(-Math.Pow(x / lmb, k));
        }

        /// <summary>
        /// Wrapper for survival function.
        /// </summary>
        /// <param name="x">The value at which the survival function is to be calculated.</param>
        /// <returns>A double which is the value of the survival function.</returns>
        public override double Survival(double x)
        {
            return this.Survival(x, this.Kappa, this.Lambda);
        }

        /// <summary>
        /// The logarithm of the survival function.
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <returns>The logrithm of the survival function at x</returns>
        public override double LogSurvival(double x, double k, double lmb)
        {
            return -Math.Pow(x / lmb, k);
        }

        /// <summary>
        /// The gradient of the survival function
        /// </summary>
        /// <param name="x">Value at function evaluation.</param>
        /// <param name="k">Shape parameter.</param>
        /// <param name="lmb">Scale parameter.</param>
        /// <param name="survive">If the survival function is pre-computed, pass it here to save some computation.</param>
        /// <returns>Gradient of the survival function (double array with 2 elements).</returns>
        public Vector<double> GradSurvival(double x, double k, double lmb, double survive = -1)
        {
            // If survival is provided, use it and save some computations. Otherwise, calculate it.
            if (survive == -1)
            {
                survive = this.Survival(x, k, lmb);
            }

            double delK = -survive * Math.Pow(x / lmb, k) * Math.Log(x / lmb);
            double delLmb = survive * Math.Pow(x / lmb, k) * (k / lmb);
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculates the gradient of the Log-likelihood function for the Weibull without features.
        /// </summary>
        /// <param name="k">The shape parameter of the Weibull distribution.</param>
        /// <param name="lmb">The scale parameter of the Weibull distribution.</param>
        /// <returns>A double array - the gradients for kappa and lambda.</returns>
        public override Vector<double> GradLL(double k, double lmb)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            double delK = (n / k) - (n * Math.Log(lmb)) + t.Sum(ti => Math.Log(ti)) - 
                t.Sum(ti => Math.Pow(ti / lmb, k) * Math.Log(ti / lmb)) - 
                x.Sum(xi => Math.Pow(xi / lmb, k) * Math.Log(xi / lmb));

            double delLmb = ((-n * k) / lmb) + ((k / Math.Pow(lmb, k + 1)) * 
                t.Sum(ti => Math.Pow(ti, k))) + x.Sum(xi => Math.Pow(xi, k));

            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculates the gradient of the log-likelihood function for Weibull with features.
        /// Provided below is the derivation in Latex. Paste it into an online 
        /// Latex editor like https://www.overleaf.com/7650945xdjrrdzrmsdd#/26769609/
        /// to read and and understand where the formulas came from.
        /// Let $f_i$ be the vector of features for the $i^{th}$ data point and $x_i$ be the duration it took to recover. 
        /// Then, the Likelihood and log-likelihood functions are given by - 
        /// \[L = \prod_{i=1}^n \log(pdf(x_i, W.f_i)) \]
        /// \[ll = \sum_{i=1}^n \log(pdf(x_i, W.f_i))\]
        /// \[\frac{\partial(ll)}{\partial W} = \frac{1}{pdf(x_i, W.f_i)} f_i(\partial \theta_i)^T \]
        /// Where,
        /// \[\theta_i = W.f_i  = [\kappa, \lambda]\] here $\kappa$, once the sigmoid is applied to it 
        /// to ensure it being positive, is the shape paramter of the Weibull and $\lambda$ once the 
        /// sigmoid is applied is the scale parameter.
        /// \[\partial \theta_i = \frac{\partial(pdf(x,\theta_i))}{\partial(\theta_i)}\]
        /// </summary>
        /// <param name="w">The matrix of parameters (2 x # of features)</param>
        /// <param name="fSamples">The features corresponding to the organic recoveries.
        /// Number of rows should be same as this.OrganicRecoveryDurations.Length</param>
        /// <param name="fCensored">The features corresponding to the reboots.
        /// Number of rows should be the same as this.InorganicRecoveryDurations.Length</param>
        /// <param name="eps">Since we divide by the PDF and survival functions, we need to make sure
        /// they don't get lower than a threshold or the gradients blow up. 
        /// This parameter is a lower bound on them.</param>
        /// <param name="bailOutSurvivalValue">The value to turn to for survival when gradients start blowing up.</param>
        /// <returns>The gradient of the log-likelihood with the matrix of parameters, w</returns>
        public Matrix<double> GradLL(Matrix<double> w, Matrix<double> fSamples, Matrix<double> fCensored, double eps = 1e-8, double bailOutSurvivalValue = 10.0)
        {
            IReadOnlyCollection<double> t = this.OrganicRecoveryDurations;
            IReadOnlyCollection<double> x = this.InorganicRecoveryDurations;
            Matrix<double> gradW = Matrix<double>.Build.Dense(w.RowCount, w.ColumnCount);
            Sigmoid sShape = new Sigmoid(this.ShapeUpperBound);
            Sigmoid sScale = new Sigmoid(this.ScaleUpperBound);

            for (int i = 0; i < fSamples.RowCount; i++)
            {
                Vector<double> currentRow = fSamples.Row(i);
                Vector<double> theta = w.Multiply(currentRow); // A 2 dim vector that will be converted to the 2 Weibull params.
                double shape = sShape.Transform(theta[0]); // To prevent the parameters from becoming negative, we use sigmoids.
                double scale = sScale.Transform(theta[1]);
                double pdf = this.PDF(t.ElementAt(i), shape, scale);
                Vector<double> pdfGrad = this.GradPDF(t.ElementAt(i), shape, scale, pdf);
                Vector<double> sigmoidGrad = Vector<double>.Build.DenseOfArray(new double[] { sShape.Grad(theta[0]), sScale.Grad(theta[1]) });
                Vector<double> delTheta = pdfGrad.PointwiseMultiply(sigmoidGrad); // Since we used the sigmoids, the product rule dictates that we point-wise multiply the derivatives.

                // Since we are dividing the matrix by the pdf, we need to be careful it doesn't blow up.
                if (pdf > eps)
                {
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(pdf)); // del/del(W) log(lik(x,W.f)) = 1/lik(x,W.f) * f. (del(lik(x))/del(W.f))^T
                }
                else
                {
                    double survival = this.Survival(bailOutSurvivalValue, shape, scale);
                    Vector<double> survivalGrad = this.GradSurvival(bailOutSurvivalValue, shape, scale, survival);
                    delTheta = survivalGrad.PointwiseMultiply(sigmoidGrad);
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(survival));
                }
            }

            for (int i = 0; i < fCensored.RowCount; i++)
            {
                Vector<double> currentRow = fCensored.Row(i);
                Vector<double> theta = w.Multiply(currentRow);
                double shape = sShape.Transform(theta[0]);
                double scale = sScale.Transform(theta[1]);
                double survival = this.Survival(x.ElementAt(i), shape, scale);
                Vector<double> survivalGrad = this.GradSurvival(x.ElementAt(i), shape, scale, survival);
                Vector<double> sigmoidGrad = Vector<double>.Build.DenseOfArray(new double[] { sShape.Grad(theta[0]), sScale.Grad(theta[1]) });
                Vector<double> delTheta = survivalGrad.PointwiseMultiply(sigmoidGrad);

                // Since we are dividing the matrix by the survival, we need to be careful it doesn't blow up.
                if (survival > eps)
                {
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(survival));
                }
                else
                {
                    survival = this.Survival(bailOutSurvivalValue, shape, scale);
                    survivalGrad = this.GradSurvival(bailOutSurvivalValue, shape, scale, survival);
                    delTheta = survivalGrad.PointwiseMultiply(sigmoidGrad);
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(survival));
                }
            }

            return gradW;
        }

        /// <summary>
        /// Hessian for the log-likelihood of the Weibull distribution. Useful in Newton Raphson method.
        /// To be used to verify the results from bisection method. Or when it fails.
        /// </summary>
        /// <param name="k">The shape parameter of the Weibull</param>
        /// <param name="lmb">The scale parameter of the Weibull</param>
        /// <returns>The matrix which is the hessian (matrix of second derivatives)</returns>
        public Matrix<double> HessianLL(double k, double lmb)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;

            double delKSq = -(n / Math.Pow(k, 2)) - t.Sum(ti => Math.Pow(ti / lmb, k) * Math.Pow(Math.Log(ti / lmb), 2))
                - x.Sum(xi => Math.Pow(xi / lmb, k) * Math.Pow(Math.Log(xi / lmb), 2));

            double delLmbSq = ((n * k) / Math.Pow(lmb, 2)) +
                ((t.Sum(ti => Math.Pow(ti, k)) + x.Sum(xi => Math.Pow(xi, k))) * (-(k * (k + 1)) / Math.Pow(lmb, k + 2)));

            double delLmbDelK = -(n / lmb) + (
                (1 / lmb) *
                (t.Sum(ti => ((k * Math.Pow(ti / lmb, k)) * Math.Log(ti / lmb)) + Math.Pow(ti / lmb, k))
               + x.Sum(xi => ((k * Math.Pow(xi / lmb, k)) * Math.Log(xi / lmb)) + Math.Log(xi / lmb, k))));

            Matrix<double> hessian = Matrix<double>.Build.Dense(2, 2);
            hessian[0, 0] = delKSq;
            hessian[1, 1] = delLmbSq;
            hessian[0, 1] = delLmbDelK;
            hessian[1, 0] = delLmbDelK;
            return hessian;
        }

        /// <summary>
        /// Calculate metrics for a given intervention cost.
        /// </summary>
        /// <param name="interventionCost">Intervention cost.</param>
        public void CalculateMetrics(double interventionCost)
        {
            this.OptimumThreshold = GetOptimumThreshold(this.Kappa, this.Lambda, interventionCost);
        }

        /// <summary>
        /// Fits the model data with and without regularizers.
        /// </summary>
        /// <param name="alphaParam">Regularizer strength.</param>
        /// <param name="kappa0Param">Start point for Kappa.</param>
        /// <param name="lambda0Param">Start point for Lambda.</param>
        public void FitParameters(double alphaParam = 0, double kappa0Param = 0, double lambda0Param = 0)
        {
            // First assign 0 to kappa and lambda as a sign of this method being called.
            this.Kappa = 0.0;
            this.Lambda = 0.0;

            // Try to find the root, if cannot assign Lambda and Kappa to -1.
            try
            {
                this.Kappa = Bisection.FindRoot(this.DelKappa, 0.0001, 20.0, 1e-6, 10000);
                this.Lambda = this.OptimumLambda(this.Kappa);
            }
            catch
            {
                this.Kappa = -1.0;
                this.Lambda = -1.0;
            }
        }

        /// <summary>
        /// Calculate the gradient of the log of PDF at a given point.
        /// </summary>
        /// <param name="t">The point to evaluate the gradient.</param>
        /// <param name="k">The shape parameter</param>
        /// <param name="lmb">The scale parameter</param>
        /// <returns>The gradient of log PDF at t.</returns>
        public override Vector<double> GradLPDF(double t, double k, double lmb)
        {
            double delK = (1 / k) - Math.Log(lmb) + Math.Log(t) - Math.Pow(t / lmb, k) * Math.Log(t / lmb);
            double delLmb = -k / lmb + (k * Math.Pow(t, k)) / Math.Pow(lmb, k + 1);
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Calculate the gradient of log of the survival function at a given point.
        /// </summary>
        /// <param name="x">The point to evaluate the gradient.</param>       
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>The gradient of the log survival function at x.</returns>
        public override Vector<double> GradLSurvival(double x, double k, double lmb)
        {
            double delK = -Math.Pow(x / lmb, k) * Math.Log(x / lmb);
            double delLmb = (k * Math.Pow(x, k)) / Math.Pow(lmb, k + 1);
            Vector<double> res = Vector<double>.Build.DenseOfArray(new double[] { delK, delLmb });
            return res;
        }

        /// <summary>
        /// Kappa function for Weibull. This is what will be set to zero through bisection.
        /// </summary>
        /// <param name="kappa">The value of kappa at which the derivative is calculated.</param>
        /// <returns>The derivative of the log likelihood for Weibull distribution w.r.t. kappa param.</returns>
        private double DelKappa(double kappa)
        {
            var logsum = this.OrganicRecoveryDurations.Sum(t => Math.Log(t));
            var n = this.OrganicRecoveryCount;
            var tiPowklogx = this.OrganicRecoveryDurations.Sum(t => Math.Pow(t, kappa) * Math.Log(t));
            var tiPowk = this.OrganicRecoveryDurations.Sum(t => Math.Pow(t, kappa));

            if (this.InorganicRecoveryDurations == null)
            {
                var m = this.InorganicRecoveryCount;
                var tau0 = this.InterventionThreshold;
                return (n / kappa) + logsum - ((n * (tiPowklogx + ((m * Math.Pow(tau0, kappa)) * Math.Log(tau0)))) / (tiPowk + (m * Math.Pow(tau0, kappa))));
            }
            else
            {
                var tau0PowkLogtau0 = this.InorganicRecoveryDurations.Sum(t => Math.Pow(t, kappa) * Math.Log(t));
                var tau0Powk = this.InorganicRecoveryDurations.Sum(t => Math.Pow(t, kappa));
                return (n / kappa) + logsum - ((n * (tiPowklogx + tau0PowkLogtau0)) / (tiPowk + tau0Powk));
            }
        }

        /// <summary>
        /// Lambda function for Weibull given kappa.
        /// </summary>
        /// <param name="kappa">The kappa at which lambda is to be calculated.</param>
        /// <returns>The value of the scale parameter, lambda of the Weibull distribution.</returns>
        private double OptimumLambda(double kappa)
        {
            var tiPowk = this.OrganicRecoveryDurations.Sum(t => Math.Pow(t, kappa));
            var n = this.OrganicRecoveryCount;

            if (this.InorganicRecoveryDurations == null)
            {
                var m = this.InorganicRecoveryCount;
                var tau0 = this.InterventionThreshold;
                return Math.Pow((tiPowk + (m * Math.Pow(tau0, kappa))) / n, 1 / kappa);
            }
            else
            {
                var tau0Powk = this.InorganicRecoveryDurations.Sum(t => Math.Pow(t, kappa));
                return Math.Pow((tiPowk + tau0Powk) / n, 1 / kappa);
            }
        }
    }
}
