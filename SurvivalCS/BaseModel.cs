using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SurvivalCS.Functions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace SurvivalCS
{
    public abstract class BaseModel
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="BaseModel"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryCount">Number of interventions.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public BaseModel(
            List<double> organicRecoveryDurations,
            double inorganicRecoveryCount,
            double interventionThreshold)
        {
            this.OrganicRecoveryDurations = organicRecoveryDurations; // data
            this.OrganicRecoveryCount = organicRecoveryDurations.Count; // n
            this.InorganicRecoveryCount = inorganicRecoveryCount; // m
            this.InterventionThreshold = interventionThreshold; // tau0
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseModel"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">The list of thresholds Fabric actually waited before rebooting.</param>
        public BaseModel(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations)
        {
            this.OrganicRecoveryDurations = organicRecoveryDurations; // data
            this.OrganicRecoveryCount = organicRecoveryDurations.Count; // n
            this.InorganicRecoveryDurations = inorganicRecoveryDurations; // tau0_j
            this.InorganicRecoveryCount = inorganicRecoveryDurations.Count; // m
        }

        public BaseModel(List<double> organicRecovveryDurations)
        {
            this.OrganicRecoveryDurations = organicRecovveryDurations;
            this.OrganicRecoveryCount = organicRecovveryDurations.Count; //n
            this.InorganicRecoveryDurations = new List<double>();
            this.InorganicRecoveryCount = 0;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseModel"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">The list of thresholds Fabric actually waited before rebooting.</param>
        /// <param name="organicRecoveryDurationsByContainer">Organic recovery durations multiplied by container count</param>
        /// <param name="inorganicRecoveryCountByContainer">Inorganic recovery count multiplied by container count</param>
        public BaseModel(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations,
            List<double> organicRecoveryDurationsByContainer,
            double inorganicRecoveryCountByContainer)
        {
            this.OrganicRecoveryDurations = organicRecoveryDurations; // data
            this.OrganicRecoveryCount = organicRecoveryDurations.Count; // n
            this.InorganicRecoveryDurations = inorganicRecoveryDurations; // tau0_j
            this.InorganicRecoveryCount = inorganicRecoveryDurations.Count; // m
            this.OrganicRecoveryDurationsByContainer = organicRecoveryDurationsByContainer;
            this.InorganicRecoveryCountByContainer = inorganicRecoveryCountByContainer;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseModel"/> class.
        /// </summary>
        /// <param name="kappa">The shape parameter.</param>
        /// <param name="lambda">The scale parameter.</param>
        /// <param name="interventionThreshold">The current threshold for interventions.</param>
        public BaseModel(double kappa, double lambda, double interventionThreshold = 600)
        {
            this.Kappa = kappa;
            this.Lambda = lambda;
            this.InterventionThreshold = interventionThreshold;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BaseModel"/> class.
        /// </summary>
        /// <param name="organicRecoveryDurations">An array of durations of transitions.</param>
        /// <param name="inorganicRecoveryDurations">The list of thresholds Fabric actually waited before rebooting.</param>
        /// <param name="organicRecoveryFeatureData">The features corresponding to samples of organic recovery in a matrix.
        /// The number of rows should be the same as organicRecoveryDurations.</param>
        /// <param name="inorganicRecoveryFeatureData">The features corresponding to samples of inorganic transitions to rebooting in a matrix.
        /// The number of rows should be the same as inorganicRecoveryDurations.</param>
        public BaseModel(
            List<double> organicRecoveryDurations,
            List<double> inorganicRecoveryDurations,
            Matrix<double> organicRecoveryFeatureData,
            Matrix<double> inorganicRecoveryFeatureData)
        {
            this.OrganicRecoveryDurations = organicRecoveryDurations; // data
            this.OrganicRecoveryCount = organicRecoveryDurations.Count; // n
            this.InorganicRecoveryDurations = inorganicRecoveryDurations; // tau0_j
            this.InorganicRecoveryCount = inorganicRecoveryDurations.Count; // m
            this.OrganicRecoveryFeatureData = organicRecoveryFeatureData; // f_samples
            this.InorganicRecoveryFeatureData = inorganicRecoveryFeatureData; // f_censored
        }

        /// <summary>
        /// Gets or sets Kappa (shape parameter - if any) of given distribution.
        /// </summary>
        public double Kappa { get; set; } = 1;

        /// <summary>
        /// Gets or sets Lambda (scale parameter - if any) of given distribution.
        /// </summary>
        public double Lambda { get; set; } = 1;

        /// <summary>
        /// Gets or sets the maximum value the shape parameter can take.
        /// Only valid for model with features.
        /// </summary>
        public double ShapeUpperBound { get; set; } = 10;

        /// <summary>
        /// Gets or sets the maximum value the scale parameter can take.
        /// Only valid for model with features.
        /// </summary>
        public double ScaleUpperBound { get; set; } = 1;

        /// <summary>
        /// Gets or sets optimum threshold.
        /// </summary>
        public double OptimumThreshold { get; set; }

        /// <summary>
        /// Gets or sets downtime reduction at optimum threshold compared to current threshold.
        /// </summary>
        public double DowntimeReduction { get; set; }

        /// <summary>
        /// Gets durations.
        /// </summary>
        public List<double> OrganicRecoveryDurations { get; }

        /// <summary>
        /// Gets durations by container.
        /// </summary>
        public List<double> OrganicRecoveryDurationsTimesContainerCount { get; }

        /// <summary>
        /// Gets durations by container. Does this by adding each duration n times if 
        /// there are n containers on the node.
        /// </summary>
        public List<double> OrganicRecoveryDurationsByContainer { get; }

        /// <summary>
        /// Gets actual times waited before reboot.
        /// </summary>
        public List<double> InorganicRecoveryDurations { get; }

        /// <summary>
        /// Gets the count of durations.
        /// </summary>
        public double OrganicRecoveryCount { get; }

        /// <summary>
        /// Gets the intervention counts.
        /// </summary>
        public double InorganicRecoveryCount { get; }

        /// <summary>
        /// Gets the intervention counts by container.
        /// </summary>
        public double InorganicRecoveryCountByContainer { get; }

        /// <summary>
        /// Gets reboot durations by container. Does this by adding each duration n times 
        /// if there are n containers on the node.
        /// </summary>
        public double InorganicRecoveryDurationsByContainer { get; }

        /// <summary>
        /// Gets the current intervention threshold.
        /// </summary>
        public double InterventionThreshold { get; }

        /// <summary>
        /// Gets the feature data corresponding to organic recoveries.
        /// </summary>
        public Matrix<double> OrganicRecoveryFeatureData { get; }

        /// <summary>
        /// Gets the feature data corresponding to inorganic reboots.
        /// </summary>
        public Matrix<double> InorganicRecoveryFeatureData { get; }

        
        /// <summary>
        /// Every distribution has a Probability Density Function.
        /// </summary>
        /// <param name="x">The value at which the density function is to be evaluated.</param>
        /// <returns>A double, the probability density.</returns>
        public abstract double PDF(double x);

        /// <summary>
        /// The log of the PDF of the distribution.
        /// </summary>
        /// <param name="x">The value at which the log pdf is to be calculated.</param>
        /// <param name="shape">The shape parameter of the distribution.</param>
        /// <param name="scale">The scale parameter of the distribution.</param>
        /// <returns>A double which is the log of the pdf.</returns>
        public abstract double LogPdf(double x, double shape, double scale);

        /// <summary>
        /// The inverse of the CDF of the distribution if it exists.
        /// </summary>
        /// <param name="p">The value at which the inverse cdf is to be calculated. Must be between 0 and 1.</param>
        /// <returns>A double which is the inverse of the cdf.</returns>
        public abstract double InverseCDF(double p);

        /// <summary>
        /// Calculates the survival function, or the probability that the distribution will exceed a certain value.
        /// </summary>
        /// <param name="x">The value at which the survival function is to be calculated (usually time in seconds).</param>
        /// <returns>A double which is the survival function</returns>
        public abstract double Survival(double x);

        /// <summary>
        /// The log of the survival function of the distribution (area under curve beyond a certain value)
        /// </summary>
        /// <param name="x">The value at which the log pdf is to be calculated</param>
        /// <param name="shape">The shape parameter of the distribution.</param>
        /// <param name="scale">The scale parameter of the distribution</param>
        /// <returns>A double which is the log of the survival.</returns>
        public abstract double LogSurvival(double x, double shape, double scale);

        /// <summary>
        /// Calculate the gradient of the log of PDF at a given point.
        /// </summary>
        /// <param name="t">The point to evaluate the gradient.</param>
        /// <param name="k">The shape parameter</param>
        /// <param name="lmb">The scale parameter</param>
        /// <returns>The gradient of log PDF at t.</returns>
        public abstract Vector<double> GradLPDF(double t, double k, double lmb);

        /// <summary>
        /// Calculate the gradient of log of the survival function at a given point.
        /// </summary>
        /// <param name="x">The point to evaluate the gradient.</param>       
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>The gradient of the log survival function at x.</returns>
        public abstract Vector<double> GradLSurvival(double x, double k, double lmb);

        /// <summary>
        /// The log-likelihood of the Weibull distribution on censored and uncensored arrays with features.
        /// </summary>
        /// <param name="w">The matrix of parameters.</param>
        /// <param name="fSamples">The features corresponding to the organic recoveries.
        /// Number of rows should be same as this.OrganicRecoveryDurations.Length</param>
        /// <param name="fCensored">The features corresponding to the reboots.
        /// Number of rows should be the same as this.InorganicRecoveryDurations.Length</param>
        /// <returns>The log-likelihood of the data along with features.</returns>
        public double LogLikelihood(Matrix<double> w, Matrix<double> fSamples, 
                                    Matrix<double> fCensored)
        {
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            double lik = 0;
            Sigmoid sShape = new Sigmoid(this.ShapeUpperBound);
            Sigmoid sScale = new Sigmoid(this.ScaleUpperBound);

            for (int i = 0; i < fSamples.RowCount; i++)
            {
                Vector<double> currentRow = fSamples.Row(i);
                Vector<double> theta = w.Multiply(currentRow);
                double shape = sShape.Transform(theta[0]);
                double scale = sScale.Transform(theta[1]);
                lik += this.LogPdf(t.ElementAt(i), shape, scale);
            }

            for (int i = 0; i < fCensored.RowCount; i++)
            {
                Vector<double> currentRow = fCensored.Row(i);
                Vector<double> theta = w.Multiply(currentRow);
                double shape = sShape.Transform(theta[0]);
                double scale = sScale.Transform(theta[1]);
                lik += this.LogSurvival(x.ElementAt(i), shape, scale);
            }

            return lik;
        }

        /// <summary>
        /// This is a wrapper for LogLikelihood which fills in some values for the 
        /// feature sets for the sampled and censored data.
        /// Useful especially for validating gradients.
        /// </summary>
        /// <param name="w">The matrix of parameters.</param>
        /// <returns>The log likelihood as a double.</returns>
        public double LogLikelihood(Matrix<double> w)
        {
            if (this.OrganicRecoveryFeatureData == null || this.InorganicRecoveryFeatureData == null)
            {
                Tuple<Matrix<double>, Matrix<double>> features 
                    = this.PopulateFeatureMatrices(w.ColumnCount);
                return this.LogLikelihood(w, features.Item1, features.Item2);
            }
            else
            {
                return this.LogLikelihood(w, this.OrganicRecoveryFeatureData, 
                                             this.InorganicRecoveryFeatureData);
            }
        }

        /// <summary>
        /// Generates samples from the distribution.
        /// </summary>
        /// <param name="numSamples">The number of samples to be generated.</param>
        /// <returns>A list of doubles containing the generated samples.</returns>
        public List<double> GenerateSample(int numSamples)
        {
            Random r = new Random(Guid.NewGuid().GetHashCode());
            List<double> samples = new List<double>();

            for (int i = 0; i < numSamples; i++)
            {
                samples.Add(this.InverseCDF(r.NextDouble()));
            }

            return samples;
        }

        /// <summary>
        /// Calculating gradient of the loglikelihood using GradLPDF and GradLSurvival. 
        /// </summary>
        /// <param name="w">The current weight where we want to evaluate the gradient.</param>
        /// <param name="fSamples">Organic recover time data.</param>
        /// <param name="fCensored">Inorganic recover time data.</param>
        /// <param name="eps">Threshold for small pdf.</param>
        /// <param name="bailOutSurvivalValue">Threshold for small survival function.</param>
        /// <returns>The gradient of the log-likelihood function with respect to w (the weights).
        /// </returns>
        public Matrix<double> GradLL2(Matrix<double> w, Matrix<double> fSamples, 
                                    Matrix<double> fCensored, 
                                    double eps = 1e-8, double bailOutSurvivalValue = 10.0)
        {
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            Matrix<double> gradW = Matrix<double>.Build.Dense(w.RowCount, w.ColumnCount);
            Sigmoid sShape = new Sigmoid(this.ShapeUpperBound);
            Sigmoid sScale = new Sigmoid(this.ScaleUpperBound);

            for (int i = 0; i < fSamples.RowCount; i++)
            {
                Vector<double> currentRow = fSamples.Row(i);
                // A 2 dim vector that will be converted to the 2 Weibull params.
                Vector<double> theta = w.Multiply(currentRow);
                // To prevent the parameters from becoming negative, we use sigmoids.
                double shape = sShape.Transform(theta[0]);
                double scale = sScale.Transform(theta[1]);
                //// double pdf = this.PDF(t.ElementAt(i), kappa, lambda);

                Vector<double> lpdfGrad = this.GradLPDF(t.ElementAt(i), shape, scale);
                Vector<double> sigmoidGrad 
                    = Vector<double>.Build.DenseOfArray(new double[] {sShape.Grad(theta[0]),
                                                        sScale.Grad(theta[1])});
                // Since we used the sigmoids, the product rule dictates that 
                // we point-wise multiply the derivatives.
                Vector<double> delTheta = lpdfGrad.PointwiseMultiply(sigmoidGrad);

                gradW = gradW.Add(delTheta.OuterProduct(currentRow)); //// currentRow is just feature vector.

                if (double.IsNaN(gradW[0, 0]) || double.IsPositiveInfinity(gradW[0, 0])
                    || double.IsNegativeInfinity(gradW[0, 0]))
                {
                    // Hopefully, we will never enter this code path.
                    throw new Exception("The moment we feared has arrived, gradient has blown"+
                        "up due to samples data." +
                        "First, try tightening the upper bounds for the shape and scale parameters" +
                        "and if that doesn't work, add a break point here. My suspicion in delTheta");
                }

                /* Since we are dividing the matrix by the pdf, we need to be careful it doesn't blow up.
                if (pdf > eps)
                {
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(pdf));
                    // del/del(W) log(lik(x,W.f)) = 1/lik(x,W.f) * f. (del(lik(x))/del(W.f))^T
                }
                else
                {
                    double survival = this.Survival(bailOutSurvivalValue, kappa, lambda);
                    Vector<double> survivalGrad = this.GradSurvival(bailOutSurvivalValue, kappa, lambda, survival);
                    delTheta = survivalGrad.PointwiseMultiply(sigmoidGrad);
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(survival));
                }
                */
            }

            for (int i = 0; i < fCensored.RowCount; i++)
            {
                Vector<double> currentRow = fCensored.Row(i);
                Vector<double> theta = w.Multiply(currentRow);
                // To prevent the parameters from becoming negative, we use sigmoids.
                double shape = sShape.Transform(theta[0]); 
                double scale = sScale.Transform(theta[1]);
                //// double survival = this.Survival(t.ElementAt(i), kappa, lambda);

                Vector<double> lsurvivalGrad = this.GradLSurvival(x.ElementAt(i), shape, scale);
                Vector<double> sigmoidGrad = 
                    Vector<double>.Build.DenseOfArray(
                        new double[] { sShape.Grad(theta[0]), sScale.Grad(theta[1])});
                Vector<double> delTheta = lsurvivalGrad.PointwiseMultiply(sigmoidGrad);

                gradW = gradW.Add(delTheta.OuterProduct(currentRow));

                if (double.IsNaN(gradW[0, 0]) || double.IsPositiveInfinity(gradW[0, 0]) 
                    || double.IsNegativeInfinity(gradW[0, 0]))
                {
                    throw new Exception("The moment we feared has arrived, gradient has blown"+ 
                        "up due to censored data." +
                        "First, try tightening the upper bounds for the shape and scale parameters" +
                        "and if that doesn't work, add a break point here. My suspicion in delTheta");
                }

                /*
                Since we are dividing the matrix by the survival, we need to be careful it doesn't blow up.
                if (survival > eps)
                {
                }
                else
                {
                    survival = this.Survival(bailOutSurvivalValue, kappa, lambda);
                    survivalGrad = this.GradSurvival(bailOutSurvivalValue, kappa, lambda, survival);
                    delTheta = survivalGrad.PointwiseMultiply(sigmoidGrad);
                    gradW = gradW.Add(delTheta.OuterProduct(currentRow).Divide(survival));
                }
                */
            }

            return gradW;
        }

        /// <summary>
        /// This is a wrapper for GradLL which fills in some values for the feature sets for the sampled 
        /// and censored data.
        /// Useful especially for validating gradients.
        /// </summary>
        /// <param name="w">The matrix of parameters.</param>
        /// <returns>The gradient as a dense matrix.</returns>
        public Matrix<double> GradLL(Matrix<double> w)
        {
            if (this.OrganicRecoveryFeatureData == null || this.InorganicRecoveryFeatureData == null)
            {
                Tuple<Matrix<double>, Matrix<double>> features 
                    = this.PopulateFeatureMatrices(w.ColumnCount);
                return this.GradLL2(w, features.Item1, features.Item2);
            }
            else
            {
                return this.GradLL2(w, this.OrganicRecoveryFeatureData, 
                                    this.InorganicRecoveryFeatureData);
            }
        }

        /// <summary>
        /// The gradient of the log likelihood function with respect to shape and scale parameters.
        /// </summary>
        /// <param name="shape">The shape parameter of the distribution.</param>
        /// <param name="scale">The scale parameter of the distribution.</param>
        /// <returns>A two dimensional vector with the derivatives of likelihood function with 
        /// respect to scale and shape.</returns>
        public abstract Vector<double> GradLL(double shape, double scale);

        /// <summary>
        /// The log of the likelihood for the model given two parameters, typically shape and scale.
        /// </summary>
        /// <param name="shape">The shape parameter of the distribution.</param>
        /// <param name="scale">The scale parameter of the distribution.</param>
        /// <returns>A double which is the log of the likelihood function.</returns>
        public double LogLikelihood(double shape, double scale)
        {
            double n = this.OrganicRecoveryCount;
            List<double> t = this.OrganicRecoveryDurations;
            List<double> x = this.InorganicRecoveryDurations;
            return t.Sum(ti => this.LogPdf(ti, shape, scale)) 
                + x.Sum(xi => this.LogSurvival(xi, shape, scale));
        }

        /// <summary>
        /// Calculates the numerical gradient of the log likelihood with respect to 
        /// a matrix of parameters.
        /// </summary>
        /// <param name="w">The matrix of input parameters.</param>
        /// <param name="h">The delta by which the parameters will be perturbed for 
        /// finding the derivative.
        /// Should be small.</param>
        /// <returns>A matrix containing the values of the numerical derivative.</returns>
        public Matrix<double> NumericalGradLL(Matrix<double> w, double h = 1e-5)
        {
            Matrix<double> gradW = Matrix<double>.Build.Dense(w.RowCount, w.ColumnCount);
            Matrix<double> hMatrix = Matrix<double>.Build.Dense(w.RowCount, w.ColumnCount);

            for (int i = 0; i < w.RowCount; i++)
            {
                for (int j = 0; j < w.ColumnCount; j++)
                {
                    hMatrix[i, j] = h;
                    gradW[i, j] = (this.LogLikelihood(w + hMatrix) 
                                 - this.LogLikelihood(w - hMatrix)) / (2 * h);
                    hMatrix[i, j] = 0;
                }
            }

            return gradW;
        }

        /// <summary>
        /// Calculates the hessian, which is the matrix of second derivatices numerically 
        /// for models with two parameters.
        /// Assumes that loglikelihood has been implemented.
        /// </summary>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>A MathNet 2x2 matrix which is the second derivatives calculated 
        /// numerically.</returns>
        public Vector<double> NumericalGradLL(double k, double lmb)
        {
            double eps = 1e-5;
            Vector<double> grad = Vector<double>.Build.Dense(2);

            double dellmb = (this.LogLikelihood(k, lmb + eps) 
                    - this.LogLikelihood(k, lmb - eps)) / (2 * eps);
            double delk = (this.LogLikelihood(k + eps, lmb) 
                    - this.LogLikelihood(k - eps, lmb)) / (2 * eps);

            grad[0] = delk;
            grad[1] = dellmb;

            return grad;
        }

        /// <summary>
        /// Calculates the hessian, which is the matrix of second derivatices numerically for models with two parameters.
        /// Assumes that loglikelihood has been implemented.
        /// </summary>
        /// <param name="k">The shape parameter.</param>
        /// <param name="lmb">The scale parameter.</param>
        /// <returns>A MathNet 2x2 matrix which is the second derivatives calculated numerically.</returns>
        public Matrix<double> NumericalHessianLL(double k, double lmb)
        {
            Matrix<double> hessian = Matrix<double>.Build.Dense(2, 2);
            double eps = 1e-4;
            double dellmbsq = ((this.LogLikelihood(k, lmb + (2 * eps)) 
                + this.LogLikelihood(k, lmb - (2 * eps)))
                - (2 * this.LogLikelihood(k, lmb))) / (4 * (eps * eps));
            double delksq = ((this.LogLikelihood(k + (2 * eps), lmb) 
                + this.LogLikelihood(k - (2 * eps), lmb))
                - (2 * this.LogLikelihood(k, lmb))) / (4 * (eps * eps));
            double dellmbk = (this.LogLikelihood(k + eps, lmb + eps) 
                + this.LogLikelihood(k - eps, lmb - eps)
                - this.LogLikelihood(k - eps, lmb + eps) 
                - this.LogLikelihood(k + eps, lmb - eps)) / (4 * (eps * eps));

            hessian[1, 1] = dellmbsq;
            hessian[0, 0] = delksq;
            hessian[0, 1] = dellmbk;
            hessian[1, 0] = dellmbk;

            return hessian;
        }

        /// <summary>
        /// Calculates the hazard rate corresponding to the current instance of the model.
        /// </summary>
        /// <param name="t">The value (in seconds) at which the hazard rate is to be calculated</param>
        /// <returns>The hazard rate at t.</returns>
        public double HazardRate(double t)
        {
            double p = this.PDF(t);
            double s = this.Survival(t);
            return p / s;
        }

        /// <summary>
        /// Calculates the optimal reboot threshold if one exists assuming it is between 3 and 10 minutes 
        /// by simply looping over candidates.
        /// </summary>
        /// <param name="interventionCost">The amount of time in seconds it takes to recover 
        /// from a reboot operation.</param>
        /// <returns>The optimal reboot threshold as a double.</returns>
        public double GetOptimalRebootThreshold(double interventionCost)
        {
            double rebootHazardRate = 1 / interventionCost;
            double previousHazardDelta = rebootHazardRate - this.HazardRate(0.01);
            double hazardDelta = previousHazardDelta;
            double worstThreshold;
            bool doesMaximaExist = false;

            for (double step = 3; step <= 1200; step += 3)
            {
                double hazard = this.HazardRate(step);
                hazardDelta = rebootHazardRate - hazard;
                if (hazardDelta * previousHazardDelta < 0)
                {
                    if (hazardDelta > previousHazardDelta)
                    {
                        return step;
                    }
                    else
                    {
                        doesMaximaExist = true;
                        worstThreshold = step; // This was a maxima in terms of expected downtime.
                    }
                }
            }

            // We'll only get here if no minima exists
            if (doesMaximaExist)
            {
                return 600.0; // Only maxima exists, no minima. So, model says high threshold is always good.
            }
            else
            {
                // Neither minima nor maxima exists.
                if (hazardDelta > 0)
                {
                    return 0.0; // Reboot hazard rate has always been higher, so reboot as soon as possible.
                }
                else
                {
                    return 600.0; // Reboot hazard rate has always been lower than organic, so reboot as late as possible.
                }
            }
        }

        /// <summary>
        /// Executes Gradient Descent on the LogLikelihood function using the analytic gradient. To be used for 
        /// debugging purposes.
        /// </summary>
        /// <param name="parameters">The matrix of starting parameters.</param>
        /// <param name="iterations">The number of iterations.</param>
        /// <returns>The parameters which correspond to a minima for the log lihelihood.</returns>
        public Matrix<double> GradientDescent(Matrix<double> parameters, int iterations=4001)
        {
            Matrix<double> parameters1 = Matrix<double>.Build.Dense(parameters.RowCount, parameters.ColumnCount);
            Matrix<double> parameters2 = Matrix<double>.Build.Dense(parameters.RowCount, parameters.ColumnCount);
            double oldLL, newLL = 0, t = 0;
            List<double> alphas = new List<double> { 10, 2, 1, 1e-1, 1e-2, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };

            Console.WriteLine("##################Beginning Gradient Ascent###################");
            for (int i = 0; i < iterations; i++)
            {
                Matrix<double> direction = this.GradLL(parameters), temp = parameters;
                //// direction = direction / direction.L2Norm(); //// Normalizing the gradient produces an alternate approach to gradient ascent/descent.
                
                // Always take a very small step to avoid having the algorithm stuck.
                parameters2 = parameters + (1e-10 * direction); 

                oldLL = this.LogLikelihood(parameters);

                foreach (double alpha in alphas)
                { // Rudimentary line search.
                    parameters1 = parameters + (alpha * direction);
                    newLL = this.LogLikelihood(parameters1);
                    if (newLL > oldLL && !double.IsInfinity(newLL))
                    {
                        oldLL = newLL;
                        parameters2 = parameters1;
                        t = alpha;
                    }
                }

                parameters = parameters2;

                if (i%25 == 0)
                {
                    Console.Out.WriteLine("Iteration: " + i + ", objective function: " 
                        + oldLL + "\n,parameters: " + parameters.ToString() + "\n,gradient with norm " 
                        + direction.L2Norm() + " " + direction.ToString());
                }
            }

            return parameters;
        }

        /// <summary>
        /// Runs gradient descent method to find the parameters that maximize the likelihood of seeing the data.
        /// </summary>
        /// <param name="initialGuess">The initial guess for the scale and shape parameters.</param>
        /// <param name="numIter">The maximum number of iterations the routine should run.</param>
        /// <returns>A vector with the optimum parameters that maximize likelihood.</returns>
        public Vector<double> GradientDescent(Vector<double> initialGuess, int numIter = 4001)
        {
            // These are candidate step lengths. Maintained as a 
            // dictionary since we might want to see how often the step lengths are 
            // called for a given use case and fine - tune.
            Dictionary<double, int> stepLengths = new Dictionary<double, int>
            {
                { 0.00001, 0 }, { 0.0001, 0 }, { 0.001, 0 }, { 0.01, 0 }, { 0.1, 0 }, { 1.0, 0 }, { 2.0, 0 }, { 2.5, 0 }, { 3.0, 0 }, { 7.0, 0 }, { 100.0, 0 }, { 200, 0 }, { 1000.0, 0 }
            };

            Vector<double> parameters1 = initialGuess;
            Vector<double> parameters2 = initialGuess;
            Vector<double> parameters3 = initialGuess;
            Random rand = new Random();
            for (int i = 0; i < numIter; i++)
            {
                Vector<double> direction = this.GradLL(parameters1[0], parameters1[1]);
                if (Math.Abs(direction[0]) + Math.Abs(direction[1]) < 1e-5)
                {
                    this.Kappa = parameters1[0]; // Set shape parameter.
                    this.Lambda = parameters1[1]; // Set scale parameter.
                    return parameters1;
                }

                if (i % 100 > 60)
                {
                    int x = rand.Next(2);
                    direction[x] = 0.0;
                }

                double likelihood = this.LogLikelihood(parameters1[0], parameters1[1]);
                Vector<double> step = direction;

                parameters3 = parameters1 + (1e-6 * step);

                foreach (double stepLength in stepLengths.Keys)
                {
                    parameters2 = parameters1 + (stepLength * step);
                    if (Math.Min(parameters2[0], parameters2[1]) > 0)
                    {
                        double likelihood1 = this.LogLikelihood(parameters2[0], parameters2[1]);
                        if (likelihood1 > likelihood)
                        {
                            likelihood = likelihood1;
                            parameters3 = parameters2;
                        }
                    }
                }

                parameters1 = parameters3;
            }

            this.Kappa = parameters1[0]; // Set shape parameter.
            this.Lambda = parameters1[1]; // Set scale parameter.
            return parameters1;
        }

        /// <summary>
        /// Runs the Newton Rhapson method to find the parameters that maximize the likelihood of seeing the data.
        /// </summary>
        /// <param name="initialGuess">The initial guess for the scale and shape parameters.</param>
        /// <param name="numIter">The maximum number of iterations the routine should run.</param>
        /// <returns>A vector with the optimum parameters that maximize likelihood.</returns>
        public Vector<double> NewtonRaphson(Vector<double> initialGuess, int numIter = 201)
        {
            // These are candidate step lengths. Maintained as a dictionary since we might want to see how often the step lengths are 
            // called for a given use case and fine - tune.
            Dictionary<double, int> stepLengths = new Dictionary<double, int>
            {
                { 0.001, 0 }, { 0.01, 0 }, { 0.1, 0 }, { 1.0, 0 }, { 2.0, 0 }, { 2.5, 0 }, { 3.0, 0 }, { 7.0, 0 }, { 100.0, 0 }, { 200, 0 }, { 1000.0, 0 }
            };

            Vector<double> parameters1 = initialGuess;
            Vector<double> parameters2 = initialGuess;
            Vector<double> parameters3 = initialGuess;
            for (int i = 0; i < numIter; i++)
            {
                Vector<double> direction = this.GradLL(parameters1[0], parameters1[1]);
                if (Math.Abs(direction[0]) + Math.Abs(direction[1]) < 1e-5)
                {
                    this.Kappa = parameters1[0]; // Set shape parameter.
                    this.Lambda = parameters1[1]; // Set scale parameter.
                    return parameters1;
                }

                double likelihood = this.LogLikelihood(parameters1[0], parameters1[1]);
                Matrix<double> hess = this.NumericalHessianLL(parameters1[0], parameters1[1]);
                Vector<double> step = hess.Solve(direction);

                parameters3 = parameters1 - (1e-6 * step);

                foreach (double stepLength in stepLengths.Keys)
                {
                    parameters2 = parameters1 - (stepLength * step);
                    if (Math.Min(parameters2[0], parameters2[1]) > 0)
                    {
                        double likelihood1 = this.LogLikelihood(parameters2[0], parameters2[1]);
                        if (likelihood1 > likelihood)
                        {
                            likelihood = likelihood1;
                            parameters3 = parameters2;
                        }
                    }
                }

                parameters1 = parameters3;
            }

            this.Kappa = parameters1[0]; // Set shape parameter.
            this.Lambda = parameters1[1]; // Set scale parameter.
            return parameters1;
        }

        /// <summary>
        /// Generates dummy feature matrices when not available.
        /// </summary>
        /// <param name="numFeatures">The number of features on which the matrices will be based.</param>
        /// <returns>A matrix of features.</returns>
        protected Tuple<Matrix<double>, Matrix<double>> PopulateFeatureMatrices(int numFeatures)
        {
            Matrix<double> fSamples = Matrix<double>.Build.Dense(this.OrganicRecoveryDurations.Count, numFeatures);
            Matrix<double> fCensored = Matrix<double>.Build.Dense(this.OrganicRecoveryDurations.Count, numFeatures);

            for (int i = 0; i < fSamples.RowCount; i++)
            {
                fSamples[i, 0] = 1.0;
            }

            for (int i = 0; i < fCensored.RowCount; i++)
            {
                fCensored[i, 0] = 1.0;
            }

            return new Tuple<Matrix<double>, Matrix<double>>(fSamples, fCensored);
        }
    }
}
