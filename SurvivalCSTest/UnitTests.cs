using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Linq;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SurvivalCS;
using SurvivalCS.Distributions;
using SurvivalCS.Functions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SurvivalCSTest
{
    [TestClass]
    public class UnitTests
    {
        [TestMethod]
        public void TestLomaxFitting()
        {
            Lomax lx = new Lomax(1.01, 100);
            var lxSamples = lx.GenerateSample(10000);

            Lomax lx2 = new Lomax(lxSamples);
            lx2.FitParameters();
            Console.WriteLine(lx2.Kappa);

            //The un-biased estimator of the Lomax distribution
            //has a high variance. So, we have to set the bar low.
            Assert.IsTrue(Math.Abs(lx2.Kappa - 1.01) < 1e-1);
        }

        [TestMethod]
        public void TestLogLogisticFitting()
        {
            LogLogistic lg = new LogLogistic(1.2, 300);
            var lgSamples = lg.GenerateSample(3000);
            var lgCensored = new List<double> {0.1, 0.1};

            double[] arr = new double[] { 0.5, 12.0 };
            Vector<double> init = Vector<double>.Build.DenseOfArray(arr);

            LogLogistic lg2 = new LogLogistic(lgSamples, lgCensored);
            lg2.GradientDescent(init);

            Console.WriteLine(lg2.Kappa);

            //The un-biased estimator of the LogLogistic distribution
            //has a high variance. So, we have to set the bar low.
            Assert.IsTrue(Math.Abs(lg2.Kappa - 1.2) < 1e-1);
        }

        [TestMethod]
        public void TestMixedLogLogisticFitting()
        {
            LogLogistic lg = new LogLogistic(1.2,300.0);
            List<double> llSamples = lg.GenerateSample(5000);
            LogLogistic lg1 = new LogLogistic(0.7, 80.0);
            List<double> llSamples1 = lg1.GenerateSample(5000);
            DataGen.AppendToLst(llSamples, llSamples1);
            var ti = llSamples;
            List<double> xi = new List<double> {.1};
            LogLogistic llModel = new LogLogistic(ti, xi);
            double[] arr = new double[] { 0.8, 150.0 };
            Vector<double> init = Vector<double>.Build.DenseOfArray(arr);
            llModel.GradientDescent(init);
            Debug.WriteLine("Model completed!");
            Assert.IsTrue(Math.Abs(llModel.Kappa-0.83)<1e-1);
        }

        [TestMethod]
        public void TestLogLogisticWFeatures()
        {
            DataGen dg = new DataGen();
            dg.GenLogLogisticWFeatures();

            double shapeMax = 5.0;
            double scaleMax = 500.0;

            double[] arr = new double[] { 1.0, 150.0 };
            Vector<double> init = Vector<double>.Build.DenseOfArray(arr);
            LogLogistic modelLogLogistic = new LogLogistic(dg.organicRecoveryDurations, dg.inorganicRecoverydurations);
            modelLogLogistic.GradientDescent(init);
            
            Console.WriteLine("LL without features is " + 
                modelLogLogistic.LogLikelihood(modelLogLogistic.Kappa, modelLogLogistic.Lambda) + 
                " with Kappa " + modelLogLogistic.Kappa + " and Lambda " + modelLogLogistic.Lambda);

            double[,] warr = new double[2, dg.fCensored.ColumnCount];
            warr[0, 0] = Sigmoid.InverseSigmoid(modelLogLogistic.Kappa, shapeMax);
            warr[1, 0] = Sigmoid.InverseSigmoid(modelLogLogistic.Lambda, scaleMax);
            Matrix<double>  w = Matrix<double>.Build.DenseOfArray(warr);

            LogLogistic modelLogLogisticFeatured = new LogLogistic(dg.organicRecoveryDurations, 
                dg.inorganicRecoverydurations,
                dg.fSamples, dg.fCensored);
            modelLogLogisticFeatured.ShapeUpperBound = shapeMax;
            modelLogLogisticFeatured.ScaleUpperBound = scaleMax;
            Matrix<double> logLogisticParameters = modelLogLogisticFeatured.GradientDescent(w, 2001);
            Vector<double> frstSample = Vector<double>.Build.DenseOfArray(new double[] {1.0, 2.0, 3.0});
            Vector<double> scndSample = Vector<double>.Build.DenseOfArray(new double[] {1.0, 4.0, 2.0});

            Vector<double> res = logLogisticParameters.Multiply(frstSample);
            var alpha_shape = res[0];
            var shape = Sigmoid.Transform(alpha_shape, shapeMax);
            var alpha_scale = res[1];
            var scale = Sigmoid.Transform(alpha_scale, scaleMax);

            res = logLogisticParameters.Multiply(scndSample);
            alpha_shape = res[0];
            shape = Sigmoid.Transform(alpha_shape, shapeMax);
            alpha_scale = res[1];
            scale = Sigmoid.Transform(alpha_scale, scaleMax);
            Assert.IsTrue(Math.Abs(scale-80.0)<2.0);
        }
    }
}
