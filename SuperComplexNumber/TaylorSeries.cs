using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.Differentiation;

/// <summary>
/// Taylor series that could be used to generalize real values functions to complex and supercomplex plane
/// </summary>
public static class TaylorSeries
{
    public static double[] BuildSeries(Func<double, double> func, double approximationPoint, int order)
    {
        var points = order + (order.IsEven() ? 1 : 2);
        var coefficients = new FiniteDifferenceCoefficients(points);
        var dif = new MyNumericalDerivative(points, points / 2);
        var (pointsEvaluation,_,stepSize) = dif.ComputeStaff(func,approximationPoint,order,func(approximationPoint));
        
        // coefficients.GetCoefficientsForAllOrders
        var A =
            Enumerable.Range(1, order)
            .Select(r => dif.EvaluateDerivative(pointsEvaluation,r,stepSize))
            .Prepend(func(approximationPoint))
            .Select((val, i) => val / SpecialFunctions.Factorial(i))
            .ToArray();
        return A;
    }

    public static double EvaluateSeries(double[] coefficients, double approximationPoint, double x)
    {
        return Series.Evaluate(coefficients.Select((c, i) => Math.Pow(x-approximationPoint, i) * c));
    }

    public static SuperComplex EvaluateSeriesOnSuperComplex(double[] coefficients, double approximationPoint, SuperComplex x)
    {
        return x.ApplyComplex(c => EvaluateSeriesOnComplex(coefficients, approximationPoint, c));
    }

    public static Complex EvaluateSeriesOnComplex(double[] coefficients, double approximationPoint, Complex x)
    {
        x -= approximationPoint;
        var power = new Complex(1, 0);
        var sum = new Complex(0, 0);
        foreach (var c in coefficients)
        {
            sum += c * power;
            power *= x;
        }
        return sum;
    }

}