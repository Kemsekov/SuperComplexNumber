using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.Differentiation;

/// <summary>
/// Taylor series that could be used to generalize real values functions to complex and supercomplex plane
/// </summary>
public static class TaylorSeries
{
    
    /// <param name="func">Function to take derivative</param>
    /// <param name="approximationPoint">Where is approximation point</param>
    /// <param name="epsilon">Finite difference epsilon approximation</param>
    /// <returns>
    /// Function that is capable of accumulating derivatives values to speed up 
    /// computation of higher order derivative based on previous values calculated
    /// </returns>
    public static Func<int,double,double> BuildDerivativesAccumulator(Func<double, double> func, double approximationPoint, double epsilon){
        Dictionary<int, Dictionary<double, double>> derivatives = new();

        var precision = (int)Math.Abs(Math.Log10(epsilon));
        //derivatives[n][x] = n'th derivative of func(x)
        double readDer(int n, double x)
        {
            if(!derivatives.ContainsKey(n))
                derivatives[n] = new();
            x = double.Round(x, precision);
            if(n<=0) return func(x);
            var der = derivatives[n];
            if (der.ContainsKey(x))
                return der[x];
            else{
                var d = (readDer(n-1,x+epsilon)-readDer(n-1,x-epsilon))*0.5/epsilon;
                der[x]=d;
                return d;
            }
        }
        return readDer;
    }
    public static double[] BuildSeriesCustom(Func<double, double> func, double approximationPoint, int order, double epsilon)
    {
        var der = BuildDerivativesAccumulator(func,approximationPoint,epsilon);

        var A =
            Enumerable.Range(1, order)
            .Select(r => der(r, approximationPoint))
            .Prepend(func(approximationPoint))
            .Select((val, i) => val / SpecialFunctions.Factorial(i))
            .ToArray();
        return A;


    }
    /// <summary>
    /// Etalon series builder. The slowest, but most precise.
    /// </summary>
    public static double[] BuildSeriesEtalon(Func<double, double> func, double approximationPoint, int order)
    {
        // coefficients.GetCoefficientsForAllOrders
        var A =
            Enumerable.Range(1, order)
            .Select(r => Differentiate.Derivative(func, approximationPoint, r))
            .Prepend(func(approximationPoint))
            .Select((val, i) => val / SpecialFunctions.Factorial(i))
            .ToArray();
        return A;
    }
    /// <summary>
    /// Works a lot faster than etalon, but diverges from it at around 30 terms -> cuz floating point numbers works bad on pc =(
    /// </summary>
    public static double[] BuildSeries(Func<double, double> func, double approximationPoint, int order)
    {
        // order++;
        var points = order + (order.IsEven() ? 1 : 2);
        var coefficients = new FiniteDifferenceCoefficients(points);
        var dif = new MyNumericalDerivative(points, points / 2);
        var (pointsEvaluation, _, stepSize) = dif.ComputeStaff(func, approximationPoint, order, func(approximationPoint));

        // coefficients.GetCoefficientsForAllOrders
        var A =
            Enumerable.Range(1, order)
            .Select(r => dif.EvaluateDerivative(pointsEvaluation, r, stepSize))
            // .Select(r => new NumericalDerivative(r + (r.IsEven() ? 1 : 2),(r + (r.IsEven() ? 1 : 2))/2).EvaluateDerivative(pointsEvaluation, r, stepSize))
            // .Select(r => dif.EvaluateDerivative(func,approximationPoint,r))
            .Prepend(func(approximationPoint))
            .Select((val, i) => val / SpecialFunctions.Factorial(i))
            .ToArray();
        return A;
    }

    public static double EvaluateSeries(double[] coefficients, double approximationPoint, double x)
    {
        return Series.Evaluate(coefficients.Select((c, i) => Math.Pow(x - approximationPoint, i) * c));
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