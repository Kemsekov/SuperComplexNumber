using System;
using System.Linq;
#pragma warning disable
namespace MathNet.Numerics.Differentiation;

// I had to do whole copy of MathNet.Numerics.Differentiation.NumericalDerivative because
// I needed to optimize my taylor series computation

//
// Summary:
//     Class to evaluate the numerical derivative of a function using finite difference
//     approximations. Variable point and center methods can be initialized MathNet.Numerics.Differentiation.FiniteDifferenceCoefficients.
//     This class can also be used to return function handles (delegates) for a fixed
//     derivative order and variable. It is possible to evaluate the derivative and
//     partial derivative of univariate and multivariate functions respectively.
public class MyNumericalDerivative
{
    private readonly int _points;

    private int _center;

    private double _stepSize = Math.Pow(2.0, -10.0);

    private double _epsilon = Precision.PositiveMachineEpsilon;

    private double _baseStepSize = Math.Pow(2.0, -26.0);

    private readonly FiniteDifferenceCoefficients _coefficients;

    //
    // Summary:
    //     Sets and gets the finite difference step size. This value is for each function
    //     evaluation if relative step size types are used. If the base step size used in
    //     scaling is desired, see MathNet.Numerics.Differentiation.NumericalDerivative.Epsilon.
    //
    //
    // Remarks:
    //     Setting then getting the StepSize may return a different value. This is not unusual
    //     since a user-defined step size is converted to a base-2 representable number
    //     to improve finite difference accuracy.
    public double StepSize
    {
        get
        {
            return _stepSize;
        }
        set
        {
            double a = Math.Log(Math.Abs(value)) / Math.Log(2.0);
            _stepSize = Math.Pow(2.0, Math.Round(a));
        }
    }

    //
    // Summary:
    //     Sets and gets the base finite difference step size. This assigned value to this
    //     parameter is only used if MathNet.Numerics.Differentiation.NumericalDerivative.StepType
    //     is set to RelativeX. However, if the StepType is Relative, it will contain the
    //     base step size computed from MathNet.Numerics.Differentiation.NumericalDerivative.Epsilon
    //     based on the finite difference order.
    public double BaseStepSize
    {
        get
        {
            return _baseStepSize;
        }
        set
        {
            double a = Math.Log(Math.Abs(value)) / Math.Log(2.0);
            _baseStepSize = Math.Pow(2.0, Math.Round(a));
        }
    }

    //
    // Summary:
    //     Sets and gets the base finite difference step size. This parameter is only used
    //     if MathNet.Numerics.Differentiation.NumericalDerivative.StepType is set to Relative.
    //     By default this is set to machine epsilon, from which MathNet.Numerics.Differentiation.NumericalDerivative.BaseStepSize
    //     is computed.
    public double Epsilon
    {
        get
        {
            return _epsilon;
        }
        set
        {
            double a = Math.Log(Math.Abs(value)) / Math.Log(2.0);
            _epsilon = Math.Pow(2.0, Math.Round(a));
        }
    }

    //
    // Summary:
    //     Sets and gets the location of the center point for the finite difference derivative.
    public int Center
    {
        get
        {
            return _center;
        }
        set
        {
            if (value >= _points || value < 0)
            {
                throw new ArgumentOutOfRangeException("value", "Center must lie between 0 and points -1");
            }

            _center = value;
        }
    }

    //
    // Summary:
    //     Number of times a function is evaluated for numerical derivatives.
    public int Evaluations { get; private set; }

    //
    // Summary:
    //     Type of step size for computing finite differences. If set to absolute, dx =
    //     h. If set to relative, dx = (1+abs(x))*h^(2/(order+1)). This provides accurate
    //     results when h is approximately equal to the square-root of machine accuracy,
    //     epsilon.
    public StepType StepType { get; set; } = StepType.Relative;


    //
    // Summary:
    //     Initializes a NumericalDerivative class with the default 3 point center difference
    //     method.
    public MyNumericalDerivative()
        : this(3, 1)
    {
    }

    //
    // Summary:
    //     Initialized a NumericalDerivative class.
    //
    // Parameters:
    //   points:
    //     Number of points for finite difference derivatives.
    //
    //   center:
    //     Location of the center with respect to other points. Value ranges from zero to
    //     points-1.
    public MyNumericalDerivative(int points, int center)
    {
        if (points < 2)
        {
            throw new ArgumentOutOfRangeException("points", "Points must be two or greater.");
        }

        _center = center;
        _points = points;
        Center = center;
        _coefficients = new FiniteDifferenceCoefficients(points);
    }

    //
    // Summary:
    //     Evaluates the derivative of equidistant points using the finite difference method.
    //
    //
    // Parameters:
    //   points:
    //     Vector of points StepSize apart.
    //
    //   order:
    //     Derivative order.
    //
    //   stepSize:
    //     Finite difference step size.
    //
    // Returns:
    //     Derivative of points of the specified order.
    public double EvaluateDerivative(double[] points, int order, double stepSize)
    {
        if (points == null)
        {
            throw new ArgumentNullException("points");
        }

        if (order >= _points || order < 0)
        {
            throw new ArgumentOutOfRangeException("order", "Order must be between zero and points-1.");
        }

        return _coefficients.GetCoefficients(Center, order).Select((double t, int i) => t * points[i]).Sum() / Math.Pow(stepSize, order);
    }

    //
    // Summary:
    //     Evaluates the derivative of a scalar univariate function.
    //
    // Parameters:
    //   f:
    //     Function handle.
    //
    //   x:
    //     Point at which to compute the derivative.
    //
    //   order:
    //     Derivative order.
    //
    //   currentValue:
    //     Current function value at center.
    //
    // Returns:
    //     Function derivative at x of the specified order.
    //
    // Remarks:
    //     Supplying the optional argument currentValue will reduce the number of function
    //     evaluations required to calculate the finite difference derivative.
    public double EvaluateDerivative(Func<double, double> f, double x, int order, double? currentValue = null)
    {
        var (array,_,num) = ComputeStaff(f,x,order,currentValue);
        return EvaluateDerivative(array, order, num);
    }

    public (double[] evaluations, int order, double num) ComputeStaff(Func<double, double> f, double x, int order, double? currentValue = null){
        double[] coefficients = _coefficients.GetCoefficients(Center, order);
        double num = CalculateStepSize(_points, x, order);
        double[] array = new double[_points];
        for (int i = 0; i < _points; i++)
        {
            if (i == Center && currentValue.HasValue)
            {
                array[i] = currentValue.Value;
            }
            else if (coefficients[i] != 0.0)
            {
                array[i] = f(x + (double)(i - Center) * num);
                Evaluations++;
            }
        }
        return (array,order,num);
    }

    //
    // Summary:
    //     Creates a function handle for the derivative of a scalar univariate function.
    //
    //
    // Parameters:
    //   f:
    //     Input function handle.
    //
    //   order:
    //     Derivative order.
    //
    // Returns:
    //     Function handle that evaluates the derivative of input function at a fixed order.
    public Func<double, double> CreateDerivativeFunctionHandle(Func<double, double> f, int order)
    {
        return (double x) => EvaluateDerivative(f, x, order);
    }

    //
    // Summary:
    //     Evaluates the partial derivative of a multivariate function.
    //
    // Parameters:
    //   f:
    //     Multivariate function handle.
    //
    //   x:
    //     Vector at which to evaluate the derivative.
    //
    //   parameterIndex:
    //     Index of independent variable for partial derivative.
    //
    //   order:
    //     Derivative order.
    //
    //   currentValue:
    //     Current function value at center.
    //
    // Returns:
    //     Function partial derivative at x of the specified order.
    public double EvaluatePartialDerivative(Func<double[], double> f, double[] x, int parameterIndex, int order, double? currentValue = null)
    {
        double num = x[parameterIndex];
        double[] coefficients = _coefficients.GetCoefficients(Center, order);
        double num2 = CalculateStepSize(_points, x[parameterIndex], order);
        double[] array = new double[_points];
        for (int i = 0; i < _points; i++)
        {
            if (i == Center && currentValue.HasValue)
            {
                array[i] = currentValue.Value;
            }
            else if (coefficients[i] != 0.0)
            {
                x[parameterIndex] = num + (double)(i - Center) * num2;
                array[i] = f(x);
                Evaluations++;
            }
        }

        x[parameterIndex] = num;
        return EvaluateDerivative(array, order, num2);
    }
    //
    // Summary:
    //     Evaluates the partial derivatives of a multivariate function array.
    //
    // Parameters:
    //   f:
    //     Multivariate vector function array handle.
    //
    //   x:
    //     Vector at which to evaluate the derivatives.
    //
    //   parameterIndex:
    //     Index of independent variable for partial derivative.
    //
    //   order:
    //     Derivative order.
    //
    //   currentValue:
    //     Current function value at center.
    //
    // Returns:
    //     Vector of functions partial derivatives at x of the specified order.
    //
    // Remarks:
    //     This function assumes the input vector x is of the correct length for f.
    public double[] EvaluatePartialDerivative(Func<double[], double>[] f, double[] x, int parameterIndex, int order, double?[] currentValue = null)
    {
        double[] array = new double[f.Length];
        for (int i = 0; i < f.Length; i++)
        {
            if (currentValue != null && currentValue[i].HasValue)
            {
                array[i] = EvaluatePartialDerivative(f[i], x, parameterIndex, order, currentValue[i].Value);
            }
            else
            {
                array[i] = EvaluatePartialDerivative(f[i], x, parameterIndex, order);
            }
        }

        return array;
    }

    //
    // Summary:
    //     Creates a function handle for the partial derivative of a multivariate function.
    //
    //
    // Parameters:
    //   f:
    //     Input function handle.
    //
    //   parameterIndex:
    //     Index of the independent variable for partial derivative.
    //
    //   order:
    //     Derivative order.
    //
    // Returns:
    //     Function handle that evaluates partial derivative of input function at a fixed
    //     order.
    public Func<double[], double> CreatePartialDerivativeFunctionHandle(Func<double[], double> f, int parameterIndex, int order)
    {
        return (double[] x) => EvaluatePartialDerivative(f, x, parameterIndex, order);
    }

    //
    // Summary:
    //     Creates a function handle for the partial derivative of a vector multivariate
    //     function.
    //
    // Parameters:
    //   f:
    //     Input function handle.
    //
    //   parameterIndex:
    //     Index of the independent variable for partial derivative.
    //
    //   order:
    //     Derivative order.
    //
    // Returns:
    //     Function handle that evaluates partial derivative of input function at fixed
    //     order.
    public Func<double[], double[]> CreatePartialDerivativeFunctionHandle(Func<double[], double>[] f, int parameterIndex, int order)
    {
        return (double[] x) => EvaluatePartialDerivative(f, x, parameterIndex, order);
    }

    //
    // Summary:
    //     Evaluates the mixed partial derivative of variable order for multivariate functions.
    //
    //
    // Parameters:
    //   f:
    //     Multivariate function handle.
    //
    //   x:
    //     Points at which to evaluate the derivative.
    //
    //   parameterIndex:
    //     Vector of indices for the independent variables at descending derivative orders.
    //
    //
    //   order:
    //     Highest order of differentiation.
    //
    //   currentValue:
    //     Current function value at center.
    //
    // Returns:
    //     Function mixed partial derivative at x of the specified order.
    //
    // Remarks:
    //     This function recursively uses MathNet.Numerics.Differentiation.NumericalDerivative.EvaluatePartialDerivative(System.Func{System.Double[],System.Double},System.Double[],System.Int32,System.Int32,System.Nullable{System.Double})
    //     to evaluate mixed partial derivative. Therefore, it is more efficient to call
    //     MathNet.Numerics.Differentiation.NumericalDerivative.EvaluatePartialDerivative(System.Func{System.Double[],System.Double},System.Double[],System.Int32,System.Int32,System.Nullable{System.Double})
    //     for higher order derivatives of a single independent variable.
    public double EvaluateMixedPartialDerivative(Func<double[], double> f, double[] x, int[] parameterIndex, int order, double? currentValue = null)
    {
        if (parameterIndex.Length != order)
        {
            throw new ArgumentOutOfRangeException("parameterIndex", "The number of parameters must match derivative order.");
        }

        if (order == 1)
        {
            return EvaluatePartialDerivative(f, x, parameterIndex[0], order, currentValue);
        }

        int num = order - 1;
        int[] array = new int[num];
        Array.Copy(parameterIndex, 0, array, 0, num);
        double[] array2 = new double[_points];
        int num2 = parameterIndex[order - 1];
        double num3 = CalculateStepSize(_points, x[num2], order);
        double num4 = x[num2];
        for (int i = 0; i < _points; i++)
        {
            x[num2] = num4 + (double)(i - Center) * num3;
            array2[i] = EvaluateMixedPartialDerivative(f, x, array, num);
        }

        x[num2] = num4;
        return EvaluateDerivative(array2, 1, num3);
    }

    //
    // Summary:
    //     Evaluates the mixed partial derivative of variable order for multivariate function
    //     arrays.
    //
    // Parameters:
    //   f:
    //     Multivariate function array handle.
    //
    //   x:
    //     Vector at which to evaluate the derivative.
    //
    //   parameterIndex:
    //     Vector of indices for the independent variables at descending derivative orders.
    //
    //
    //   order:
    //     Highest order of differentiation.
    //
    //   currentValue:
    //     Current function value at center.
    //
    // Returns:
    //     Function mixed partial derivatives at x of the specified order.
    //
    // Remarks:
    //     This function recursively uses MathNet.Numerics.Differentiation.NumericalDerivative.EvaluatePartialDerivative(System.Func{System.Double[],System.Double}[],System.Double[],System.Int32,System.Int32,System.Nullable{System.Double}[])
    //     to evaluate mixed partial derivative. Therefore, it is more efficient to call
    //     MathNet.Numerics.Differentiation.NumericalDerivative.EvaluatePartialDerivative(System.Func{System.Double[],System.Double}[],System.Double[],System.Int32,System.Int32,System.Nullable{System.Double}[])
    //     for higher order derivatives of a single independent variable.
    public double[] EvaluateMixedPartialDerivative(Func<double[], double>[] f, double[] x, int[] parameterIndex, int order, double?[] currentValue = null)
    {
        double[] array = new double[f.Length];
        for (int i = 0; i < f.Length; i++)
        {
            if (currentValue != null && currentValue[i].HasValue)
            {
                array[i] = EvaluateMixedPartialDerivative(f[i], x, parameterIndex, order, currentValue[i].Value);
            }
            else
            {
                array[i] = EvaluateMixedPartialDerivative(f[i], x, parameterIndex, order);
            }
        }

        return array;
    }

    //
    // Summary:
    //     Creates a function handle for the mixed partial derivative of a multivariate
    //     function.
    //
    // Parameters:
    //   f:
    //     Input function handle.
    //
    //   parameterIndex:
    //     Vector of indices for the independent variables at descending derivative orders.
    //
    //
    //   order:
    //     Highest derivative order.
    //
    // Returns:
    //     Function handle that evaluates the fixed mixed partial derivative of input function
    //     at fixed order.
    public Func<double[], double> CreateMixedPartialDerivativeFunctionHandle(Func<double[], double> f, int[] parameterIndex, int order)
    {
        return (double[] x) => EvaluateMixedPartialDerivative(f, x, parameterIndex, order);
    }

    //
    // Summary:
    //     Creates a function handle for the mixed partial derivative of a multivariate
    //     vector function.
    //
    // Parameters:
    //   f:
    //     Input vector function handle.
    //
    //   parameterIndex:
    //     Vector of indices for the independent variables at descending derivative orders.
    //
    //
    //   order:
    //     Highest derivative order.
    //
    // Returns:
    //     Function handle that evaluates the fixed mixed partial derivative of input function
    //     at fixed order.
    public Func<double[], double[]> CreateMixedPartialDerivativeFunctionHandle(Func<double[], double>[] f, int[] parameterIndex, int order)
    {
        return (double[] x) => EvaluateMixedPartialDerivative(f, x, parameterIndex, order);
    }

    //
    // Summary:
    //     Resets the evaluation counter.
    public void ResetEvaluations()
    {
        Evaluations = 0;
    }

    private double CalculateStepSize(int points, double x, double order)
    {
        if (StepType == StepType.RelativeX)
        {
            StepSize = BaseStepSize * (1.0 + Math.Abs(x));
        }
        else if (StepType == StepType.Relative)
        {
            double num = (double)points - order;
            BaseStepSize = Math.Pow(Epsilon, 1.0 / (num + order));
            StepSize = BaseStepSize * (1.0 + Math.Abs(x));
        }

        return StepSize;
    }
}