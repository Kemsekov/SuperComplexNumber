
global using CMatrix = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix;
global using DMatrix = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
global using CVector = MathNet.Numerics.LinearAlgebra.Complex.DenseVector;
global using DVector = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq.Expressions;
using System.Security.AccessControl;
using System.Net.Quic;
public struct SuperComplex
{
    public double[] Coefficients { get; }
    double this[int index]
    {
        get => Coefficients[index];
        set => Coefficients[index] = value;
    }
    public SuperComplex()
    {
        Coefficients = new double[4];
    }
    public SuperComplex(double[] coefficients)
    {
        Coefficients = coefficients.Take(4).ToArray();
    }
    public SuperComplex(Complex[] coefficients)
    {
        Coefficients = castToCoefficients(coefficients);
    }
    public SuperComplex(SuperComplex[] coefficients)
    {
        Coefficients = castToCoefficients(coefficients);
    }
    public SuperComplex(Matrix<double> value)
    {
        Coefficients = castToCoefficients(value);
    }
    public SuperComplex(Matrix<Complex> value)
    {
        Coefficients = castToCoefficients(value);
    }
    public Matrix<double> ToMatrix()
    {
        return real * Coefficients[0] + imaginary * Coefficients[1] + theta * Coefficients[2] + mu * Coefficients[3];
    }
    public Complex ToComplex() => new(Coefficients[0], Coefficients[1]);
    public static implicit operator Matrix<double>(SuperComplex val)
    {
        return val.ToMatrix();
    }
    public static implicit operator SuperComplex(double[] coef)
    {
        return new(coef);
    }
    public static implicit operator SuperComplex(SuperComplex[] coef)
    {
        return new(coef);
    }
    public static implicit operator SuperComplex(double real)
    {
        return new(new[] { real, 0, 0, 0 });
    }
    public static implicit operator SuperComplex(Complex complex)
    {
        return new(new[] { complex.Real, complex.Imaginary, 0, 0 });
    }
    public static implicit operator SuperComplex(Complex[] coef)
    {
        return new(coef);
    }
    public static implicit operator SuperComplex(Matrix<double> val)
    {
        return new(val);
    }
    public static implicit operator SuperComplex(Matrix<Complex> val)
    {
        return new(val);
    }
    public static implicit operator double[](SuperComplex num)
    {
        return num.Coefficients;
    }
    public static SuperComplex operator +(SuperComplex s1, SuperComplex s2)
    {
        return new SuperComplex(
            s1.Coefficients.Zip(s2.Coefficients)
            .Select(t => t.First + t.Second)
            .ToArray()
        );
    }
    public static SuperComplex operator -(SuperComplex s1, SuperComplex s2)
    {
        return new SuperComplex(
            s1.Coefficients.Zip(s2.Coefficients)
            .Select(t => t.First - t.Second)
            .ToArray()
        );
    }
    
    public static bool operator ==(SuperComplex s1,SuperComplex s2){
        return s1.Coefficients.Zip(s2.Coefficients).Sum(f=>Math.Abs(1-f.First/f.Second))<1e-12;
    }
    public static bool operator !=(SuperComplex s1,SuperComplex s2){
        return s1.Coefficients.Zip(s2.Coefficients).Sum(f=>Math.Abs(1-f.First/f.Second))>=1e-12;
    }
    public static bool operator>(SuperComplex s1,SuperComplex s2){
        return s1.Magnitude>s2.Magnitude;
    }
    public static bool operator<(SuperComplex s1,SuperComplex s2){
        return s1.Magnitude<s2.Magnitude;
    }
    public static SuperComplex operator -(SuperComplex s)
    {
        return new SuperComplex(s.Coefficients.Select(x => -x).ToArray());
    }
    public static SuperComplex operator *(SuperComplex s1, SuperComplex s2)
    {
        var result = s1.ToMatrix() * s2.ToMatrix();
        return new(result);
    }
    //TODO add function transformation and check if this operation is consistent with function 1/Complex
    public static SuperComplex operator /(SuperComplex s1, SuperComplex s2)
    {
        if (s1.Magnitude == 0) return new(new[] { 0.0, 0, 0, 0 });
        var result = s1.ToMatrix() * s2.ToMatrix().Inverse();
        return new(result);
    }
    public double Magnitude => Math.Sqrt(Coefficients.Sum(x => x * x));

    /// <summary>
    /// Applies function by analytic continuation to super complex number from complex plane.<br/>
    /// Here is a catch: you may notice that <see cref="Complex.Sqrt"/> results when squared does not
    /// yield same value that you used to take sqrt, when magnitude of value is > 1<br/>
    /// This is not mistake, this method indeed generalizes complex functions
    /// to super complex plane, but it does so ONLY IN CORRESPONDING SERIES CONVERGENCE RADIUS. <br/>
    /// That is, for sqrt magnitude of input value must be smaller than 1, so when you try to take sqrt of
    /// anything that is bigger than 1 by magnitude, it will give you wrong results. <br/>
    /// That's crazy how does math know about series convergence even here, although there is no series used at all.
    /// </summary>
    public SuperComplex ApplyComplex(Func<Complex, Complex> func)
    {
        // see https://en.wikipedia.org/wiki/Analytic_function_of_a_matrix
        var mat = ToMatrix();
        var A = CMatrix.Create(2,2,(i,j)=>mat[i,j]);
        var evd = mat.Evd();

        var lp = evd.EigenValues.MaxBy(x=>x.Real);
        var lm = evd.EigenValues.MinBy(x=>x.Real);

        var trace = mat.Trace();
        var sum = func(lp)+func(lm);
        var diff = func(lp)-func(lm);

        var identity = CMatrix.CreateIdentity(2);
        var sqrt = Complex.Sqrt(trace*trace/4-mat.Determinant());
        //I have no idea why setting it to NaN works, but it does
        if(sqrt==0) sqrt=Complex.NaN;
        var part1 = sum*identity/2;
        var part2 = (A-identity*trace/2)/sqrt*diff/2;
        var res = part1+part2;
        return res;
    }

    double[] castToCoefficients(Complex[] coefficients)
    {
        Matrix<double> total = DenseMatrix.Create(2, 2, 0);
        total += coefficients[0].Real * real + coefficients[0].Imaginary * imaginary;
        total += coefficients[1].Real * imaginary - coefficients[1].Imaginary * real;
        total += (coefficients[2].Real * real - coefficients[2].Imaginary * imaginary) * theta;
        total += (coefficients[3].Real * real - coefficients[3].Imaginary * imaginary) * mu;
        return castToCoefficients(total);
    }
    double[] castToCoefficients(SuperComplex[] coefficients)
    {
        SuperComplex total = DenseMatrix.Create(2, 2, 0);
        var part = new SuperComplex[] { real, imaginary, theta, mu };
        foreach (var (coef, complexPart) in coefficients.Zip(part))
            total += coef * complexPart;
        return total.Coefficients;
    }
    
    double[] castToCoefficients(Matrix<Complex> val)
    {
        var coefficentsMatrix = CMatrix.Create(4,4,0);
        coefficentsMatrix.SetColumn(0,real.ToRowMajorArray().Select(x=>new Complex(x,0)).ToArray());
        coefficentsMatrix.SetColumn(1,imaginary.ToRowMajorArray().Select(x=>new Complex(x,0)).ToArray());
        coefficentsMatrix.SetColumn(2,theta.ToRowMajorArray().Select(x=>new Complex(x,0)).ToArray());
        coefficentsMatrix.SetColumn(3,mu.ToRowMajorArray().Select(x=>new Complex(x,0)).ToArray());

        var valVec = CVector.OfArray(val.ToRowMajorArray().Cast<Complex>().ToArray());
        return castToCoefficients(coefficentsMatrix.Solve(valVec).ToArray());
    }
    double[] castToCoefficients(Matrix<double> val)
    {
        var coefficentsMatrix = DMatrix.Create(4,4,0);
        coefficentsMatrix.SetColumn(0,real.ToRowMajorArray());
        coefficentsMatrix.SetColumn(1,imaginary.ToRowMajorArray());
        coefficentsMatrix.SetColumn(2,theta.ToRowMajorArray());
        coefficentsMatrix.SetColumn(3,mu.ToRowMajorArray());

        var valVec = DVector.OfArray(val.ToRowMajorArray());
        return coefficentsMatrix.Solve(valVec).ToArray();
    }
    //these values should not be changed, so they are private
    static Matrix<double> imaginary = Imaginary;
    static Matrix<double> real = Real;
    static Matrix<double> theta = Theta;
    static Matrix<double> mu = Mu;
    public override string ToString()
    {
        return ToString();
    }
    public string ToString(string? numbersFormat = null)
    {
        char sign(double val) => val >= 0 ? '+' : '-';
        return $"({Coefficients[0].ToString(numbersFormat)} {sign(Coefficients[1])} {Math.Abs(Coefficients[1]).ToString(numbersFormat)}i {sign(Coefficients[2])} {Math.Abs(Coefficients[2]).ToString(numbersFormat)}θ {sign(Coefficients[3])} {Math.Abs(Coefficients[3]).ToString(numbersFormat)}μ)";
    }
    /// <summary>
    /// i^2=-1
    /// </summary>
    public static Matrix<double> Imaginary
    {
        get
        {
            var i = DenseMatrix.Create(2, 2, 0); //i*i=-1
            i[0, 1] = -1;
            i[1, 0] = 1;
            return i;
        }
    }
    /// <summary>
    /// Real 1
    /// </summary>
    public static Matrix<double> Real
    {
        get
        {
            var one = DenseMatrix.Create(2, 2, 0);
            one[0, 0] = 1;
            one[1, 1] = 1;
            return one;
        }
    }
    /// <summary>
    /// theta^2=0
    /// </summary>
    public static Matrix<double> Theta
    {
        get
        {
            var theta = DenseMatrix.Create(2, 2, 0); //theta*theta=0
            theta[0, 0] = 1;
            theta[0, 1] = 1;
            theta[1, 0] = -1;
            theta[1, 1] = -1;
            return theta;
        }
    }
    /// <summary>
    /// Sqrt of one
    /// </summary>
    public static Matrix<double> Mu
    {
        get
        {
            var mu = DenseMatrix.Create(2, 2, 0); //mu*mu = 1
            mu[0, 0] = 2;
            mu[0, 1] = -1;
            mu[1, 0] = 3;
            mu[1, 1] = -2;
            return mu;
        }
    }
}