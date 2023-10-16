
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

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
    public static implicit operator Matrix<double>(SuperComplex val)
    {
        return val.ToMatrix();
    }
    public static implicit operator SuperComplex(double[] coef)
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
        var result = s1.ToMatrix() * s2.ToMatrix().Inverse();
        return new(result);
    }
    /// <summary>
    /// Applies function to super complex number, yeah, complex valued result is sufficient.<br/>
    /// For example: to compute sin of super complex number just put Complex.Sin in this method, and
    /// it will apply sin function to whole super complex number
    /// </summary>
    public SuperComplex Apply(Func<Complex,Complex> func){
        var mat = ToMatrix();
        var a = mat[0,0];
        var b = mat[0,1];
        var c = mat[1,0];
        var d = mat[1,1];

        var D = Complex.Sqrt((a+d)*(a+d)-4*(a*d-c*b));
        
        //eigenvalues 
        var k1 = (a+d+D)/2;
        var k2 = (a+d-D)/2;

        //eigenvectors first dimension. Second is 1
        var v1 = b/(k1-a);
        var v2 = b/(k2-a);

        //eigenvectors matrix
        var V = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix.Create(2,2,1);
        V[0,0]=v1;
        V[0,1]=v2;
        var VInverse = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix.Create(2,2,0);
        VInverse[0,0]=1;
        VInverse[0,1]=-v2;
        VInverse[1,0]=-1;
        VInverse[1,1]=v1;

        VInverse*=1.0/(v1-v2);

        var identity = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix.Create(2,2,0);
        identity[0,0]=func(k1);
        identity[1,1]=func(k2);
        var result = V*identity*VInverse;
        return new(result);
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
    double[] castToCoefficients(Matrix<Complex> val)
    {
        var a = val[0, 0];
        var b = val[0, 1];
        var c = val[1, 0];
        var d = val[1, 1];

        var k1 = (a + d) / 2;
        var k2 = (a - d - 5 * b - 3 * c) / 2;
        var k3 = (a - d) / 2 - b - c;
        var k4 = (b + c) / 2;
        return castToCoefficients(new[] { k1, k2, k3, k4 });
    }
    double[] castToCoefficients(Matrix<double> val)
    {
        var a = val[0, 0];
        var b = val[0, 1];
        var c = val[1, 0];
        var d = val[1, 1];

        var k1 = (a + d) / 2;
        var k2 = (a - d - 5 * b - 3 * c) / 2;
        var k3 = (a - d) / 2 - b - c;
        var k4 = (b + c) / 2;
        return new[] { k1, k2, k3, k4 };
    }
    //these values should not be changed, so they are private
    static Matrix<double> imaginary = Imaginary;
    static Matrix<double> real = Real;
    static Matrix<double> theta = Theta;
    static Matrix<double> mu = Mu;
    public override string ToString()
    {
        char sign(double val)=>val>=0 ? '+' : '-';
        return $"({Coefficients[0]} {sign(Coefficients[1])} {Math.Abs(Coefficients[1])}i {sign(Coefficients[2])} {Math.Abs(Coefficients[2])}θ {sign(Coefficients[3])} {Math.Abs(Coefficients[3])}μ)";
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