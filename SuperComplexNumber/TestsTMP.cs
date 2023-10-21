using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

//here i temporary store tests
public class TestsTMP
{
    static Matrix<double> real = SuperComplex.Real;
    static Matrix<double> i = SuperComplex.Imaginary;
    static Matrix<double> th = SuperComplex.Theta;
    static Matrix<double> mu = SuperComplex.Mu;
    public static void Square()
    {
        SuperComplex v = 6 * th + 3 * real + 8*i;
        var a = v.Coefficients[0];
        var b = v.Coefficients[1];
        var c = v.Coefficients[2];
        SuperComplex pp = (a * a - b) * real + 2 * a * b * i + 2 * a * c * th + b * c * (i * th + th * i);
        //assert eq
        System.Console.WriteLine(pp);
        System.Console.WriteLine(v * v);
    }
    public static void Sqrt()
    {
        SuperComplex v = 2 * th + 3 * real + i;
        // v/=v.Magnitude;
        var res = v.ApplyComplex(Complex.Sqrt);
        System.Console.WriteLine(res);
        System.Console.WriteLine("Expect\t" + v.ToString("0.0000"));
        System.Console.WriteLine("Actual\t" + (res * res).ToString("0.0000"));

    }
}
