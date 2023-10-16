# SuperComplexNumber
A general implementation for super complex number, where $i^2 = -1$, $\theta^2=0$, $\mu^2=1$ (non-trivial)

There is some staff in quantum physics that allows you to have such object, that

$$
\theta^2=0
$$

meanwhile $\theta \ne 0$

And also there is a way to define such object that

$$
\mu^2=1
$$

meanwhile $\mu \ne 1$

So here I am presenting the mind blowing implementation for arithemtic for that stuff.

Thanks to some linear algebra machinery, I was able to construct such set, that this super complex number perfectly maps one
two 2 by 2 matrix to exactly one super complex number, so these set of $real, i, \theta, \mu$ is basis for two 2 by 2 matricies,
and I was able to create analytic continuation of any complex-valued function to this super complex numbers.

Look bro

```cs
var real = SuperComplex.Real;
var i = SuperComplex.Imaginary;
var th = SuperComplex.Theta;
var mu = SuperComplex.Mu;

SuperComplex v1 = 3 * real + 4 * i + 2 * th + 4 * mu;
SuperComplex v2 = 6 * real + 1.5 * i - 6 * th - mu;

System.Console.WriteLine(v1);
System.Console.WriteLine(v2);
System.Console.WriteLine(1/v2+2*v2*v2);
System.Console.WriteLine(v2.Apply(x=>1/x+2*x*x));
```

It outputs

```
(3 + 4i + 2θ + 4μ)
(6 + 1.5i - 6θ - 1μ)
(146.29999999999998 + 34.800000000000004i - 139.20000000000002θ - 23.200000000000003μ)
(146.29999999999995 + 34.799999999999955i - 139.20000000000002θ - 23.200000000000003μ)
```

Here in `v2.Apply` I only giving complex valued function, but it gives same output as value computed on super complex plane!

Isn't it insane?

Yeah, so you can take `sin`, `cos`, any function that is defined on complex plane is defined on super complex plane as well!

Code is pretty basic except for `Apply` method, which uses technique from this video to expand matrix-view of
super complex numbers to work on any function.

See https://www.youtube.com/watch?v=-1loSrioE4Y&t=527s
