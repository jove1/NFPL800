#!/usr/bin/env python3
"""
Piecewise functions and their derivatives.

The script uses NumPy to evaluate functions over arrays of values and
Matplotlib to plot the original function and compare its derivatives with
finite-difference approximations.
"""

import numpy as np

class FT:
    """Represent a piecewise-defined function.

    Arguments are supplied as alternating boundary/value pairs:
        FT(boundary_0, terms_for_interval_0, boundary_1, ...)

    Each ``terms`` argument is a dictionary mapping a basis function to its
    coefficient. For example, ``{T: 1, T2: 1e-3}`` represents
    ``T + 1e-3*T**2``.

    The intervals are selected using the boundaries in ``self.limits``.
    Values that do not match any interval are set to NaN.
    """

    def __init__(self, *args):
        # Every other argument starting at index 0 is a boundary.
        self.limits = args[0::2]
        # The arguments between boundaries describe the function on each
        # corresponding interval.
        self.parts = args[1::2]

    def __call__(self, T, d=0):
        """Evaluate the piecewise function at T.

        Parameters
        ----------
        T : scalar or array-like
            Input value(s).
        d : int, default=0
            Derivative order: 0 means the function itself

        Returns
        -------
        numpy.ndarray or scalar-like
            Values selected according to the piecewise intervals.
        """
        return np.select(
            # For each interval boundary (except the first), check whether
            # each input value lies below that boundary.
            [T < lim for lim in self.limits[1:]],
            # Build the expression for each piece as a weighted sum of its
            # basis functions. Pass d through so each basis function can
            # return the requested derivative.
            [
                sum(f(T, d=d) * coefficient
                    for f, coefficient in part.items())
                for part in self.parts
            ],
            # np.select uses this when none of the conditions is true.
            np.nan,
        )


def C(T, d=0):
    """Constant basis function, equal to 1. The argument T is included so all
    basis functions share the same calling convention; the constant's value
    does not depend on T.
    """
    if d == 0:
        return 1
#PROBLEM
#    # TODO: implement the required derivative behavior for d > 0.
#    raise NotImplementedError
#SOLUTION
    else:
        return 0
#END


def T(T, d=0):
    """Identity basis function f(T) = T."""
    if d == 0:
        return T
#PROBLEM
#    # TODO: implement derivative behavior for d > 0.
#    raise NotImplementedError
#SOLUTION
    elif d == 1:
        return 1
    else:
        return 0
#END


def Tk(k):
    """Create a power basis function f(T) = T**k.

    This is a function factory: calling Tk(k) returns a function ``f`` that
    evaluates the selected power.

    Parameters
    ----------
    k : int or numeric
        Exponent for the power function.
    """
    def f(T, d=0):
        if d == 0:
            return T**k
#PROBLEM
#        # TODO: implement derivative behavior for d > 0.
#        raise NotImplementedError
#SOLUTION
        elif d == 1:
            return k*T**(k-1)
        elif d == 2:
            return k*(k-1)*T**(k-2)
        raise NotImplementedError
#END
    return f


def TLNT(T, d=0):
    """Logarithmic basis function f(T) = T * log(T)."""
    if d == 0:
        return T * np.log(T)
#PROBLEM
#    # TODO: implement derivative behavior for d > 0.
#    raise NotImplementedError
#SOLUTION
    elif d == 1:
        return 1 + np.log(T)
    elif d == 2:
        return 1/T
    raise NotImplementedError
#END


# Named power functions used as basis functions in piecewise expressions.
T2 = Tk(2)
T3 = Tk(3)
T4 = Tk(4)
T7 = Tk(7)
Tm1 = Tk(-1)
Tm9 = Tk(-9)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # Show a grid on every plot to make curves easier to read.
    plt.rc("axes", grid=True)

    # Define one piece of a piecewise function on the interval beginning at
    # 0 and ending at 1000. The expression is:
    #   5000*C(T) + 1*T + 1e-3*T**2 + 2*T*log(T)
    # The coefficients are stored in a dictionary as basis-function keys.
    func = FT(0, {C: 5000, T: 1, T2: 1e-3, TLNT: 2}, 1000)

    # Generate 50 evenly spaced sample points from 0 to 1000.
    t = np.linspace(0, 1000)

    # Plot the function itself (d=0 by default).
    plt.figure()
    plt.plot(t, func(t))
    plt.title("Function")

    # A small step size for finite-difference approximations.
    eps = 1e-4

    # Compare a central finite-difference approximation of the first
    # derivative with the derivative computed by func(..., d=1).
    # Formula: [f(t + eps) - f(t - eps)] / (2*eps).
    plt.figure()
    plt.plot(t, (func(t + eps) - func(t - eps)) / (2 * eps), ".",
             label="Finite-difference approximation")
    plt.plot(t, func(t, d=1), label="Implemented first derivative")
    plt.title("First derivative")
    plt.legend()

    plt.figure()
#PROBLEM
#    # TODO: replace NotImplementedError below with a numerical
#    # central-difference approximation for the second derivative. Then compare
#    # that approximation with func(t, d=2).
#    #
#    # A common central-difference formula for the second derivative uses
#    # f(t + eps), f(t), and f(t - eps). Derive and implement the expression.
#    plt.plot(t, NotImplementedError, ".")
#SOLUTION
    plt.plot(t, ( func(t+eps) - 2*func(t) + func(t-eps) )/eps/eps, ".")
#END
    plt.plot(t, func(t, d=2))
    plt.title("Second derivative")

    plt.show()
