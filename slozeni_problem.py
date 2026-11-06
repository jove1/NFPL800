#!/usr/bin/env python3
"""
Converting between atomic and weight fractions.
"""

import numpy as np
# Atomic masses (relative atomic masses; numerically equivalent to molar masses
# in g/mol for the calculations in this exercise).
# Source: NIST atomic weights and isotopic compositions.
# https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
atomic_mass_table = {
    "H": 1.008,
    "O": 16.000,
    "Cu": 63.546,
    "Al": 26.982,
}


def at2wt(elements, x):
    """Convert atomic fractions into normalized weight fractions.

    Parameters
    ----------
    elements : sequence of str
        Chemical symbols identifying the components in the composition.
    x : array-like
        Atomic proportions or fractions. The first axis corresponds to
        ``elements``; any additional axes represent separate samples.

    Returns
    -------
    numpy.ndarray
        Weight fractions with the same shape as the input composition.
    """
    # Look up each element's atomic mass. Initially, w has shape
    # (number_of_elements,).
    ma = np.array([atomic_mass_table[e] for e in elements])

    # STUDENT TASK: Convert atomic proportions into weight fractions. Each
    # component's contribution to mass depends on both its atomic amount and
    # its atomic mass; normalize the contributions so they sum to one along the
    # element axis.
    #
    # IMPORTANT: x may contain NumPy arrays, not just single numbers.  For
    # example, [x, 1-x] in the plotting section becomes an array of shape
    # (number_of_elements, number_of_samples). A 1-D mass array w will not
    # automatically align with the first axis of that 2-D array under NumPy's
    # default broadcasting rules. Arrange dimensions as needed so the mass for
    # each element is applied to every sample.
    #
    # Broadcasting tutorial:
    # https://numpy.org/doc/stable/user/basics.broadcasting.html
#PROBLEM
    raise NotImplementedError
#END


def wt2at(elements, x):
    """Convert weight fractions into normalized atomic fractions.

    Parameters
    ----------
    elements : sequence of str
        Chemical symbols identifying the components in the composition.
    x : array-like
        Weight fractions. The first axis corresponds to ``elements``; any
        additional axes represent separate samples.

    Returns
    -------
    numpy.ndarray
        Atomic fractions with the same shape as the input composition.
    """
    # STUDENT TASK:
    # Convert the supplied weight fractions into atomic proportions using the
    # atomic masses in atomic_mass_table, then normalize along the element
    # axis so the atomic fractions for each sample sum to one.
    #
    # As in at2wt, x can be multidimensional. In the plot below, [x, 1-x]
    # represents two components evaluated at many sample points. Ensure that
    # the atomic-mass array broadcasts across the sample dimensions correctly;
    # remember that NumPy aligns dimensions from the right.
    #
    # Broadcasting tutorial:
    # https://numpy.org/doc/stable/user/basics.broadcasting.html
#PROBLEM
    raise NotImplementedError
#END


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    # Show a grid on every plot to make the curves easier to read.
    plt.rc("axes", grid=True)

    # Convert a simple hydrogen/oxygen atomic ratio to weight fractions.
    print(at2wt(["H", "O"], [2, 1]))
    print("should be:", [0.11190053, 0.88809947])

    # Convert the resulting weight fractions back to atomic fractions.
    print(wt2at(["H", "O"], [0.11190053, 0.88809947]))
    print("should be:", [0.66666666, 0.33333334])

    # Sweep the fraction of copper from 0 to 1; the aluminum fraction is
    # 1-x. Each call therefore evaluates many compositions at once.
    x = np.linspace(0, 1)
    plt.figure()
    plt.plot(x, at2wt(["Cu", "Al"], [x, 1 - x]))
    plt.xlabel("Input atomic fraction of Cu")
    plt.ylabel("Weight fraction of Cu and Al")

    plt.figure()
    plt.plot(x, wt2at(["Cu", "Al"], [x, 1 - x]))
    plt.xlabel("Input weight fraction of Cu")
    plt.ylabel("Atomic fraction of Cu and Al")
 
    plt.show()
