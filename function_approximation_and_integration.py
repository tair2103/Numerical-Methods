# Fill in your details
full_name1 = "Noa Cohen"        # Fill full name of 1st student here
student_ID1 = "211561105"            # Fill ID of 1st student here

full_name2 = "Tair Steinberg"        # Fill full name of 2nd student here
student_ID2 = "322736521"               # Fill ID of 2nd student here


import math

# 1. Taylor Series Approximation
def taylor_series_ln_sinx(x0, n, x):
    """
    Compute the Taylor polynomial of ln(1 + sin(x)) of order n around point x0.

    Args:
        x0 (float): Reference point.
        n (int): Order of the Taylor series.
        x (float): Point at which to evaluate the polynomial.

    Returns:
        float: Approximation of ln(1 + sin(x)) using the Taylor series.
    """
    taylor_sum = 0
    for k in range(n + 1):
        if k == 0:
            term = math.log(1 + math.sin(x0))
        elif k == 1:
            term = (math.cos(x0) / (1 + math.sin(x0))) * (x - x0)
        elif k == 2:
            term = (-math.sin(x0) / ((1 + math.sin(x0)) ** 2)) * ((x - x0) ** 2) / 2
        elif k == 3:
            term = ((math.sin(x0) ** 2 - math.cos(x0) ** 2) / ((1 + math.sin(x0)) ** 3)) * ((x - x0) ** 3) / 6
        else:
            term = 0
        taylor_sum += term
    return taylor_sum
    pass


# 2a.Legendre Polynomials
def legendre_polynomial(n, x):
    """
    Compute the Legendre polynomial of degree n at point x.

    Args:
        n (int): Degree of the polynomial.
        x (float): Point at which to evaluate the polynomial.

    Returns:
        float: Value of the Legendre polynomial at x.
    """
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return ((2 * n - 1) * x * legendre_polynomial(n - 1, x) - (n - 1) * legendre_polynomial(n - 2, x))  / n
    pass

# 2.b Legendre Approximation 
def approximate_ln_sinx_with_legendre(n, x, x0):
    """
    Approximate ln(1 + sin(x)) using Legendre polynomials.

    Args:
        n (int): Number of terms in the approximation.
        x (float): Point at which to evaluate the approximation.
        x0 (float): Reference point.

    Returns:
        float: Approximation of ln(1 + sin(x)) using Legendre polynomials.
    """
    a, b, N = -1, 1, 1000
    h = (b - a) / N
    approximation = 0

    for j in range(n + 1):
        integral = 0.5 * (math.log(1 + math.sin(a)) * legendre_polynomial(j, a) +
                          math.log(1 + math.sin(b)) * legendre_polynomial(j, b))

        for i in range(1, N):
            x_i = a + i * h
            integral += math.log(1 + math.sin(x_i)) * legendre_polynomial(j, x_i)

        integral *= h
        lambda_j = (2 * j + 1) / 2 * integral
        approximation += lambda_j * legendre_polynomial(j, x)

    return approximation
    pass

if __name__ == "__main__":
    print (f"This work is the work of:\n{full_name1} &  {full_name2}({student_ID1}, {student_ID2})")

    # Example usage for Taylor series
    x0_taylor = 0
    n_taylor = 2
    x_value = 1
    taylor_result = taylor_series_ln_sinx(x0_taylor, n_taylor, x_value)
    print(f"Taylor series approximation: {taylor_result}")

    # Example usage for Legendre polynomial approximation
    n_legendre = 3
    x0_legendre = 0
    legendre_result = approximate_ln_sinx_with_legendre(n_legendre, x_value, x0_legendre)
    print(f"Legendre polynomial approximation: {legendre_result}")
    


