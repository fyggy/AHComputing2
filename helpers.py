from sympy.core.compatibility import as_int
from sympy import oo
from numpy import zeros, complex128, nan, inf, ndarray


# send an error back to php
def send_error(msg):
    print("E" + msg)
    quit()


# generator function to traverse an expr and yield everything of an object
# for example, if head = sy.Integral, then preorder_stop will traverse expr and yield all integrals in the expression
# however, if an integral is nested within another integral, it will not yield that integral
def preorder_stop(expr, head):
    if expr.func == head:
        yield expr
    else:
        for arg in expr.args:
            for i in preorder_stop(arg, head):
                yield i


# create an array of complex numbers spaced equally by lengths of magnitude step
# different to crange as the step can be specified, and the final number will likely not be equally spaced from the
# real end point
def cdist(start, end, step):
    delta = end - start
    difference = (delta / abs(delta)) * step

    num = int(abs(end - start) / step) + 1

    out = zeros(num, dtype=complex128)
    for i in range(num):
        out[i] = start
        start += difference

    return out


# determine if a sympy number is an integer
def isint(n):
    try:
        as_int(n, strict=False)
    except ValueError:
        return False
    else:
        return True


# determine if a sympy number is a real number
def isfloat(n):
    try:
        float(n)
    except TypeError:
        return False
    else:
        return True


# wraps functions, meaning that if they throw an error, then it is caught and a sufficient replacement is outputted
def error_wrapper(func):

    # define wrapper
    def inner(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValueError:
            return nan
        except OverflowError:
            return inf
        except ZeroDivisionError:
            return nan

    return inner


# allows the comparison of complex numbers if their imaginary parts are 0. if not, then it throws an error like usual
def better_ineq(a, b, ineq):
    if type(a) == complex and a.imag == 0:
        a = a.real
    if type(b) == complex and b.imag == 0:
        b = b.real

    # print(ineq(a, b))

    return ineq(a, b)


# allows converting complex numbers to integers if their imaginary parts are 0. also doesnt throw an error on infinity,
# and instead returns it like usual
def better_int(n):
    if type(n) == complex and n.imag == 0:
        n = n.real

    # attempt to interpret as an int
    try:
        return int(n)

    # otherwise, attempt to interpret as infinity
    except (OverflowError, TypeError):
        if n == oo:
            return inf
        elif n == -oo:
            return -inf
        else:
            raise TypeError(f"{n} is not int or inf")


# allows converting complex numbers to floats if their imaginary parts are 0. also doesnt throw an error on infinity,
# and instead returns it like usual
def better_float(n):
    if type(n) == complex and n.imag == 0:
        n = n.real

    # attempt to interpret as a float
    try:
        return float(n)

    # otherwise, attempt to interpret as infinity
    except (OverflowError, TypeError):
        if n == oo:
            return inf
        elif n == -oo:
            return -inf
        else:
            raise TypeError(f"{n} is not float or inf")


# wrapper allow functions that only take single numbers as inputs to accept arrays of numbers as input
# assumes that all input arrays are of the same size, which is a reasonable assumption, since they all come from the
# same initial input arrays
def broadcast(func):

    # define inner function
    def inner(*args):
        casters = []
        constants = []

        # separate constants from lists
        for i in args:
            if isinstance(i, ndarray):
                casters.append(i)
            else:
                constants.append(i)

        # assign output array
        try:
            out = zeros(len(casters[0]), dtype=complex128)

        # this means that there are no casters, so function can be called normally
        except IndexError:
            return func(*args)

        # iterate through casters, splitting them up into singular args that the function can take as input
        for i, arg in enumerate(zip(*casters)):

            # call complex on the inputs, as mpmath often incorrectly handles numpy complex numbers
            out[i] = func(*[complex(j) for j in arg], *[complex(j) for j in constants])
        return out

    return inner


# wraps a function, making sure that all of its outputs are complex numbers
def complex_wrapper(func):

    # define inner function with arbitrary arguments
    def inner(*args, **kwargs):
        return complex(func(*args, **kwargs))

    return inner


# allow rounding of complex numbers. simply rounds the real part and the imaginary part
def better_round(x, deg=15):
    return complex(round(x.real, deg), round(x.imag, deg))


# workaround for another strange problem
# for unrecognised functions, mathquill wraps them in operatorname{<funcname>}
# for example Gamma(x) would be coded as operatorname{\Gamma}\left(x\right)
# this function removes this operatorname and its associated brackets
def remove_bracketed(latex, target):
    out = ""
    i = 0
    after = False

    # loop through each character in latex
    # uses while loop so index can be modified on the fly
    while i < len(latex):
        char = latex[i]

        # if we come across our target, skip over it, and skip over the next curly brace
        if latex[i:].startswith(target):
            i += len(target)
            after = True

        # if this is the next curly brace, skip over it
        elif after and char == "}":
            after = False

        # otherwise, add this character to the output
        else:
            out += char
        i += 1
    return out


# a workaround for a strange problem
# mathquill uses latex, and under latex, powers are supposed to be coded like
# this: x^{2}
# however, mathquill codes single digit powers like this: x^2
# this only happens with single digits, for example,
#        2x
#      x
# will be coded as x^{2x}
# this behaviour cannot be turned off, so instead, we shall use this workaround
def add_brackets(latex):
    num = 0
    out = ""

    # loop through each character in latex
    for i, char in enumerate(latex):

        # add any closing brackets if needed
        if num == 1:
            out += "}"

        # check if brackets are missing
        if char == "^" and not latex[i + 1] == "{":

            # since this only occurs on single digit powers, then 3 is sufficient
            num = 3

        # add current character to the output
        out += char

        # add opening brackets if needed
        if num == 3:
            out += "{"

        # reduce num by one if needed
        if num:
            num -= 1

    # add brackets at the very end of an expression
    if num:
        out += "}"

    return out
