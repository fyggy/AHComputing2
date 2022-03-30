#!/usr/bin/env python
import sys

from latex2sympy_custom4.process_latex import process_sympy
from sympy import Sum

from convert import *
from lines import *


# allows for the break algorithm, which is recursive, to break an entire line (with max size of 5000) without throwing
# an error
sys.setrecursionlimit(5100)

# get arguments from command line
args = sys.argv

# parse arguments
latex = args[1]
botx, topx = float(args[2]), float(args[3])
boty, topy = float(args[4]), float(args[5])
linestep = float(args[6])
precision = int(args[7])
centerx, centery = float(args[8]), float(args[9])
radius = float(args[10])
clientw, clienth = float(args[11]), float(args[12])

# assign precision settings
if precision == 0:
    step = 50
elif precision == 1:
    step = 100
elif precision == 2:
    step = 500
elif precision == 3:
    step = 1000
else:
    step = 5000

# remove \left and \right from brackets that are included by latex but are not parsed properly for some reason
latex = latex.replace(r"\right", "")
latex = latex.replace(r"\left", "")

# remove spaces from latex which are also not parsed properly for some reason
latex = latex.replace(r"\ ", "")

# remove extraneous operator names
latex = remove_bracketed(latex, "operatorname")

# add brackets around single digit powers powers
latex = add_brackets(latex)

# parse expression in sympy expression
try:
    expr = process_sympy(latex)
except:
    send_error("Invalid Syntax")

# print(expr.doit())

# initialise x as a symbol
x = Symbol("x")

# call replace function from convert.py
# replaces all unrecognised functions that the parser was not able to parse with the real sympy functions
expr = replace(expr)

# check that expression does not have too many or too few (zero) variables
if len(expr.free_symbols) < 1:
    send_error("Too Few Variables! I don't know what to do with this!")

elif len(expr.free_symbols) > 1:
    send_error(f"Too Many Variables: {str(expr.free_symbols)[1:-1]}")

# make sure that an expression's only variable is x
elif expr.free_symbols != {x}:
    send_error(f"Please use the variable x")


# traverse expression and check if variables are valid at all points
def check_bounds(expr, namespace=[x]):

    # traverse all summations
    for i in preorder_stop(expr, Sum):
        lims = i.limits

        # if this is one sum, then limits is a tuple of a variable, lower bound, and upper bound
        # if it it multiple sums, then it is a tuple of tuples
        # this converts single sums to the same format
        if type(lims[0]) == Symbol:
            lims = [lims]

        # copy the namespace so that it can be reverted to on the outside
        tmp_namespace = namespace[:]

        for lim in lims[::-1]:
            dummy = lim[0]
            lower = lim[1]
            upper = lim[2]

            # make sure bounds are evaluated
            lower_done = lower.doit()
            upper_done = upper.doit()

            # check if dummy variable has already been used
            if dummy in tmp_namespace:
                send_error(f"Dummy variable {dummy} has already been used")
                return False

            # make sure that bounds are either integers, or infinity
            elif not (isint(lower_done) or lower_done == oo or lower_done == -oo):
                send_error(f"Lower bound {lower} is not an integer or infinity")
                return False
            elif not (isint(upper_done) or upper_done == oo or upper_done == -oo):
                send_error(f"Upper bound {upper} is not an integer or infinity")
                return False

            # check the bounds of the lower bound, upper bound, and the general term
            else:
                tmp = tmp_namespace[:] + [dummy]
                check_bounds(i.args[0], namespace=tmp)
                check_bounds(lower, namespace=tmp)
                check_bounds(upper, namespace=tmp)

            tmp_namespace.append(dummy)

    # make sure that the degree of all polygammas is a constant positive integer
    for i in preorder_stop(expr, polygamma):
        degree = i.args[0]

        # make sure that degree is a constant, and check if it a positve integer
        if not (degree.free_symbols == set() or isint(degree) or degree >= 0):
            send_error(f"Polygamma degree is not a constant positive integer")
            return False

        # check the bounds of any expressions inside the polygamma
        check_bounds(i.args[1], namespace=namespace)

    # make sure that the base of all polylogarithms is a real number (which cannot be 0)
    # since the input numbers are all complex, then it must be a constant
    for i in preorder_stop(expr, polylog):
        base = i.args[0]

        # check if the base is constant, if the base is a real number, or if the base is 0
        if not(base.free_symbols == set() and isfloat(base) and base != 0):
            send_error("Polylogarithm base must be a real nonzero number")
            return False

        check_bounds(i.args[1], namespace=namespace)

    # if expression passes all checks, return true
    return True


# theoretically, this never executes, but it is here just in case
if not check_bounds(expr):
    send_error(f"Expression {expr} is not valid")

# sometimes, expressions evaluate into functions which cannot be resolved by create_func. in this case, we simply catch
# the error, and convert the unevaluated expressions, which can always be converted
try:
    done_expr = expr.doit()

    # turn off error outputs here, as they are caught and ignored if they occur
    create_func(done_expr, "f", error=False)
    create_func(done_expr.diff(x), "df", error=False)
except:
    create_func(expr, "f")
    create_func(expr.diff(x), "df")
# since create_func is in a different file, then it defines f and df in convert.py's global context, meaning they must
# be imported
from convert import f, df

# after much experimentation, 10 decimal places of precision is necessary
mp.mp.dps = 10

# set the arguments for the transform function in lines.py
args = (centerx, centery, clientw, clienth, radius)

# create horizontal lines
horizontal = []
starts = cdist(complex(botx, boty), complex(botx, topy), linestep)
for start, end in zip(starts, starts - botx + topx):
    horizontal.append(Line(start, end, step, f, df, args))

# create vertical lines
vertical = []
starts = cdist(complex(botx, boty), complex(topx, boty), linestep)
for start, end in zip(starts, starts + 1j * (topy - boty)):
    vertical.append(Line(start, end, step, f, df, args))

linepoints = {"horizontal": [], "vertical": []}

# calculate, trim, break, and convert the horizontal lines
for line in horizontal:
    line.calculate()
    line.trim()
    line.break_fully()

    # these lines are in the horizontal direction, so the derivative does not need to be rotated
    linepoints["horizontal"].append(line.convert(1))

# calculate, trim, break, and convert the vertical lines
for line in vertical:
    line.calculate()
    line.trim()
    line.break_fully()

    # these lines are in the vertical direction, so the derivative must be rotated to be in the vertical direction
    linepoints["vertical"].append(line.convert(1j))

# output converted lines to php
print(linepoints)
