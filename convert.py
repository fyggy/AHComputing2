import functools as ft
import mpmath as mp
import numpy as np
import scipy.special as sp

from mpmath import fp
from operator import lt, gt, le, ge
from sympy import Symbol, Function, pi, E, EulerGamma, I, gamma, polygamma, digamma, \
    zeta, lerchphi, airyai, airybi, Ci, Ei, Si, li, polylog, csch, sech, coth, asinh, acosh, atanh, acsch, asech, \
    acoth, nsimplify

from helpers import *


# replaces functions which latex2sympy cannot parse with their real sympy equivalents
def replace(expr):

    # replace all decimal numbers with rational numbers
    # if this is not done then decimal numbers can sometimes cause recursion errors
    expr = nsimplify(expr, rational=True)

    # polygamma and digamma is a special case, as they share a symbol
    # attempt to replace \psi with polygamma, otherwise replace it with digamma
    try:
        expr = expr.replace(func_polygamma, polygamma)
    except TypeError:
        expr = expr.replace(func_digamma, digamma)

    # replace each entry in the replace dict (at the bottom of the file) with its respective sympy function
    for query in replaces:
        expr = expr.replace(query, replaces[query])

    return expr


# dataclass to contain code snippets
# essentially a glorified array
class Snippet:
    def __init__(self, name, code):
        self.name = name
        self.code = code


# !! BINARY SEARCH CONTAINED HERE !! #
# class to contain all code snippets
class Code:
    def __init__(self, snippets):
        self.snippets = snippets

    # uses binary search to find the target snippet by name
    # error controls whether it uses the error function, or raises a normal error
    # this is important, as sometimes the error is escaped, but it would still be printed
    def get_by_name(self, target, error):
        bot = 0
        top = len(self.snippets) - 1
        while top >= bot:
            mid = (top + bot) // 2
            current = self.snippets[mid]
            if current.name == target:
                return current
            elif current.name > target:
                top = mid - 1
            else:
                bot = mid + 1

        if error:
            send_error(f"Function {target} is not recognised")
        raise ValueError


# wraps the normal zeta function removing keyword arguments
def _zeta(x, a):
    return mp.zeta(x, a=a)


# various functions from mpmath which must be wrapped
_li = broadcast(fp.li)
_lerchphi = broadcast(error_wrapper(mp.lerchphi))
_zeta = broadcast(error_wrapper(_zeta))
_psi = error_wrapper(fp.psi)


# polygamma is a special case, as m is a constant, and there is a more efficient digamma function that exists within
# scipy, but which doesnt exist for general polygamma
# it also cannot use broadcast, as m must be an integer, whereas broadcast converts all integers into complex numbers
def _polygamma(m, z):

    # convert m from a complex number to an integer, if possible
    m = better_int(m)

    # use more efficient digamma function if possible
    if m == 0:
        return sp.digamma(z)

    # otherwise use mpmath polygamma
    else:

        # if z is an array, then create an output array, otherwise, call function on z normally
        try:
            out = np.zeros(z.shape, dtype=np.complex128)
        except AttributeError:
            return _psi(m, complex(z))
        else:

            # loop through z and evaluate at each input
            for i, z0 in enumerate(z):
                out[i] = _psi(m, complex(z0))

            return out


# polylogarithm is also a special case, as it only takes real numbers as inputs for the base
# this means, like polygamma, broadcast cannot be used, as it converts all floats to complex numbers
def _polylog(s, z):

    # convert s from complex number to a real number, if possible
    s = better_float(s)

    # if z is an array, then create an output array, otherwise, call function on z normally
    try:
        out = np.zeros(z.shape, dtype=np.complex128)
    except AttributeError:
        return fp.polylog(s, complex(z))
    else:

        # loop through z and evaluate at each input
        for i, z0 in enumerate(z):
            out[i] = fp.polylog(s, complex(z0))

        return out


# global variable used in conv
last_name = 0


# function to convert an expression to python code, which can be executed and then evaluated at arbitrary points
# this function has many special cases, due to the complexity of this action
# however, at its core, it is simple
# it simply finds the code for the head, and then substitutes in the code for the arguments
# the arguments themselves have heads, and respective code for their heads, and so on and so forth
def conv(expr, error):

    # set last_name as a global variable which can be used in this function
    global last_name

    # assign variables for use within conv
    # head is the sympy function on the outside of expr
    # args is the arguments for this function, which are likely themselves expressions
    # this means that an expression has a tree structure, and this function navigates this tree structure
    # to read more go here https://docs.sympy.org/latest/tutorial/manipulation.html
    head = expr.func.__name__
    args = expr.args

    # get the code for the head of the function
    head_code = code_snippets.get_by_name(head, error).code

    # if args == (), then expr is a leaf on the tree
    # this constitutes a special case, and the following allows all possible leafs to be evaluated correctly
    if args == ():
        if head in ["Float", "Rational", "Integer", "Symbol", "Dummy"]:
            return head_code.format(str(expr))
        else:
            return head_code

    # add and mul dont list their arguments as the addition of 2 terms, but a addition of n terms
    # for example,
    #       "x + 3 + 4"
    # would be stored as
    #       Add(Symbol("x"), Integer(3), Integer(4))
    # instead of
    #       Add(symbol("x"), Add(Integer(3), Integer(4)))
    # doing this allows Add and Mul to be treated as binary operators
    if head in ["Add", "Mul"]:
        args = expr.as_two_terms()

    # differentiation is difficult
    # since not all expressions can be differentiated analytically, then numerical differentiation is required
    # this case creates 2 new functions which greatly simplify and speed up the code
    # although this doesnt seem more simple at all, this is a significant improvement over my last attempt, especially
    # in readability
    if head == "Derivative":

        # set the names of the 2 new functions
        name1 = "inner_" + str(last_name)
        name2 = "outer_" + str(last_name)
        last_name += 1

        # since a single Derivative may in fact be multiple differentials with respect to multiple variables, then
        # fix these arrays of variables being differentiated with respect to, and how many times
        dummies = [i[0] for i in args[1:]]
        orders =  [i[1] for i in args[1:]]

        # fix array of the dummy variables but in string form
        string_dummies = [str(i) for i in dummies]

        # fix array of the symbols that are *not* being differentiated
        frozen_vars = list(args[0].free_symbols.difference(set(dummies)))

        # convert that array to a string
        frozen_vars_string = [str(i) for i in frozen_vars]

        # get array of the total variables as strings. however, make sure the dummies come first
        total_vars = string_dummies + frozen_vars_string

        # use the create_func function to create a function, with name name1, out of the expression being differentiated
        create_func(args[0], name1, variables=total_vars)

        # memorise the function, as it is being evaluated at many points which are exactly the same
        # this speeds up the function slightly
        # assignment to globals like this, while slightly hacky, is completely safe
        globals()[name1] = mp.memoize(complex_wrapper(globals()[name1]))

        # create another function, this time with name name2
        # this function contains the actual differentiation function, from mpmath, which take a function as an input
        # it uses the wrapper partial from functools to freeze the variables of the function name1, defined above
        # this means that it acts as if the frozen variables are being entered, but not differentiated with respect to
        # the value h is the step used for numerical differentiation. after much testing, i have determined that this
        # value is ideal
        create_func(f"fp.diff(ft.partial({name1}, {', '.join([f'{i}={i}' for i in frozen_vars_string])}), " +
                    f"({', '.join(string_dummies)}, ), n={tuple(orders)}, relative=True, h=0.001)",
                    name2, variables=total_vars, string=True)

        # since mpmath diff does not take arrays as input then it must be broadcasted
        globals()[name2] = broadcast(globals()[name2])

        # return a string that allows the function we just defined to be evaluated
        return f"{name2}({', '.join(total_vars)})"

    # summation is also difficult
    # this case also creates 2 new functions to speed up and simplify code
    # uses a similar principle to differentiation
    if head == "Sum":

        # initialise names for the functions to be defined
        name1 = "inner_" + str(last_name)
        name2 = "outer_" + str(last_name)
        last_name += 1

        # assign arrays of the dummy variables, lower bounds, and upper bounds
        dummies = [i[0] for i in args[1:]]
        lowers  = [i[1] for i in args[1:]]
        uppers  = [i[2] for i in args[1:]]

        # fix array of the dummy variables but in string form
        string_dummies = [str(i) for i in dummies]

        # fix array of the symbols that are *not* being summed
        frozen_vars = list(args[0].free_symbols.difference(set(dummies)))

        # convert that array to a string
        frozen_vars_string = [str(i) for i in frozen_vars]

        # get array of the total variables as strings. however, make sure the dummies come first
        total_vars = string_dummies + frozen_vars_string

        # create a string containing the bounds which is in the right format to be inputted into the function
        bounds = ""
        for lower, upper in zip(lowers, uppers):
            bounds += f"[{better_int(lower)}, {better_int(upper)}]" + ", "
        bounds = bounds[:-2]

        # define a function that is the general term, with name name1, allowing for evaluation in the sum
        create_func(args[0], name1, variables=total_vars)

        # memorise the function, as it will be evaluated at the same point many times
        globals()[name1] = mp.memoize(complex_wrapper(globals()[name1]))

        # create a function using the mpmath sum function which evaluates the sum
        # like before, the partial wrapper from functools is used to freeze the variables not being summed over
        # this means that it acts as if the frozen variables are being entered, but not summed over
        # the tol and maxterms allow fine tuning of the precision without making it take an extreme amount of time to
        # evaluate. i found these values provide a happy medium
        create_func(f"fp.nsum(ft.partial({name1}, {', '.join([f'{i}={i}' for i in frozen_vars_string])}), {bounds}, tol=0.001, maxterms=100)",
                    name2, variables=frozen_vars_string, string=True)

        # like before, broadcast the new function so it can be evaluated at lists of points
        globals()[name2] = broadcast(error_wrapper(globals()[name2]))

        # return a string that allows the function we just defined to be evaluated
        return f"{name2}({', '.join(frozen_vars_string)})"

    # this case is a little strange
    # it occurs in the analytic differentiation of certain tricky functions, such as the Lerch Phi
    # sympy uses it such that it can differentiate a multivariable function with respect to just one of its variables
    # for our purposes, all we need to do is create the function inside which is to be substituted into, and then
    # make a substitution of variables
    # to read more go here https://docs.sympy.org/latest/modules/core.html#subs
    if head == "Subs":

        # set name for the function to be defined
        name = "inner_" + str(last_name)
        last_name += 1

        # hold the tuple to be substituted from, and the tuple to be substituted to
        first = args[1]
        last  = args[2]

        # contains the to's and from's in convenient dictionary form
        replace_dict = {i: j for i, j in zip(first, last)}

        # hold the variables before substitution
        free_vars = args[0].free_symbols

        # hold the variables after substitution
        changed_vars = [replace_dict[j] if j in replace_dict else j for j in free_vars]

        # create a function for the variables to be substituted into
        create_func(args[0], name, variables=[str(i) for i in free_vars])

        # return function with the variables substituted
        return f"{name}({', '.join([str(i) for i in changed_vars])})"

    # convert all the arguments
    # this is what allows the tree structure of expressions to be navigated
    arg_codes = [conv(i, error) for i in args]

    if head == "Piecewise":
        name = "inner_" + str(last_name)
        last_name += 1
        variables = [str(i) for i in expr.free_symbols]
        create_func(head_code.format(*arg_codes), name, variables=variables, string=True, error=error)
        globals()[name] = broadcast(globals()[name])
        return f"{name}({', '.join(variables)})"

    # zeta has a special case because sometimes it has 1 argument and sometimes it has 2
    # having 1 argument implicitly means the second argument is 1
    if head == "zeta" and len(arg_codes) == 1:
        arg_codes.append("1 + 0j")


    # insert the code from the arguments into the code from the head
    return head_code.format(*arg_codes)


# the wrapping code for the function, which allows it to be evaluated
wrapper = """
def {0}({1}):
    return {2}
    """


# take an expression or code and turn it into a function
# creates a function with a name of "name", using variables of "variables"
# if "string" is true, then expr is interpreted as a string of code, otherwise it is interpreted as a sympy expression
# if "error" is true, then it prints errors like normal, otherwise it doesnt print anything
def create_func(expr, name, variables=["x"], string=False, error=True):
    if not string:
        code = conv(expr, error)
    else:
        code = expr

    # check if function is a constant, and, if so, turn it into a numpy array
    if not string and expr.free_symbols == set():
        code = f"np.full(x.shape, ({code}), dtype=np.complex128)"

    # insert parameters into the wrapper
    code = wrapper.format(name, ", ".join(variables), code)

    # print(code)

    # attempt to execute function into the global variables
    try:
        exec(code, globals())

    # otherwise throw error
    except Exception as e:
        if error:
            send_error(f"Fatal Error: {e}")
        raise e


# digamma and polygamma represent special cases
func_digamma = Function("psi")
func_polygamma = Function("psi")

# dict for replacements
replaces = {
    Symbol("pi"): pi,
    Symbol("e"): E,
    Symbol("gamma"): EulerGamma,
    Symbol("i"): I,
    Function("Gamma"): gamma,
    Function("zeta"): zeta,
    Function("Phi"): lerchphi,
    Function("Ai"): airyai,
    Function("Bi"): airybi,
    Function("Ci"): Ci,
    Function("Ei"): Ei,
    Function("Si"): Si,
    Function("li"): li,
    Function("Li"): polylog,
    Function("csch"): csch,
    Function("sech"): sech,
    Function("coth"): coth,
    Function("arcsinh"): asinh,
    Function("arccosh"): acosh,
    Function("arctanh"): atanh,
    Function("arccsch"): acsch,
    Function("arcsech"): asech,
    Function("arccoth"): acoth,
}

# snippets of code
# comments on specific lines below
code_snippets = Code([
    Snippet("Abs", "np.absolute({0})"),
    Snippet("Add", "np.add(({0}), ({1}))"),
    Snippet("BooleanFalse", "False"),
    Snippet("BooleanTrue", "True"),
    Snippet("Catalan", "fp.catalan + 0j"),
    Snippet("Ci", "sp.sici({0})[1]"),
    Snippet("ComplexInfinity", "np.inf + 0j"),
    Snippet("Derivative", ""),
    Snippet("Dummy", "{0}"),
    Snippet("Ei", "sp.expi({0})"),
    Snippet("Equality", "({0}) == ({1})"),
    Snippet("EulerGamma", "fp.euler + 0j"),
    Snippet("Exp1", "fp.e + 0j"),
    Snippet("ExprCondPair", "({0}) if ({1}) else"),
    Snippet("Float", "{0} + 0j"),
    Snippet("GoldenRatio", "fp.phi + 0j"),
    Snippet("GreaterThan", "better_ineq({0}, {1}, gt)"),
    Snippet("Half", "0.5 + 0j"),
    Snippet("ImaginaryUnit", "0 + 1j"),
    Snippet("Infinity", "np.inf"),
    Snippet("Integer", "{0} + 0j"),
    Snippet("LessThan", "better_ineq({0}, {1}, lt)"),
    Snippet("Mul", "np.multiply(({0}), ({1}))"),
    Snippet("NegativeInfinity", "-np.inf"),
    Snippet("NegativeOne", "-1 + 0j"),
    Snippet("One", "1 + 0j"),
    Snippet("Pi", "fp.pi + 0j"),
    Snippet("Piecewise", "{0} {1} None"),
    Snippet("Pow", "np.power(({0}), ({1}))"),
    Snippet("Rational", "{0} + 0j"),
    Snippet("Si", "sp.sici({0})[0]"),
    Snippet("StrictGreaterThan", "better_ineq({0}, {1}, ge)"),
    Snippet("StrictLessThan", "better_ineq({0}, {1}, le)"),
    Snippet("Subs", ""),
    Snippet("Sum", ""),
    Snippet("Symbol", "{0}"),
    Snippet("TribonacciConstant", "1.839286755214161 + 0j"),
    Snippet("Tuple", ""),
    Snippet("Zero", "0 + 0j"),
    Snippet("acos", "np.arccos({0})"),
    Snippet("acosh", "np.arccosh({0})"),
    Snippet("acot", "np.arctan(np.reciprocal({0}))"),
    Snippet("acoth", "np.arctanh(np.reciprocal({0}))"),
    Snippet("acsc", "np.arcsin(np.reciprocal({0}))"),
    Snippet("acsch", "np.arcsinh(np.reciprocal({0}))"),
    Snippet("airyai", "sp.airy({0})[0]"),
    Snippet("airyaiprime", "sp.airy({0})[1]"),
    Snippet("airybi", "sp.airy({0})[2]"),
    Snippet("airybiprime", "sp.airy({0})[3]"),
    Snippet("asec", "np.arccos(np.reciprocal({0}))"),
    Snippet("asech", "np.arccosh(np.reciprocal({0}))"),
    Snippet("asin", "np.arcsin({0})"),
    Snippet("asinh", "np.arcsinh({0})"),
    Snippet("atan", "np.arctan({0})"),
    Snippet("atanh", "np.arctanh({0})"),
    Snippet("beta", "sp.beta(({0}), ({1})"),
    Snippet("cos", "np.cos({0})"),
    Snippet("cosh", "np.cosh({0})"),
    Snippet("cot", "np.reciprocal(np.tan({0}))"),
    Snippet("coth", "np.reciprocal(np.tanh({0}))"),
    Snippet("csc", "np.reciprocal(np.sin({0}))"),
    Snippet("csch", "np.reciprocal(np.sinh({0}))"),
    Snippet("exp", "np.exp({0})"),
    Snippet("factorial", "sp.gamma(({0}) + 1+0j)"),
    Snippet("gamma", "sp.gamma({0})"),
    Snippet("im", "({0}).imag"),
    Snippet("lerchphi", "_lerchphi({0}, {1}, {2})"),
    Snippet("li", "_li({0})"),
    Snippet("log", "np.log({0})"),
    Snippet("loggamma", "sp.loggamma({0})"),
    Snippet("polygamma", "_polygamma(({0}), ({1}))"),
    Snippet("polylog", "_polylog({0}, {1})"),
    Snippet("re", "({0}).real"),
    Snippet("sec", "np.reciprocal(np.cos({0}))"),
    Snippet("sech", "np.reciprocal(np.cosh({0}))"),
    Snippet("sin", "np.sin({0})"),
    Snippet("sinh", "np.sinh({0})"),
    Snippet("tan", "np.tan({0})"),
    Snippet("tanh", "np.tanh({0})"),
    Snippet("zeta", "_zeta(({0}), ({1}))")
])
