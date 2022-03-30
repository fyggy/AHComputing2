from numpy import linspace, isnan, isinf, complex128
from helpers import better_round, better_float


# class to hold a single point
class Point:
    def __init__(self, inpt, output, derivative):
        self.input = inpt
        self.output = output
        self.derivative = derivative


# class to represent a point which has been deleted
class DeletedPoint(Point):
    input = None
    output = None
    derivative = None


# class to represent part of a line after transformation which has been broken up
class LinePart:
    def __init__(self, points, args):
        self.points = points
        self.args = args

    # calculates the control point for a pair of points
    @staticmethod
    def convert_single(z0, z1, d0, d1, direction):

        # z0 = better_round(z0, deg=8)
        # z1 = better_round(z1, deg=8)

        # extract real and imaginary parts for ease of use
        x0, y0 = z0.real, z0.imag
        x1, y1 = z1.real, z1.imag

        # rotate derivatives appropriately
        d0 *= direction
        d1 *= direction

        # d0 = better_round(d0, deg=8)
        # d1 = better_round(d1, deg=8)

        # if both derivatives are exactly vertical
        if d0.real == 0 and d1.real == 0:
            tmp = (z0 + z1) / 2
            return tmp.real, tmp.imag

        # if first derivative is exactly vertical
        elif d0.real == 0:

            m1 = d1.imag / d1.real
            m1 = better_round(m1, deg=8)

            # better_float ensures that no complex numbers are entered
            return better_float(x0), better_float((m1 * (x0 - x1)) + y1)

        # if second derivative is exactly vertical
        elif d1.real == 0:

            m0 = d0.imag / d0.real
            m0 = better_round(m0, deg=8)

            # better_float ensures that no complex numbers are entered
            return better_float(x1), better_float((m0 * (x1 - x0)) + y0)

        else:

            # convert from derivatives to gradients
            m0 = d0.imag / d0.real
            m1 = d1.imag / d1.real
            m0 = better_round(m0, deg=8)
            m1 = better_round(m1, deg=8)

            # if both derivatives are parallel
            if m0 == m1:
                tmp = (z0 + z1) / 2
                return tmp.real, tmp.imag

            # otherwise calculate normally
            else:
                tmp = (m0 * x0 - m1 * x1 + y1 - y0) / (m0 - m1)

                # better_float ensures that no complex numbers are entered
                return better_float(tmp), better_float(m0 * (tmp - x0) + y0)

    # transform from "normal" coordinate system to the coordinate system on the canvas
    @staticmethod
    def transform(point, centerx, centery, clientw, clienth, radius):
        # split point into its coordinates
        w, x, y, z = tuple(point)

        # scaling factor
        k = clientw / (2 * radius)

        # scale and shift coordinates
        tw = k * (w - centerx) + (clientw / 2)
        tx = k * (x - centery) + (clienth / 2)
        ty = k * (y - centerx) + (clientw / 2)
        tz = k * (z - centery) + (clienth / 2)
        return [tw, tx, ty, tz]

    # convert entire line
    def convert(self, direction):
        # edge case if line has no points
        if len(self.points) == 0:
            return [[]]

        # initialise output array
        # output array is interweaved between normal points and control points, like so
        # [point, control point, point, control point, ... , point]
        output = [[0, 0, 0, 0]] * (2 * len(self.points) - 1)

        # initialise after, for if there is exactly 1 point
        after = self.points[0]
        w1 = after.input
        z1 = after.output

        for i in range(len(self.points) - 1):

            # initialise values into variables
            current = self.points[i]
            after = self.points[i + 1]

            w0 = current.input
            w1 = after.input
            z0 = current.output
            z1 = after.output
            d0 = current.derivative
            d1 = after.derivative

            # take average of input point
            avgw = (w0 + w1) / 2

            # assign point, and control point respectively
            output[2 * i] = LinePart.transform([w0.real, w0.imag, z0.real, z0.imag], *self.args)
            output[2 * i + 1] = LinePart.transform([avgw.real, avgw.imag] +
                                                   list(LinePart.convert_single(z0, z1, d0, d1, direction)), *self.args)

        # assign final point
        output[-1] = LinePart.transform([w1.real, w1.imag, z1.real, z1.imag], *self.args)
        return output


# class which defines an entire line
class Line:
    def __init__(self, start, end, step, function, dfunction, args):
        self.start = start
        self.end = end
        self.step = step
        self.function = function
        self.dfunction = dfunction
        self.args = args

        # initialise input array for this line
        self.input = linspace(start, end, num=step, dtype=complex128)

        # initialise output point array for this line
        self.points = [Point] * len(self.input)

    # calculate each point of the function
    def calculate(self):
        # calculate function and its derivative
        # since both function and dfunction are able to be called on arrays of numbers, then this call is quite easy
        output = self.function(self.input)
        doutput = self.dfunction(self.input)

        # insert into point array, the appropriate points
        for i, (inp, out, dout) in enumerate(zip(self.input, output, doutput)):
            self.points[i] = Point(inp, out, dout)

    # remove values from points that are infinite, or nan
    def trim(self):
        for i, point in enumerate(self.points):
            out = point.output
            dout = point.derivative
            if isinf(out) or isnan(out) or isinf(dout) or isnan(dout):

                # replace an invalid point with a deleted point
                self.points[i] = DeletedPoint
                # print("del")

    # break the line when it reaches a discontinuity
    # this is a staticmethod because it is called recursively
    @staticmethod
    def break_up(points, args):

        # returns an empty array if there are no points to break
        if len(points) == 0:
            return []

        output = []
        for i in range(len(points) - 1):

            # initialise points
            current = points[i]
            after = points[i + 1]

            # if a deleted point is encountered, then break up to the next non-deleted point
            if current == DeletedPoint or after == DeletedPoint:
                start = i

                # loop until the next non deleted point
                while current == DeletedPoint or after == DeletedPoint:
                    i += 1
                    current = points[i]
                    try:
                        after = points[i + 1]
                    except IndexError:
                        break

                # if line starts with a deleted point, then add an empty line part
                if start == 0:
                    output.append(LinePart([], args))

                # otherwise add a line part containing all points before the last deleted point
                else:
                    output.append(LinePart(points[:start+1], args))

                # then break the next points
                output += Line.break_up(points[i+1:], args)
                return output

            # define the small variable s that is used in the check expression
            s = after.input - current.input

            # calculate the check expression
            check = ((abs(after.output - (current.output + (s * current.derivative)))) ** 2) / (abs(current.output) + s)

            # compare the check expression with the average of the derivative
            # obviously this is not the average of the derivative, but i have found these parameters to work very
            # reliably
            if abs(check) >= (abs((current.derivative + after.derivative) + 10) / 20):
                output.append(LinePart(points[:i+1], args))
                output += Line.break_up(points[i+1:], args)
                return output

        # otherwise line is unbroken, therefore, add the entire line as a line part
        output.append(LinePart(points, args))
        return output

    # a non-staticmethod version of break_up, allowing it to be called as a mutator
    def break_fully(self):
        self.lineparts = Line.break_up(self.points, self.args)

    # converts all the contained line parts
    def convert(self, direction):
        output = []
        for i in self.lineparts:
            output.append(i.convert(direction))

        return output
