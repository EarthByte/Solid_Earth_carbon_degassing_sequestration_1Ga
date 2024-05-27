
"""
    Copyright (C) 2017 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


###################################################################################################################################
# Generate a bilinear, biquadratic or bicubic regular grid of z values modelled from a sparse input xyz file using GMT 'trend2d'. #
###################################################################################################################################


from __future__ import print_function
from call_system_command import call_system_command
import math
import sys


def read_input_file(input_filename, log_scale_input_x, log_scale_input_y):
    try:
        with open(input_filename, 'r') as input_file:
            input_file_data = input_file.read()
    except IOError:
        print('Unable to open input data file {0}.'.format(input_filename), file=sys.stderr)
        raise
    
    # Return data unmodified if not using log scaling of xy.
    if not log_scale_input_x and not log_scale_input_y:
        return input_file_data
    
    # ...otherwise convert x,y to log(x),log(y) as requested...
    
    log_scaled_xy_lines = []
    
    for line in input_file_data.splitlines(True):
        if line.strip().startswith('#'):
            continue
        
        line_data = line.split(None, 2)
        num_values = len(line_data)
        
        # If just a line containing white-space then skip to next line.
        if num_values == 0:
            continue
        
        if num_values < 3:
            print('Ignoring line "{0}" - has fewer than 3 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            x = float(line_data[0])
            y = float(line_data[1])
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point x and y values.'.format(line), file=sys.stderr)
            continue
        
        rest_of_line = line_data[2]
        
        transformed_x = math.log10(x) if log_scale_input_x else x
        transformed_y = math.log10(y) if log_scale_input_y else y
        
        log_scaled_xy_lines.append('{0} {1} {2}'.format(transformed_x, transformed_y, rest_of_line))
    
    return ' '.join(log_scaled_xy_lines)


# Shift/scale the data into the range of the original input data (the data that was modelled).
def transform_2d(x, y, offset_x, offset_y, scale_x, scale_y, log_scale_input_x, log_scale_input_y):
    if log_scale_input_x:
        x = math.log10(x)
    if log_scale_input_y:
        y = math.log10(y)
        
    # We need to offset and scale the input data into the range [-1,1].
    # This is what GMT 'trend2d' does to get the data in the range of Chebyshev polynomials.
    x = (x - offset_x) * scale_x
    y = (y - offset_y) * scale_y
    
    return (x, y)


def model_function(x, y, offset_x, offset_y, scale_x, scale_y, log_scale_input_x, log_scale_input_y, m):
    # Shift/scale the data into the range of the input data passed to GMT 'trend2d'.
    x, y = transform_2d(x, y, offset_x, offset_y, scale_x, scale_y, log_scale_input_x, log_scale_input_y)
    
    # The Chebyshev linear polynomial is T1(x) = x
    # The Chebyshev quadratic polynomial is T2(x) = 2*x*x - 1
    # The Chebyshev cubic polynomial is T3(x) = 4*x*x*x - 3*x
    #
    # Note: When the model number is less than 10 the last (10 - model_number) coefficients
    # will be zero and have no effect - we do this to avoid 10 'if' statements.
    model = (m[0] +
        m[1] * x +
        m[2] * y +
        m[3] * x*y +
        m[4] * (2*x*x - 1) +
        m[5] * (2*y*y - 1) +
        m[6] * (4*x*x*x - 3*x) +
        m[7] * (2*x*x - 1) * (y) +
        m[8] * (x) * (2*y*y - 1) +
        m[9] * (4*y*y*y - 3*y))
    
    return model


def gmt_trend2d(input_data, use_weights, model_number):
    
    # The command-line strings to execute GMT 'trend2d'.
    trend2d_command_line = ["gmt", "trend2d", "-Fxym", "-N{0}r".format(model_number), "-V"]
    if use_weights:
        trend2d_command_line.append("-W")
    
    stdout_data, stderr_data = call_system_command(trend2d_command_line, stdin=input_data, return_stdout=True, return_stderr=True)
    
    # On windows, with some versions of GMT, stderr does not get properly redirected for some reason and hence
    # 'stderr_data' is empty (which is where the model coefficients are).
    # In this situation you can run GMT manually on the command (eg, 'gmt trend2d -Fxym -N4r -V') and then
    # copy'n'paste the model coefficients here. It's dodgy, but works as a temporary hack.
    #
    #stderr_data = 'trend2d: Model Coefficients: 1.46527976937      1.09140373664   -0.606654695159 0.4803232879'
    
    # With some versions of GMT on Windows the stderr output is lost when redirecting to a file/pipe for some reason.
    if not stderr_data.strip():
        # See workaround above.
        raise RuntimeError('Stderr output lost when redirected from GMT "trend2d" command - this appears to happen on Windows with some GMT versions.')
    
    #print('Stdout: {0}'.format(stdout_data))
    #print('Stderr: {0}'.format(stderr_data))
    
    stdout_lines = stdout_data.splitlines()
    
    # Find the bilinear model coefficients by searching for the 'Model Coefficients:' string in the GMT verbose output.
    model_coefficients_string = 'Model Coefficients:'
    stderr_lines = stderr_data.splitlines()
    for line in stderr_lines:
        model_coefficients_string_index = line.find(model_coefficients_string)
        if model_coefficients_string_index >= 0:
            model_coefficients_data_index = model_coefficients_string_index + len(model_coefficients_string)
            model_coefficients_data_string = line[model_coefficients_data_index:]
            model_coefficients_data = model_coefficients_data_string.split()
            
            if len(model_coefficients_data) != model_number:
                raise RuntimeError('Expected {0} model coefficients from GMT "trend2d" command.'.format(model_number))
            
            try:
                # Convert strings to numbers.
                model_coefficients = [float(coeff_str) for coeff_str in model_coefficients_data]
            except ValueError:
                raise RuntimeError('Could not extract model coefficients from GMT "trend2d" command.')
            
            # Extend the list with zeroes if it's not the maximum length of 10 (supported by trend2d).
            # This is because we always use 10 coefficients in our model function.
            if model_number < 10:
                model_coefficients.extend([0.0] * (10 - model_number))
            
            #print('Model coefficients: {0}'.format(model_coefficients))
            
            return (stdout_lines, model_coefficients)
    
    raise RuntimeError('Could not find model coefficients from GMT "trend2d" command.')


def calc_input_data_min_max(input_data):
    have_input_data = False
    
    for line in input_data.splitlines():
        
        if line.strip().startswith('#'):
            continue
        
        line_data = line.split()
        num_values = len(line_data)
        
        # If just a line containing white-space then skip to next line.
        if num_values == 0:
            continue
        
        if num_values < 2:
            print('Ignoring line "{0}" - has fewer than 2 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            x = float(line_data[0])
            y = float(line_data[1])
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point x and y values.'.format(line), file=sys.stderr)
            continue
        
        if not have_input_data:
            # First line of values.
            x_min = x
            x_max = x
            y_min = y
            y_max = y
            have_input_data = True
        else:
            if x < x_min:
                x_min = x
            if x > x_max:
                x_max = x
            if y < y_min:
                y_min = y
            if y > y_max:
                y_max = y
    
    if not have_input_data:
        raise RuntimeError('Input file does not contain any x and y values.')
    
    if (x_max - x_min) == 0 or (y_max - y_min) == 0:
        raise RuntimeError('Input file contains x or y values whose min and max are equal.')
    
    return x_min, y_min, x_max, y_max


def model_grid_trend2d(
        input_filename,
        log_scale_input_x,
        log_scale_input_y,
        use_weights,
        model_number):
    
    input_data = read_input_file(input_filename, log_scale_input_x, log_scale_input_y)
    
    min_x, min_y, max_x, max_y = calc_input_data_min_max(input_data)
    
    offset_x = 0.5 * (min_x + max_x)
    offset_y = 0.5 * (min_y + max_y)
    scale_x  = 2.0 / (max_x - min_x)
    scale_y  = 2.0 / (max_y - min_y)
    
    _, model_coefficients = gmt_trend2d(input_data, use_weights, model_number)
    
    return model_coefficients, offset_x, offset_y, scale_x, scale_y, min_x, min_y, max_x, max_y


def generate_regular_grid(
        output_min_x, output_min_y,
        output_max_x, output_max_y,
        output_num_x, output_num_y,
        offset_x, offset_y,
        scale_x, scale_y,
        log_scale_input_x, log_scale_input_y,
        model_coefficients):

    #print('Range: {0}->{1}:{2} {3}->{4}:{5}'.format(output_min_x, output_max_x, output_num_x, output_min_y, output_max_y, output_num_y))
    
    step_x = (output_max_x - output_min_x) / (output_num_x - 1)
    step_y = (output_max_y - output_min_y) / (output_num_y - 1)
    
    regular_data = []
    
    for x_index in range(output_num_x):
        x = output_min_x + x_index * step_x
        
        for y_index in range(output_num_y):
            y = output_min_y + y_index * step_y
            
            z = model_function(x, y, offset_x, offset_y, scale_x, scale_y, log_scale_input_x, log_scale_input_y, model_coefficients)
            
            regular_data.append((x, y, z))
    
    return regular_data


def write_output_file(data_filename, data):
    with open(data_filename, 'w') as output_file:
        for x, y, z in data:
            output_file.write('{0} {1} {2}\n'.format(x, y, z))


def grid_trend2d(input_filename, log_scale_input_x, log_scale_input_y, use_weights, model_number, output_filename, xrange, yrange, num_x, num_y):
    
    model_coefficients, offset_x, offset_y, scale_x, scale_y, x_min, y_min, x_max, y_max = model_grid_trend2d(
        input_filename,
        log_scale_input_x,
        log_scale_input_y,
        use_weights,
        model_number)
    
    # Use specified x range (if requested) otherwise use input data range.
    if xrange:
        output_x_min = xrange[0]
        output_x_max = xrange[1]
        # Avoid logarithmic errors.
        if log_scale_input_x and (output_x_min <= 0 or output_x_max <= 0):
            raise RuntimeError('Output x range must be positive and non-zero when using logarithm.')
    elif log_scale_input_x:
        # Undo the logarithm.
        output_x_min = math.pow(10, x_min)
        output_x_max = math.pow(10, x_max)
    else:
        output_x_min = x_min
        output_x_max = x_max
    
    # Use specified y range (if requested) otherwise use input data range.
    if yrange:
        output_y_min = yrange[0]
        output_y_max = yrange[1]
        # Avoid logarithmic errors.
        if log_scale_input_y and (output_y_min <= 0 or output_y_max <= 0):
            raise RuntimeError('Output y range must be positive and non-zero when using logarithm.')
    elif log_scale_input_y:
        # Undo the logarithm.
        output_y_min = math.pow(10, y_min)
        output_y_max = math.pow(10, y_max)
    else:
        output_y_min = y_min
        output_y_max = y_max
    
    output_data = generate_regular_grid(
        output_x_min, output_y_min,
        output_x_max, output_y_max,
        num_x, num_y,
        offset_x, offset_y,
        scale_x, scale_y,
        log_scale_input_x, log_scale_input_y,
        model_coefficients)

    write_output_file(output_filename, output_data)


if __name__ == '__main__':
    
    import argparse
    
    __description__ = \
    """Generate a bilinear, biquadratic or bicubic regular grid of z values modelled from a sparse input xyz file.
    
    The output file is also xyz but the xy values are sampled in a regular grid and the z values
    are obtained from a bilinear/biquadratic/bicubic model of the input data obtained from GMT 'trend2d'.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -lx -ly -m 3 -xc 100 -yc 100 -xr 0 200 -yr 0 100 -- input.xy output.xy
     """

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-lx', '--log_scale_input_x',
            action='store_true',
            help="Model using the logarithm (log10) of the input x data (defaults to modelling the original input x data). "
                "NOTE: The input range must be positive and non-zero if using logarithm.")
    parser.add_argument('-ly', '--log_scale_input_y',
            action='store_true',
            help="Model using the logarithm (log10) of the input y data (defaults to modelling the original input y data). "
                "NOTE: The input range must be positive and non-zero if using logarithm.")
    
    parser.add_argument('-w', '-W', dest='use_weights',
            action='store_true',
            help="Use weights supplied in the 4th column of input file (see GMT trend2d docs on '-W' option).")
    
    parser.add_argument('-m', '--model_number', type=int, choices=[1,2,3,4,5,6,7,8,9,10], default=4,
            help="The model number to fit with (in range [1,10] - see GMT trend2d docs on '-N' option). Default is 4 (bilinear).")
    
    parser.add_argument('-xr', '--xrange', type=float, nargs=2,
            metavar=('xmin', 'xmax'),
            help="The range of x values to grid over (note: this is not affected by logarithmic scaling). "
            "If not specified then the min/max x range of input data is used.")
    parser.add_argument('-yr', '--yrange', type=float, nargs=2,
            metavar=('ymin', 'ymax'),
            help="The range of y values to grid over (note: this is not affected by logarithmic scaling). "
            "If not specified then the min/max y range of input data is used.")
    
    def positive_integer_ge_2(string):
        count = int(string)
        if count < 2:
            msg = '{0} is not an integer >= 2'.format(string)
            raise argparse.ArgumentTypeError(msg)
        return count
    
    parser.add_argument('-xc', '--xcount', type=positive_integer_ge_2, default=20,
            metavar='xcount',
            help="The regular grid count along the x direction (defaults to 20). Must be >= 2.")
    
    parser.add_argument('-yc', '--ycount', type=positive_integer_ge_2, default=20,
            metavar='ycount',
            help="The regular grid count along the y direction (defaults to 20). Must be >= 2.")
    
    parser.add_argument('input_filename', type=str,
            metavar='input_filename',
            help='The input xyz.')
    
    parser.add_argument('output_filename', type=str,
            metavar='output_filename',
            help='The output xyz.')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    grid_trend2d(
            args.input_filename,
            args.log_scale_input_x, args.log_scale_input_y,
            args.use_weights,
            args.model_number,
            args.output_filename,
            args.xrange, args.yrange,
            args.xcount, args.ycount)
        
    sys.exit(0)
