from cmath import exp as ec
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.integrate import quad
from pylab import *
from matplotlib.pyplot import title
n = 50

def lin_zero(x1, n=n):
    x = [x1*s/n for s in range(n)]
    x.append(x1)
    x = [float(s) for s in x]
    return x

def lin_mirror(x1, n=n):
    x = lin_zero(n, x1)
    y = x[::-1]
    y = [-1*s for s in y]
    y = y[:-1]
    y.extend(x)
    y = [float(s) for s in y]
    return y

def lin_space_zero(x1, x2, n=n):
    delta = (x2 - x1)/n
    x = [x1 + s*delta for s in range(n)]
    x.append(x2)
    x = [float(s) for s in x]
    return x

def the_unit_circle(s=n, a=2, b=1, plots=True):
    from sympy import pi, cos, sin, sqrt
    theta = lin_space_zero(0, 2*pi, s-1)

    if plots is True:
        e = [ec(1j*s) for s in theta]
        plot([a*s.imag for s in e], [b*s.real for s in e]) # If you want to rotate it.
        plot([b*s.real for s in e], [a*s.imag for s in e])
        plot([-a*s.imag for s in e], [-b*s.real for s in e])
        plot([-b*s.real for s in e], [-a*s.imag for s in e])
        title('The Unit Circle - Parameterization')
        show()

    x = [b*cos(s) for s in theta]
    y = [a*sin(s) for s in theta]
    if plots is True:
        plot(x, y)
        title('Alternative Unit Circle')
        show()

    x = [b*cos(s) for s in theta]
    y = [a*sin(s) for s in theta]
    r = [sqrt(s**2 + t**2) for (s, t) in zip(x, y)]

    if plots is True:
        plot(theta, r)
        title('Radius as a Function of Angle')
        show()

    p = [s*t for (s, t) in zip(r, theta)]

    q = elliptical_integral_numerical_confirmation(a, b, theta)

    if plots is True:
        plot(theta, r, c='0.85')
        plot(theta, p, 'k')
        plot(theta, q, 'b--')
        title('Normalized U-Substitution')
        legend([r'Radius: $r(\theta) = \sqrt{a^2sin^{2}\theta + b^2cos^{2}\theta}$',
                r'Normalized U-Substitution: $p(\theta) = \theta r(\theta)$',
                r'Numerical Confirmation: $p(\theta)$'])
        xlabel(r'$\theta_{0}^{2\pi}$')
        ylabel(r'$f(\theta)$')

        show()
    return a, b, p, theta

def elliptical_integral_numerical_confirmation(a, b, theta):
    q = []
    for s in theta:
        I = quad(lambda t: np.sqrt((b * np.cos(s)) ** 2 + (a * np.sin(s)) ** 2), 0, s)
        I = I[0]
        q.append(I)
    return q

def normalized_gaussian(n=n):
    std = pi
    theta = lin_space_zero(-std, std, n-1)
    gauss = [(exp(-1/2*(std**2)) - (exp(-1/2*(s**2))))/(exp(-1/2*(std)**2)-1) for s in theta]
    plot(theta, gauss, 'k')
    title('Normalized Gaussian')
    xlabel(r'Standard Deviation ($\sigma$)')
    ylabel(r'Normalized ($f_{n}(\theta)$)')
    legend([r'$\frac{e^{-\sigma^2/2} - e^{-\theta^2/2}}{e^{-\sigma^2/2} - 1}$'])
    show()

def elliptical_coordinate_system(a=2, b=1, n=n):
    from math import pi, sin, cos, sqrt
    theta = lin_zero(2*pi, n)
    x_coord = []
    y_coord = []
    r_coord = []
    for t in theta:
        x = [cos(t)/pi * a * s for s in theta]
        y = [sin(t)/pi * b * s for s in theta]
        r = [sqrt((u**2 + v**2)) for (u, v) in zip(x, y)]
        x_coord.append(x)
        y_coord.append(y)
        r_coord.append(r)
    plot(x_coord, y_coord, 'k')
    title('Elliptical Domain')
    xlabel('a')
    ylabel('b')
    show()
    fig = figure()
    ax = Axes3D(fig)
    ax.set_title('Elliptical Domain')
    ax.set_xlabel('a')
    ax.set_ylabel('b')
    ax.set_zlabel('r')

    ax.plot_surface(x_coord, y_coord, r_coord, cmap='gray')
    show()
    return theta, x_coord, y_coord, r_coord

def modes(modes=12):
    cos_summation = []
    theta = lin_space_zero(-pi/2, pi/2, 200)
    for m in range(1, modes):
        print('m' + str(m))
        wave = [cos(m * t) for t in theta]
        plot(theta, wave,  c='0.85')
        cos_summation.append(wave)
    cos_summation = [sum(s) for s in zip(*cos_summation)]
    plot(theta, cos_summation, 'k')
    plot(theta, [1]*len(theta), 'k')
    ylabel('Un-Normalized Intensity (#m)')
    xlabel('Forced (Constrained) Wavelength Span ($\pi$)')
    title('Idealistic Finite Mode (m=' + str(modes-1) + ') \n'
          'Single Slit Diffraction Propogation')
    show()


if __name__ == '__main__':
    the_unit_circle()
    normalized_gaussian()
    elliptical_coordinate_system()
    modes()



